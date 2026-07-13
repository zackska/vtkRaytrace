#include "vtkRaytrace.h"
#include <vtkGenericCell.h>
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

//Mesh input file reader
void vtkRaytrace::readMesh(char *inputFilename, int nContours, bool sharp, double isoValue)
{


	if(!volume) // if the data is composed of volumes with uniform refractive index
	{
		// Read in STL mesh using vtkSTLReader
		vtkSTLReader *reader = vtkSTLReader::New();
		reader->SetFileName(inputFilename);
		reader->Update();
		mesh->DeepCopy(reader->GetOutput());
		reader->Delete();

		if(!calcIndex)
		{
			vtkDoubleArray* index =  vtkDoubleArray::New();
			index->SetNumberOfComponents(1);
			index->SetName("index ratio");
			index->SetNumberOfTuples(mesh->GetNumberOfCells());
			index->FillComponent(0, glassIndex);
			// The data is added to FIELD data (rather than POINT data as usual)
			mesh->GetCellData()->AddArray(index);

		}
	}

	if(volume) // if the data is a field of non-uniform refractive index
	{
		vtkStructuredGrid *grid = vtkStructuredGrid::New();

		if(subset) // read every XML file in the directory specifice by inputFilename and combine the ones within a certain subspace
		{
			double distance = 40;
			meshTools::extractMultiblockSubset(inputFilename, distance, grid, "subset.vts");
		}
		else // just read the XML file, inputFilename is just a filename
		{
			vtkXMLStructuredGridReader *reader = vtkXMLStructuredGridReader::New();
			reader->SetFileName(inputFilename);
			reader->Update();
			grid->DeepCopy(reader->GetOutput());
			double *bounds = new double [6];
			grid->GetBounds(bounds);
			meshTools::printVec(bounds, "Bounds", 6);
		}

		// convert the structured grid data into polydata surface contours with refractive index ratios assigned.
		// sharp -> single physical interface (alpha=0.5) at the full liquid index (glassIndex);
		// diffuse -> original nested iso-index shells (for continuous fields, e.g. gas density).
		meshTools::gridScalarContours(grid, mesh, nContours, sharp, glassIndex, isoValue);
	}

	//Populate mesh normals for access in loop
	vtkPolyDataNormals *normal =
	vtkPolyDataNormals::New();

	normal->SetInputData(mesh);
	normal->ComputePointNormalsOn();   // smoothed (vertex-averaged) normals -> less facet scatter
	normal->ComputeCellNormalsOn();
	normal->SplittingOff();            // no feature splitting => fully smoothed across the surface
	normal->ConsistencyOn();
	normal->FlipNormalsOff();
	normal->AutoOrientNormalsOn();
	normal->Update();

	mesh->DeepCopy(normal->GetOutput());
	normal->Delete();

	// Cell locator for ray/mesh intersection. vtkStaticCellLocator is immutable after
	// BuildLocator() and its IntersectWithLine(...,vtkGenericCell*) is thread-safe — unlike
	// vtkOBBTree, whose traversal keeps mutable internal state and segfaults under OpenMP.
	// MUST be built AFTER the normals/DeepCopy: VTK 9.2 vtkContourFilter emits triangle
	// STRIPS, which vtkPolyDataNormals de-strips into individual triangles (~2x the cell
	// count). Building the locator on the pre-normals mesh left its cell ids referencing the
	// strip topology while trace() read the de-stripped mesh -> 2x array overrun / SIGSEGV.
	meshOBBTree = vtkStaticCellLocator::New();
	meshOBBTree->SetDataSet(mesh);
	meshOBBTree->CacheCellBoundsOn();
	meshOBBTree->BuildLocator();
	//delete[] stlName;

	// Cache mesh metrics now (serial) — GetLength()/GetCenter() write shared member
	// buffers and must not be called concurrently from the parallel trace() loop.
	m_meshLength = mesh->GetLength();
	double *cc = mesh->GetCenter();
	m_meshCenter[0]=cc[0]; m_meshCenter[1]=cc[1]; m_meshCenter[2]=cc[2];
}

void vtkRaytrace::plotNormals()
{
// Create a 'dummy' 'vtkCellCenters' to force the glyphs to the cell-centers
vtkCellCenters* centers = vtkCellCenters::New();
centers->VertexCellsOn();
centers->SetInputData(mesh);
centers->Update();

// Create a new 'default' arrow to use as a glyph
vtkArrowSource* arrow = vtkArrowSource::New();
arrow->Update();

vtkGlyph3D* glyph = vtkGlyph3D::New();
// Set its 'input' as the cell-center normals calculated at the sun's cells
glyph->SetInputData(centers->GetOutput());
// Set its 'source', i.e., the glyph object, as the 'arrow'
glyph->SetSourceData(arrow->GetOutput());
// Enforce usage of normals for orientation
glyph->SetVectorModeToUseNormal();
// Set scale for the arrow object
glyph->SetScaleFactor(1);
glyph->Update();
vtkPolyData*glyphMesh = glyph->GetOutput();

// Create an actor for the arrow-glyphs
meshTools::show(glyphMesh);
}




// ray auxiliary functions

void vtkRaytrace::addRay(vtkRenderer *renderer, double *p1, double *p2)
{
	vtkPolyData * linesPolyData = vtkPolyData::New();

	vtkPoints *pts = vtkPoints::New();
	pts->InsertNextPoint(p1);
	pts->InsertNextPoint(p2);

	linesPolyData->SetPoints(pts);

	vtkLine *line = vtkLine::New();
    line->GetPointIds()->SetId(0, 0); // the second 0 is the index of the Origin in linesPolyData's points
    line->GetPointIds()->SetId(1, 1);

	vtkCellArray *lines = vtkCellArray::New();
	lines->InsertNextCell(line);

	linesPolyData->SetLines(lines);

	unsigned char blue[3] = { 0, 0, 255 };

	vtkUnsignedCharArray *colors = vtkUnsignedCharArray::New();
	colors->SetNumberOfComponents(3);
	colors->InsertNextTypedTuple(blue);  // InsertNextTupleValue removed in VTK 9

	linesPolyData->GetCellData()->SetScalars(colors);

    vtkPolyDataMapper *mapper =
    vtkPolyDataMapper::New();
	mapper->SetInputData(linesPolyData);

    vtkActor *actor =
    vtkActor::New();
    actor->SetMapper(mapper);

    renderer->AddActor(actor);

	actor->Delete();
	mapper->Delete();
	colors->Delete();
	linesPolyData->Delete();
	line->Delete();
	lines->Delete();
	pts->Delete();
}

void vtkRaytrace::refract(double * norm, double * ray, double * newray, double indexRatio)
{
		vtkMath::Normalize(ray);

		//Calculate new refracted ray
		double refindex;
		double negnorm[3];   // was `new double` (1 elem) but used as a 3-vector → heap overflow + leak
		std::memcpy(negnorm, norm, sizeof(double)*3);
		vtkMath::MultiplyScalar(negnorm,-1);
		////cout << negnorm [0] << " " << negnorm[1] << " " << negnorm[2] << endl;
		if(vtkMath::AngleBetweenVectors(negnorm, ray)>vtkMath::Pi()/2)// normal doesn't point toward ray
		{
			std::memcpy(norm, negnorm, sizeof(double)*3); // reverse direction
			refindex = indexRatio; // means a hit on the interior
		}
		else
			refindex = 1/indexRatio;


		// perform refraction calculation
		if(refindex>1)
		{

			std::memcpy(negnorm, norm, sizeof(double)*3);
			vtkMath::MultiplyScalar(negnorm,-1);
			if(vtkMath::AngleBetweenVectors(negnorm, ray)>asin(1/refindex)) // total internal reflection
			{
				////cout << "TIR" << endl;
				double *reflect = new double [3];
				std::memcpy(reflect, norm, sizeof(double)*3);
				vtkMath::MultiplyScalar(reflect, 2*vtkMath::Dot(ray, norm));
				vtkMath::Subtract(ray, reflect, newray);
				delete[] reflect;
			}
			else
			{
				double *negnorm = new double [3];
				std::memcpy(negnorm, norm, sizeof(double)*3);
				vtkMath::MultiplyScalar(negnorm, -1);

				double *crossray = new double [3];
				vtkMath::Cross(norm, ray, crossray);

				double *crossnegray = new double [3];
				vtkMath::Cross(negnorm, ray, crossnegray);

				double *doublecross = new double [3];
				vtkMath::Cross(norm, crossnegray, doublecross);

				vtkMath::MultiplyScalar(doublecross, refindex);

				double *finalvec = new double [3];
				std::memcpy(finalvec, norm, sizeof(double)*3);
				vtkMath::MultiplyScalar(finalvec, sqrt(1-pow(refindex, 2)*vtkMath::Dot(crossray, crossray)));

				vtkMath::Subtract(doublecross, finalvec, newray);

				delete[] doublecross;
				delete[] finalvec;
				delete[] crossnegray;
				delete[] crossray;
				delete[] negnorm;
			}
		}
		else
		{
			double *negnorm = new double  [3];
			std::memcpy(negnorm, norm, sizeof(double)*3);
			vtkMath::MultiplyScalar(negnorm, -1);

			double *crossray = new double  [3];
			vtkMath::Cross(norm, ray, crossray);

			double *crossnegray  = new double  [3];
			vtkMath::Cross(negnorm, ray, crossnegray);


			double *doublecross  = new double  [3];
			vtkMath::Cross(norm, crossnegray, doublecross);

			vtkMath::MultiplyScalar(doublecross, refindex);

			double *finalvec = new double [3];
			std::memcpy(finalvec, norm, sizeof(double)*3);

			double factor = sqrt(1-pow(refindex, 2)*vtkMath::Dot(crossray, crossray));
			vtkMath::MultiplyScalar(finalvec, factor);

			vtkMath::Subtract(doublecross, finalvec, newray);

			delete[] doublecross;
			delete[] finalvec;
			delete[] crossnegray;
			delete[] crossray;
			delete[] negnorm;
		}
}

void vtkRaytrace::getRay(double *p1, double *p2, double * ray)
{
	vtkMath::Subtract(p2, p1, ray);
	vtkMath::Normalize(ray);
}

void vtkRaytrace::writeBMP()
{
FILE *f;
int w = nPixels[0];
int h = nPixels[1];
unsigned char *img = NULL;
int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int
if( img )
    free( img );
img = (unsigned char *)malloc(3*w*h);
memset(img,0,sizeof(img));

for(int i=0; i<w; i++)
{
    for(int j=0; j<h; j++)
{
    int x=i; int y=(h-1)-j;
	int intensity = image[i][j]*255;
    int r = intensity;
    int g = intensity;
    int b = intensity;
    if (r > 255) r=255;
    if (g > 255) g=255;
    if (b > 255) b=255;
    img[(x+y*w)*3+2] = (unsigned char)(r);
    img[(x+y*w)*3+1] = (unsigned char)(g);
    img[(x+y*w)*3+0] = (unsigned char)(b);
}
}

unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
unsigned char bmppad[3] = {0,0,0};

bmpfileheader[ 2] = (unsigned char)(filesize    );
bmpfileheader[ 3] = (unsigned char)(filesize>> 8);
bmpfileheader[ 4] = (unsigned char)(filesize>>16);
bmpfileheader[ 5] = (unsigned char)(filesize>>24);

bmpinfoheader[ 4] = (unsigned char)(       w    );
bmpinfoheader[ 5] = (unsigned char)(       w>> 8);
bmpinfoheader[ 6] = (unsigned char)(       w>>16);
bmpinfoheader[ 7] = (unsigned char)(       w>>24);
bmpinfoheader[ 8] = (unsigned char)(       h    );
bmpinfoheader[ 9] = (unsigned char)(       h>> 8);
bmpinfoheader[10] = (unsigned char)(       h>>16);
bmpinfoheader[11] = (unsigned char)(       h>>24);

const char *outName = getenv("VTKRT_OUT"); // lets concurrent processes write distinct files
f = fopen(outName ? outName : "img.bmp", "wb");
fwrite(bmpfileheader,1,14,f);
fwrite(bmpinfoheader,1,40,f);
for(int i=0; i<h; i++)
{
    fwrite(img+(w*(h-i-1)*3),3,w,f);
    fwrite(bmppad,1,(4-(w*3)%4)%4,f);
}
fclose(f);
}

void vtkRaytrace::makeCamera(bool confine)
{
	double *cameraZ = new double [3];
	double *cameraX = new double [3];
	for(int i=0; i<3;i++)
	{
		cameraZ[i] = 1;
		cameraX[i] = 0;
	}

	vtkMath::Subtract(cameraZ, opticalAxis, cameraZ);

	bool found = 1;
	int i = 0;
	while(found)
	{
		if(cameraZ[i])
		{
				cameraX[i] = 1;
				cameraZ[i] = 0;
				found = !found;
		}
		else
			i++;
	}


	// Framing reference: per-mesh bounds by default, or a FIXED box via env
	// VTKRT_FRAME_BOUNDS="xmin xmax ymin ymax zmin zmax" so an image SEQUENCE shares
	// one camera (no inter-frame jitter / zoom as the jet length changes).
	double frameBounds[6];
	double *mb = mesh->GetBounds();
	for(int i=0;i<6;i++) frameBounds[i] = mb[i];
	const char *fb = getenv("VTKRT_FRAME_BOUNDS");
	if(fb)
		sscanf(fb, "%lf %lf %lf %lf %lf %lf",
		       &frameBounds[0],&frameBounds[1],&frameBounds[2],
		       &frameBounds[3],&frameBounds[4],&frameBounds[5]);
	double frameCenter[3] = { 0.5*(frameBounds[0]+frameBounds[1]),
	                          0.5*(frameBounds[2]+frameBounds[3]),
	                          0.5*(frameBounds[4]+frameBounds[5]) };

	 // define a 3-dimensional point for the centroid of the canvas
	cameraPos = new double [3];
	std::memcpy(cameraPos, opticalAxis, sizeof(double)*3);
	vtkMath::MultiplyScalar(cameraPos, focalPlane + focalLength);
	vtkMath::Add(cameraPos, frameCenter, cameraPos); // add centroid coordinates along directions orthogonal to optical axis

	//calculate perspective point from camera parameters
	perspectivePoint = new double [3];
	std::memcpy(perspectivePoint, opticalAxis, sizeof(double)*3);
	vtkMath::MultiplyScalar(perspectivePoint, magnification*focalLength/(1-magnification));
	vtkMath::Add(perspectivePoint, cameraPos, perspectivePoint);
	//printVec(perspectivePoint, "Perspective Point", 3);

	// Thin-lens objective: cache the in-plane axes, the focus-plane centre (jet axis),
	// and the aperture radius implied by the f-number. aperture = focusDist/(2*fNumber);
	// a large fNumber -> ~0 aperture -> the original pinhole camera.
	std::memcpy(m_cameraX, cameraX, 3*sizeof(double));
	std::memcpy(m_cameraZ, cameraZ, 3*sizeof(double));
	std::memcpy(m_focusCenter, frameCenter, 3*sizeof(double));
	double pdiff[3];
	vtkMath::Subtract(perspectivePoint, frameCenter, pdiff);
	double focusDistAxial = fabs(vtkMath::Dot(pdiff, opticalAxis));
	m_apertureRadius = (m_fNumber > 1e-6) ? focusDistAxial/(2.0*m_fNumber) : 0.0;


	// pixelSize = physical length per pixel. The original used mesh->GetMaxCellSize(),
	// which returns the max POINTS-per-cell (==3 for triangles), NOT a length — at SI
	// (meter) scale that floored the pixel count to 0 (blank image). Instead we frame
	// the mesh: square pixels sized so the larger transverse extent spans the requested
	// resolution, perpendicular to the optical axis.
	double *bounds = frameBounds; // per-mesh bounds, or the fixed sequence frame (see above)
	double extent[3];
	for(int i=0;i<3;i++)
		extent[i] = bounds[2*i+1] - bounds[2*i]; // frame size along each axis (>=0)

	double extentX = abs(vtkMath::Dot(cameraX, extent)); // transverse extents (perp to optical axis)
	double extentZ = abs(vtkMath::Dot(cameraZ, extent));
	double maxT    = (extentX > extentZ) ? extentX : extentZ;
	int    reqMax  = (nPixels[0] > nPixels[1]) ? nPixels[0] : nPixels[1];

	double pixelSize = maxT / reqMax; // square pixel; long transverse axis = reqMax pixels

	if(confine) // shrink the chip to the mesh's silhouette (no wasted border)
	{
		nPixels[0] = (int)vtkMath::Floor(extentX/pixelSize); if(nPixels[0]<1) nPixels[0]=1;
		nPixels[1] = (int)vtkMath::Floor(extentZ/pixelSize); if(nPixels[1]<1) nPixels[1]=1;
	}

	//Set chip size (physical) = nPixels * pixelSize  ->  spans the mesh extent
	chipSize = new double [2];
	chipSize[0] = double(nPixels[0]);
	chipSize[1] = double(nPixels[1]);
	vtkMath::MultiplyScalar(chipSize, pixelSize);
	//printVec(chipSize, "Chip size", 2);

	// Define light source centroid
	lightCenter = new double[3];
	std::memcpy(lightCenter, opticalAxis, sizeof(double)*3);
	vtkMath::MultiplyScalar(lightCenter, lightLoc);
	vtkMath::Add(lightCenter,  frameCenter, lightCenter);


	world2Camera = vtkMatrix4x4::New();
	world2Camera->Identity();

	for(int i=0; i<3;i++)
	{
		world2Camera->SetElement(0, i, cameraX[i]);
		world2Camera->SetElement(1, i, opticalAxis[i]);
		world2Camera->SetElement(2, i, cameraZ[i]);
		world2Camera->SetElement(3, i, cameraPos[i]);
	}

	//vtkIndent *indent = vtkIndent::New();
	//world2Camera->PrintSelf(cout, *indent);
	world2Camera->Transpose();// transpose the matrix since this is right-hand multiplication
}

//ray tracing calculation
double vtkRaytrace::trace(double *source, double *target, int depth, vtkRenderer *&renderer, bool insideMedium)
{
	if(depth<depthLimit)
	{
		// Thread-safe nearest intersection along source->target. The vtkPoints/vtkIdList
		// overload, GetTuple(double*) and GetCenter()/GetLength() all share internal
		// buffers and corrupt/crash under multithreading; the vtkGenericCell overload
		// gives each thread its own scratch cell.
		// One scratch cell PER THREAD (not per ray): constructing/destroying VTK objects
		// in the hot loop corrupts VTK's global object bookkeeping across threads.
		static thread_local vtkSmartPointer<vtkGenericCell> gcell = vtkSmartPointer<vtkGenericCell>::New();
		double tline = 0, xhit[3], pcoords[3];
		int subId = 0;
		vtkIdType cellId = -1;
		int hit = meshOBBTree->IntersectWithLine(source, target, 1e-6, tline, xhit, pcoords, subId, cellId, gcell);

		double intensity = 0;
		double point[3];
		double ray[3];                 // local scratch (was the class member -> leak + non-reentrant)
		getRay(source, target, ray);

		if(hit)// if the ray intersects the mesh
		{

			point[0]=xhit[0]; point[1]=xhit[1]; point[2]=xhit[2];
			//printVec(point, "Intersection");
			if(plot)
				addRay(renderer, source, point);

			// --- Beer-Lambert: attenuate over the segment just traversed, but only if it
			//     lay inside the liquid (ray is still unit here from getRay). ---
			double segLen = sqrt(vtkMath::Distance2BetweenPoints(source, point));
			double atten  = insideMedium ? exp(-absorption*segLen) : 1.0;

			// Smooth (Phong-interpolated) surface normal: blend the hit cell's point
			// normals by the parametric coords of the intersection, instead of the faceted
			// per-cell normal. Far less spurious scatter -> more realistic transmission.
			double norm[3] = {0.0, 0.0, 0.0};
			double wts[8];
			gcell->InterpolateFunctions(pcoords, wts);
			vtkDataArray *pnorm = mesh->GetPointData()->GetNormals();
			int npc = gcell->GetNumberOfPoints();
			for(int kk=0; kk<npc; kk++)
			{
				double n3[3];
				pnorm->GetTuple(gcell->GetPointId(kk), n3); // thread-safe copy overload
				norm[0]+=wts[kk]*n3[0]; norm[1]+=wts[kk]*n3[1]; norm[2]+=wts[kk]*n3[2];
			}
			vtkMath::Normalize(norm);
			double nU[3] = { norm[0], norm[1], norm[2] };
			vtkMath::Normalize(nU);
			double cosI = fabs(vtkMath::Dot(ray, nU));
			if(cosI > 1.0) cosI = 1.0;
			double sinI = sqrt(1.0 - cosI*cosI);

			// --- Fresnel transmittance (unpolarised) across the interface. n1->n2 depends
			//     on whether we are entering or leaving the liquid. ---
			double n1 = insideMedium ? glassIndex : 1.0;
			double n2 = insideMedium ? 1.0 : glassIndex;
			double sinT = (n1/n2)*sinI;
			double T;
			if(sinT >= 1.0)
				T = 0.0;                                   // total internal reflection
			else
			{
				double cosT = sqrt(1.0 - sinT*sinT);
				double rs = (n1*cosI - n2*cosT)/(n1*cosI + n2*cosT);
				double rp = (n1*cosT - n2*cosI)/(n1*cosT + n2*cosI);
				T = 1.0 - 0.5*(rs*rs + rp*rp);
			}

			if(T <= 0.0)                                   // TIR: transmitted path is extinguished
				return 0.0;

			// refract ray. GetComponent is thread-safe; the double* GetTuple shares an
			// internal per-array buffer and races across threads.
			double newray[3];
			double idx = mesh->GetCellData()->GetArray("index ratio")->GetComponent(cellId, 0);
			refract(norm, ray, newray, idx);

			//calculate new target point: step a small, MESH-RELATIVE epsilon past the
			//interface (the old fixed 0.1 was absolute and jumped clean off mm-scale meshes).
			double eps = 1e-3 * m_meshLength;
			vtkMath::MultiplyScalar(newray, eps);
			vtkMath::Add(point, newray, point);
			vtkMath::MultiplyScalar(newray, 1e5);
			vtkMath::Add(point, newray, target);

			//increment recursion depth; crossing the interface flips the medium
			depth++;
			intensity = atten * T * trace(point, target, depth, renderer, !insideMedium);
		}

		else //if the ray missed the mesh, test if it hits the light source. If not, there is no light hitting the source point, so return 0 intensity.
		{
			vtkMath::MultiplyScalar(ray, abs(vtkMath::Dot(source, opticalAxis)-(lightLoc+vtkMath::Dot(m_meshCenter, opticalAxis)))); // trace the ray to the plane of the light source
			vtkMath::Add(source, ray, point);

			// Direction from the scene toward the back-light. The ray travels camera->light,
			// i.e. roughly along sign(lightLoc)*opticalAxis (here lightLoc<0 => -opticalAxis).
			// The original compared against +opticalAxis, so unobstructed background rays came
			// out at ~180deg and were rejected -> inverted (dark-background) shadowgram.
			double lightDir[3] = { opticalAxis[0], opticalAxis[1], opticalAxis[2] };
			if(lightLoc < 0) vtkMath::MultiplyScalar(lightDir, -1.0);

			if(sqrt(vtkMath::Distance2BetweenPoints(point, lightCenter))<lightRadius && abs(vtkMath::AngleBetweenVectors(ray, lightDir))<acceptanceAngle)
			{
				intensity = 1; // it hit the light, return 1 for intensity
				if(plot)
					addRay(renderer, source, point);
			}
			else
				intensity = 0; // missed the light, return 0
		}

		return intensity; // gcell auto-frees (vtkSmartPointer)
	}

	else // if recursion depth is exceeded, return 0
		return 0;

}

void vtkRaytrace::visibilityTrace()
{

	// Pre-allocate the image rows (serial, cheap).
	image = new double*[nPixels[0]];
	for(int i = 0; i<nPixels[0]; i++)
		image[i] = new double[nPixels[1]];

	time_t start = std::time(0);// start timer
	int rowsDone = 0;

	// The pixel grid is embarrassingly parallel — distribute rows over threads with
	// OpenMP. Each thread owns its RNG and stack scratch; the OBBTree is read-only.
	// Requires plot==0 (drawing debug rays into a shared renderer is not thread-safe).
	#pragma omp parallel
	{
#ifdef _OPENMP
		std::mt19937 rng(0x9E3779B9u ^ (unsigned)(2654435761u * omp_get_thread_num()));
#else
		std::mt19937 rng(0x9E3779B9u);
#endif
		std::uniform_real_distribution<double> uni(0.0, 1.0);
		vtkRenderer* renderer = NULL; // unused while plot==0

		#pragma omp for schedule(dynamic)
		for(int i = 0; i<nPixels[0]; i++)
		{
			for(int j = 0; j<nPixels[1]; j++)
			{
				double intensity = 0.0;
				std::vector<double> intensities;
				intensities.reserve(initialSample);

				int k = 0;
				double mean = 0.0, variance = 0.0, stdErr = 1e30;
				double target[4], newTarget[3], ray[3], tc[4];

				// Sample until we have >= initialSample rays AND the standard error of
				// the mean drops below varianceLimit (the correct MC convergence test;
				// the old code stopped on |mean-oldMean| and never used the variance).
				while(k < initialSample || stdErr > varianceLimit)
				{
					// jittered sub-pixel coordinate (antialiasing)
					target[0] = chipSize[0]*((i+uni(rng))/nPixels[0] - .5);
					target[1] = 0;
					target[2] = chipSize[1]*((j+uni(rng))/nPixels[1] - .5);
					target[3] = 1;

					world2Camera->MultiplyPoint(target, tc); // pixel position in world
					std::memcpy(newTarget, tc, sizeof(double)*3);
					getRay(perspectivePoint, newTarget, ray); // chief-ray direction through pixel

					double origin[3], fwd[3];
					if(m_apertureRadius > 0.0)
					{
						// THIN-LENS objective: each sample is a ray from a random point on the
						// finite aperture, through the pixel's focus-plane conjugate. Integrating
						// over the aperture = finite collection NA (more transillumination + DOF).
						double pf[3];
						vtkMath::Subtract(m_focusCenter, perspectivePoint, pf);
						double denom = vtkMath::Dot(ray, opticalAxis);
						double sFocus = (fabs(denom) > 1e-12) ? vtkMath::Dot(pf, opticalAxis)/denom : 1.0;
						double F[3] = { perspectivePoint[0]+sFocus*ray[0],
						                perspectivePoint[1]+sFocus*ray[1],
						                perspectivePoint[2]+sFocus*ray[2] };
						double rr = m_apertureRadius*sqrt(uni(rng));     // concentric disk sample
						double ph = 2.0*vtkMath::Pi()*uni(rng);
						for(int c=0;c<3;c++)
							origin[c] = perspectivePoint[c] + rr*cos(ph)*m_cameraX[c] + rr*sin(ph)*m_cameraZ[c];
						getRay(origin, F, fwd);                          // lens point -> focus point
					}
					else
					{
						// pinhole fallback (m_fNumber huge)
						origin[0]=perspectivePoint[0]; origin[1]=perspectivePoint[1]; origin[2]=perspectivePoint[2];
						fwd[0]=ray[0]; fwd[1]=ray[1]; fwd[2]=ray[2];
					}
					vtkMath::Normalize(fwd);
					double farPt[3] = { origin[0]+1000.0*fwd[0], origin[1]+1000.0*fwd[1], origin[2]+1000.0*fwd[2] };

					double s = trace(origin, farPt, 1, renderer);
					intensities.push_back(s);
					intensity += s;
					k++;
					mean = intensity/k;

					variance = 0.0;
					for(int l = 0; l<k; l++)
						variance += (intensities[l]-mean)*(intensities[l]-mean);
					variance /= k;
					stdErr = (k>1) ? sqrt(variance/k) : 1e30;
				}

				image[i][j] = mean;
			}

			// coarse progress: atomic row counter, printed ~every 10% of rows
			int done;
			#pragma omp atomic capture
			done = ++rowsDone;
			int step = (nPixels[0] >= 10) ? nPixels[0]/10 : 1;
			if(done % step == 0)
			{
				double secs = difftime(std::time(0), start);
				#pragma omp critical
				cout << (100*done/nPixels[0]) << "% (" << done << "/" << nPixels[0]
				     << " rows) at " << secs << " s" << endl;
			}
		}
	}

	// NOTE: debug ray visualisation (meshTools::show) was removed here — it is
	// incompatible with the parallel loop (no shared renderer accumulates rays).
}




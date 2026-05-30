#include "vtkRaytrace.h"

//Mesh input file reader
void vtkRaytrace::readMesh(char *inputFilename, int nContours)
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
		
		// convert the structured grid data into polydata surface contours with refractive index ratios assigned
		meshTools::gridScalarContours(grid, mesh, nContours);
	}

	//OBBTree objects for chip and mesh for intersection calculations
	meshOBBTree = vtkOBBTree::New();
	meshOBBTree->SetDataSet(mesh);
	meshOBBTree->CacheCellBoundsOn();
	meshOBBTree->LazyEvaluationOff();
	meshOBBTree->UseExistingSearchStructureOn();
	meshOBBTree->SetNumberOfCellsPerNode(16);
	meshOBBTree->RetainCellListsOn();
	meshOBBTree->BuildLocator();
	
	//Populate mesh normals for access in loop
	vtkPolyDataNormals *normal = 
	vtkPolyDataNormals::New();
	
	normal->SetInputData(mesh);
	normal->ComputePointNormalsOff();
	normal->ComputeCellNormalsOn();
	normal->SplittingOff();
	normal->FlipNormalsOff();
	normal->AutoOrientNormalsOn();
	normal->Update();

	mesh->DeepCopy(normal->GetOutput());
	normal->Delete();
	//delete[] stlName;
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
	colors->InsertNextTupleValue(blue);

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
		double *negnorm = new double;
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

f = fopen("img.bmp","wb");
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


	 // define a 3-dimensional point for the centroid of the canvas
	cameraPos = new double [3];
	std::memcpy(cameraPos, opticalAxis, sizeof(double)*3);
	vtkMath::MultiplyScalar(cameraPos, focalPlane + focalLength);
	vtkMath::Add(cameraPos, mesh->GetCenter(), cameraPos); // add centroid coordinates along directions orthogonal to optical axis

	//calculate perspective point from camera parameters
	perspectivePoint = new double [3];
	std::memcpy(perspectivePoint, opticalAxis, sizeof(double)*3);
	vtkMath::MultiplyScalar(perspectivePoint, magnification*focalLength/(1-magnification));
	vtkMath::Add(perspectivePoint, cameraPos, perspectivePoint);
	//printVec(perspectivePoint, "Perspective Point", 3);
	
	
	double pixelSize = mesh->GetMaxCellSize()*cellSizeFactor*magnification; // calculate the pixelSize from the mesh's cell size

	if(confine) // minimize the number of pixels so that the chip is not larger than the mesh's image
	{
	double *bounds = mesh->GetBounds();
	double *extent = new double [3];

	for(int i =0;i<3;i++)
	{
		extent[i] = bounds[2*i] - bounds[2*i+1]; // calculate the mesh's size
	}

	double extentX =  abs(vtkMath::Dot(cameraX, extent));
	double extentZ =  abs(vtkMath::Dot(cameraZ, extent));

	magnification = .9; // define magnification of the scene on the image

	// Limit the number of pixels to the extent of the mesh's image on the chip divided by the decided pixel size
	if(nPixels[0] > vtkMath::Floor(extentX*magnification/pixelSize))
		nPixels[0] = vtkMath::Floor(extentX*magnification/pixelSize);
	if(nPixels[1] > vtkMath::Floor(extentZ*magnification/pixelSize))
		nPixels[1] = vtkMath::Floor(extentZ*magnification/pixelSize);

	delete[] extent;
	}

	//Set chip size
	chipSize = new double [2];
	chipSize[0] = double(nPixels[0]);
	chipSize[1] = double(nPixels[1]);
	vtkMath::MultiplyScalar(chipSize, pixelSize); 
	//printVec(chipSize, "Chip size", 2);

	// Define light source centroid
	lightCenter = new double[3];
	std::memcpy(lightCenter, opticalAxis, sizeof(double)*3);
	vtkMath::MultiplyScalar(lightCenter, lightLoc);
	vtkMath::Add(lightCenter,  mesh->GetCenter(), lightCenter);

	
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
double vtkRaytrace::trace(double *source, double *target, int depth, vtkRenderer *&renderer)
{
	if(depth<depthLimit)
	{
		// test for mesh intersections along the ray defined by source and target
		vtkPoints* x = vtkPoints::New();
		vtkIdList* cellID = vtkIdList::New();
		int hit = meshOBBTree->IntersectWithLine(source, target, x, cellID); 

		double intensity;
		double *point = new double[3];
		ray = new double [3];
		getRay(source, target, ray);

		if(hit)// if the ray intersects the mesh
		{

			x->GetPoint(0, point);
			//printVec(point, "Intersection");
			if(plot)
				addRay(renderer, source, point);

			// get cell normal at intersection point
			double *norm  = new double [3];
			meshTools::getCellNormal(mesh, cellID->GetId(0), norm);

			// refract ray
			double *newray  = new double [3];
			double *indexRatio = new double [1];
			indexRatio = mesh->GetCellData()->GetArray("index ratio")->GetTuple(cellID->GetId(0)); // get indexRatio at intersection cell
			refract(norm, ray, newray, indexRatio[0]);
			
			//calculate new target point
			vtkMath::MultiplyScalar(newray, .1);
			vtkMath::Add(point, newray, point);
			vtkMath::MultiplyScalar(newray, 100000);
			vtkMath::Add(point, newray, target);

			//increment recursion depth
			depth++;

			// recurse
			intensity = trace(point, target, depth, renderer);
			
			//free memory
			delete[] norm;
			delete[] newray;
	
		}

		else //if the ray missed the mesh, test if it hits the light source. If not, there is no light hitting the source point, so return 0 intensity.
		{
			vtkMath::MultiplyScalar(ray, abs(vtkMath::Dot(source, opticalAxis)-(lightLoc+vtkMath::Dot(mesh->GetCenter(), opticalAxis)))); // trace the ray to the plane of the light source
			vtkMath::Add(source, ray, point);
			
			if(sqrt(vtkMath::Distance2BetweenPoints(point, lightCenter))<lightRadius && abs(vtkMath::AngleBetweenVectors(ray, opticalAxis))<acceptanceAngle) 
			{
				intensity = 1; // it hit the light, return 1 for intensity
				cout << "HIT!" << endl;
				if(plot)
					addRay(renderer, source, point);
			}
			else
				intensity = 0; // missed the light, return 0
			
			
		}
	
		return intensity;

		//free memory
		cellID->Delete();
		x->Delete();
		delete[] ray;
		delete[] point;
	}
	
	else // if recursion depth is exceeded, return 0
		return 0;

}

void vtkRaytrace::visibilityTrace()
{
	
	// renderer for debugging
	vtkRenderer* renderer = vtkRenderer::New();

	int percentageCounter = 0;

	double *target = new double [4];
	double *ray = new double [3];
	double *newTarget = new double [3];
	
	//Calculate pixel intensities using ray tracing
	image = new double*[nPixels[0]];
	time_t start = std::time(0);// start timer
	srand((unsigned)time(NULL)); // seed random number generator
	for(int i = 0; i<nPixels[0];i++)
	{
		image[i] = new double[nPixels[1]];
		for(int j = 0; j<nPixels[1];j++)
		{
			double intensity = 0;
			std::vector<double> intensities(initialSample, 0);

			int k = 0;
			double mean;
			double oldMean = 0;
			double variance = 0;
			double tol = 0;

			while(k<initialSample || tol > varianceLimit)
			{
				if(k+1>initialSample) // reallocate intensity vector to increase size
				{
					intensities.resize(k+1);
				}	

				// define pixel coordinate
				target[0] = chipSize[0]*((i+((double)rand()/(double)RAND_MAX))/nPixels[0] - .5);
				target[1] = 0;
				target[2] = chipSize[1]*((j+((double)rand()/(double)RAND_MAX))/nPixels[1] - .5);

				target[3] = 1; //extend target point definition to allow translation during camera-to-world transformation
			
				// transform pixel coordinate to scene coordinate for initial ray calculation
				world2Camera->MultiplyPoint(target, target);
				std::memcpy(newTarget, target, sizeof(double)*3);
				getRay(perspectivePoint, newTarget, ray);
				vtkMath::MultiplyScalar(ray, 1000);
				vtkMath::Add(perspectivePoint, ray, newTarget);
				//printVec(perspectivePoint, "Perspective Point", 3);
			
				// determine intensity of this sample
				intensities[k] = trace(perspectivePoint, newTarget, 1, renderer);

				 

				intensity += intensities[k];

				mean = intensity/(k+1);

				//cout << "mean = " << mean << endl;

				if(k>0)
				tol = abs(mean - oldMean);

				oldMean = mean;


				// calculate variance
				variance = 0;
				for(int l = 0; l<k+1;l++)
				{
					variance += pow(intensities[l] - mean, 2)/(k+1);
				}
				
				if(tol!=0)
				{
				//cout << "tol = " << tol << endl;
				}
				k++; // increment iterator

			}

			image[i][j] = mean;

			// progress monitoring
			double percentage = 100*(i*nPixels[1]+j)/(nPixels[0]*nPixels[1]);
			int percent = vtkMath::Floor(percentage);
			if(percent%5==0 && percent!=0 && percent > percentageCounter*5)
			{
				percentageCounter++;
				time_t thisPercentage = std::time(0);
				double time = difftime(thisPercentage, start);

				if(time>3600)
				cout << percent << "%" << " at " << time/3600 << " hours." << endl;
				if(time<3600 && time>60)
				cout << percent << "%" << " at " << time/60 << " minutes." << endl;
				else
				cout << percent << "%" << " at " << time << " seconds." << endl;
			}
		}
	}


	if(plot)
		meshTools::show(renderer, mesh);
	
	
}




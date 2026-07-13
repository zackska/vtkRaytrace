#include <iostream>
#include <cstring>
#include <fstream>
#include "meshTools.h"
#include "vtkRaytrace.h"
using namespace std;

// Headless batch driver for vtkRaytrace.
// Usage: ReadSTL [input.stl] [nPixels] [output.bmp-handled-internally]
//   argv[1] : path to triangulated STL surface (default below)
//   argv[2] : square image resolution in pixels (default 400)

int main (int argc, char **argv)
{
	// This is a back-lit illumination ray tracing code designed to handle transmissive objects.
	vtkRaytrace r;                       // central object
    vtkRaytrace *rayTrace = &r;

	// Interface model:
	//   VTKRT_INTERFACE=sharp   -> one refracting surface at the alpha=0.5 iso-level (resolved
	//                              two-phase liquid/gas interface). Physically correct for VoF.
	//   VTKRT_INTERFACE=diffuse -> (default) original nested iso-index shells, for a continuous
	//                              index field (e.g. compressible gas density shadowgraphy).
	//   VTKRT_NCONTOURS=<n>     -> shell count in the DIFFUSE case (default 10).
	//   VTKRT_ISOVALUE=<v>      -> SHARP iso-level in rho_s units (default: midpoint of the field
	//                              range, i.e. alpha=0.5 for an [n_gas, n_liquid] field).
	int nContours = (getenv("VTKRT_NCONTOURS")) ? atoi(getenv("VTKRT_NCONTOURS")) : 10;
	const char *ifaceEnv = getenv("VTKRT_INTERFACE");
	bool sharp = (ifaceEnv && strcmp(ifaceEnv, "sharp") == 0);
	double isoValue = (getenv("VTKRT_ISOVALUE")) ? atof(getenv("VTKRT_ISOVALUE")) : -1.0;
	if(nContours < 3) nContours = 3;   // diffuse loop runs k=1..nContours-2; needs >=3

	// Read in mesh file and convert to density iso-contours, then calculate normals and OBBTree
	rayTrace->mesh = vtkPolyData::New();
	// VTKRT_VOLUME=1 -> read a single .vts structured grid (cell array "rho_s" = refractive
	// index n) and trace through nested iso-index shells. Default 0 = STL surface mode.
	rayTrace->volume = (getenv("VTKRT_VOLUME")) ? atoi(getenv("VTKRT_VOLUME")) : 0;
	rayTrace->subset = 0;   // single .vts file (not a directory of blocks)

	char *inputFilename = new char [512];
	if (argc > 1)
		strcpy(inputFilename, argv[1]);
	else
		strcpy(inputFilename,
			"/home/zfalg/OpenFOAM/zfalg-13/run/rayleighJet_sweep_freq/stl_out/jet_U0.5_f1.35x_A0.05_t0.750.stl");

	int res = (argc > 2) ? atoi(argv[2]) : 400;
	double absorptionCoeff = (argc > 3) ? atof(argv[3]) : 80.0; // 1/length; ~half-attenuation per jet diameter

    ifstream ifile(inputFilename);
    if(!ifile){ printf("File does not exist: %s\n", inputFilename); return 1; }
	ifile.close();

	// Liquid refractive index. Must be set BEFORE readMesh, which bakes the
	// "index ratio" cell array from glassIndex. VTKRT_GLASSINDEX -> ~1.0 makes the
	// jet index-matched (no refraction/TIR) so the image is pure Beer-Lambert
	// absorption -> graded transmission through the liquid.
	rayTrace->glassIndex = (getenv("VTKRT_GLASSINDEX")) ? atof(getenv("VTKRT_GLASSINDEX")) : 1.33;

	rayTrace->readMesh(inputFilename, nContours, sharp, isoValue);
	cout << "Loaded mesh: " << rayTrace->mesh->GetNumberOfCells() << " cells" << endl;
	// NOTE: meshTools::show() opens a blocking render window — disabled for headless batch.
	//meshTools::show(rayTrace->mesh);

	//Define camera parameters
	rayTrace->focalPlane = 0; // define image plane distance from mesh centroid along optical axis
	rayTrace->focalLength = 60; // define the focal length of the camera system
	rayTrace->magnification = .9; // define magnification of the scene on the image

	// Camera chip parameters
	rayTrace->nPixels = new int [2];
	rayTrace->nPixels[0] = res;
	rayTrace->nPixels[1] = res;
	rayTrace->cellSizeFactor = .5; // ratio between pixel size and the largest mesh cell's image

	// define the Optical axis for the setup. VTKRT_OPTAXIS="x y z" to view perpendicular to a
	// jet (e.g. "0 0 1" for a jet running along x -> side-profile shadowgram). Default +x.
	rayTrace->opticalAxis = new double [3];
	rayTrace->opticalAxis[0] = 1;
	rayTrace->opticalAxis[1] = 0;
	rayTrace->opticalAxis[2] = 0;
	if(getenv("VTKRT_OPTAXIS"))
		sscanf(getenv("VTKRT_OPTAXIS"), "%lf %lf %lf",
		       &rayTrace->opticalAxis[0], &rayTrace->opticalAxis[1], &rayTrace->opticalAxis[2]);

	// plotting rays and mesh for debugging
	rayTrace->plot = 0;

	// if appropriate, confine image to mesh bounds
	bool confineToBounds = 1;

	// Variance limit for Monte Carlo
	rayTrace->varianceLimit = .1;
	rayTrace->initialSample = 5; // number of rays to trace before evaluating variance

	// Light source parameters
	rayTrace->lightLoc = -50; // location of light along optical axis
	rayTrace->lightRadius = 200; // size of light source
	rayTrace->acceptanceAngle = 45*vtkMath::Pi()/180; // spreading angle of the light source in radians

	// Optional env overrides to sweep the ILLUMINATION numerical aperture / geometry:
	//   VTKRT_ACCEPT_DEG  - light acceptance half-angle in degrees (collimated <-> diffuse)
	//   VTKRT_LIGHTLOC    - light plane distance along optical axis (negative = behind jet)
	//   VTKRT_LIGHTRADIUS - light source disk radius
	{ const char *e;
	  if((e=getenv("VTKRT_ACCEPT_DEG")))  rayTrace->acceptanceAngle = atof(e)*vtkMath::Pi()/180.0;
	  if((e=getenv("VTKRT_LIGHTLOC")))    rayTrace->lightLoc        = atof(e);
	  if((e=getenv("VTKRT_LIGHTRADIUS"))) rayTrace->lightRadius     = atof(e); }

	// Objective: thin-lens collection aperture. f-number sets the collection NA
	// (large f-number -> pinhole). Default f/4, focus plane on the jet axis.
	rayTrace->m_fNumber = (getenv("VTKRT_FNUMBER")) ? atof(getenv("VTKRT_FNUMBER")) : 4.0;

	// Construct camera representation
	rayTrace->makeCamera(confineToBounds);

	// recursion limit for individual ray traces
	rayTrace->depthLimit = 7;

	// Light transport parameters: STL = uniform-index volume, so calcIndex MUST be 0
	// (calcIndex=1 skips creating the "index ratio" cell array -> null deref in trace()).
	rayTrace->calcIndex = 0;
	// glassIndex already set from VTKRT_GLASSINDEX (default 1.33) before readMesh.
	rayTrace->K_GD = .00023; // gladstone-dale coefficient for air
	rayTrace->absorption = absorptionCoeff; // Beer-Lambert coefficient (1/length) inside the liquid

	// Simulate image
	rayTrace->visibilityTrace();

	// clear memory
	rayTrace->meshOBBTree->Delete();

	//write image data array to .bmp
	rayTrace->writeBMP();
	cout << "Done." << endl;
	return 0;
}

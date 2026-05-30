#include <iostream>
#include "meshTools.h"
#include "vtkRaytrace.h"

int main ()
{
	// This is a back-lit illumination ray tracing code designed to handle transmissive objects.
	vtkRaytrace r = vtkRaytrace::vtkRaytrace(); // central object
    vtkRaytrace *rayTrace = &r;

	int nContours = 10;

	// Read in mesh file and convert to density iso-contours, then calculate normals and OBBTree
	rayTrace->mesh = vtkPolyData::New();
	rayTrace->volume = 0;
	rayTrace->subset = 1;

	char *inputFilename = new char [100];
	strcpy(inputFilename, "/Users/zack/Desktop/vtkImage/OpticalNozzle.stl");

    ifstream ifile(inputFilename);
    if(!ifile){printf("File does not exist!");}
    else{
	rayTrace->readMesh(inputFilename, nContours);
	meshTools::show(rayTrace->mesh);
	//rayTrace->plotNormals();
    }

	//cout << "Number of cells = " << rayTrace->mesh->GetNumberOfCells() << endl;

	//Define camera parameters
	rayTrace->focalPlane = 0; // define image plane distance from mesh centroid along optical axis
	rayTrace->focalLength = 60; // define the focal length of the camera system
	rayTrace->magnification = .9; // define magnification of the scene on the image
	
	// Camera chip parameters
	rayTrace->nPixels = new int [2];
	rayTrace->nPixels[0] = 1000;
	rayTrace->nPixels[1] = 1000;
	rayTrace->cellSizeFactor = .5; // ratio between pixel size and the largest mesh cell's image
	
	// define the Optical axis for the setup
	rayTrace->opticalAxis = new double [3];
	rayTrace->opticalAxis[0] = 1; 
	rayTrace->opticalAxis[1] = 0; 
	rayTrace->opticalAxis[2] = 0; 

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

	// Construct camera representation
	rayTrace->makeCamera(confineToBounds);

	// recursion limit for individual ray traces
	rayTrace->depthLimit = 7;

	// Light transport parameters
	rayTrace->calcIndex = 1;
	rayTrace->glassIndex = 1.5;
	rayTrace->K_GD = .00023; // gladstone-dale coefficient for air
	rayTrace->absorption = 0;

	// Simulate image
	rayTrace->visibilityTrace();

	// clear memory
	rayTrace->meshOBBTree->Delete();

	//write image data array to .bmp
	rayTrace->writeBMP();
}
	

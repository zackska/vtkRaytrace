#include <vtkCellCenters.h>
#include <vtkCellArray.h>
#include <vtkRenderer.h>
#include <vtkPolyData.h>
#include <vtkUnsignedCharArray.h>
#include <vtkSTLReader.h>
#include <vtkSmartPointer.h>
#include <vtkOBBTree.h>
#include <vtkStaticCellLocator.h>
#include <vtkAbstractCellLocator.h>
#include <vtkPointLocator.h>
#include "meshTools.h"
#include <vtkPoints.h>
#include <vtkLine.h>
#include <vtkMath.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkMatrix4x4.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkFloatArray.h>
#include <vtkPolyDataNormals.h>
#include <vtkSTLWriter.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>

#ifndef vtkRaytrace_H
#define vtkRaytrace_H

class  vtkRaytrace
{

//member functions
public:

vtkRaytrace(){};

~vtkRaytrace(){};

void makeCamera(bool confine);

void addRay(vtkRenderer *renderer, double *p1, double *p2);

void refract(double * norm, double * ray, double * newray, double indexRatio);

void getRay(double p1[3], double p2[3],double * ray);

void writeBMP();

void readMesh(char *inputFilename, int nContours);


void printProgress (double percentage);

void plotSurfaces(double (*p1)[3], int  n_points1, double (*p2)[3], int n_points2);


//functions for ray tracing calculation

double trace(double *source, double *target, int depth, vtkRenderer *&render, bool insideMedium = false);

void visibilityTrace();

void rasterize();

void plotNormals();



//"Camera" parameters
int *nPixels;
double chipLoc;
double *chipSize;
double **image;
vtkMatrix4x4 *world2Camera;
double magnification; // define magnification of the scene on the image
double focalLength; // define the focal length of the camera system
double focalPlane; // define image plane
double cellSizeFactor;
double *cameraPos; // define a 3-dimensional point for the centroid of the canvas
double *perspectivePoint;
char *fileType;

// Thin-lens objective: finite collection aperture + focus plane (vs the default pinhole).
double m_cameraX[3];      // in-plane (transverse) camera axes
double m_cameraZ[3];
double m_focusCenter[3];  // a point on the focus plane (normal = opticalAxis)
double m_apertureRadius;  // lens aperture radius (0 => pinhole)
double m_fNumber;         // f-number (focusDist/(2*aperture)); large => pinhole


// auxiliary variables used in tracing algorithm
int depthLimit;
double *ray;
bool plot;
double varianceLimit;
int initialSample;

// variables signifying file import options
bool subset;
bool volume;
bool calcIndex; // specify whether or not to calculate the refractive index of the mesh from field data

//scene and mesh parameters
vtkPolyData *mesh;
vtkAbstractCellLocator *meshOBBTree; // vtkStaticCellLocator (thread-safe queries)
double m_meshLength;      // cached mesh->GetLength() (GetLength/GetCenter write shared
double m_meshCenter[3];   // buffers -> not thread-safe to call inside the parallel loop)
vtkScalarsToColors *colorLookupTable;	
double *opticalAxis;// define the distance between the simulated image plane and the perspective point;
double glassIndex; // refractive index of transmissive material.
double lightLoc;
double lightRadius;
double *lightCenter;
double acceptanceAngle;
double K_GD;
double absorption;


private:

};

#endif

#include <vtkMath.h>
#include <vtkQuad.h>
#include <vtkTriangle.h>
#include <vtkIdTypeArray.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkPoints.h>
#include <vtkCellData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkCellDataToPointData.h>
#include <vtkStructuredGrid.h>
#include <vtkMarchingContourFilter.h>
#include <vtkFieldDataToAttributeDataFilter.h>
#include <vtkScalarsToColors.h>
#include <vtkAppendPolyData.h>
#include "dirent.h"
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkStructuredGridAppend.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLMultiBlockDataReader.h>
#include <vtkMultiBlockDataSet.h>

#ifndef meshTools_H
#define meshTools_H

class meshTools
{

public:

	
static void printVec(double* vec, const char *name, int length);
	
static void printVec(int* vec, const char *name, int length);

static void extractMultiblockSubset(char *inputFilename, double distance, vtkStructuredGrid *grid, char *outputFilename);

static void gridScalarContours(vtkStructuredGrid *grid, vtkPolyData *mesh, int nContours);

static void show(vtkSmartPointer<vtkRenderer> renderer,  vtkPolyData * mesh);

static void show(vtkPolyData * mesh);

static void show(vtkPolyData * mesh, vtkScalarsToColors *table);

static void displayAll(vtkPolyData * mesh, double chip_size, double imageLoc, int n_pixels);

static void readList(vtkIdTypeArray * array, int floor);

static void getCellNormal(vtkPolyData * mesh, int cellNum, double * norm);

static void getCellPoints(vtkPolyData * mesh, int cellNum, double (*&p1)[3], int  n_points);

static void cellFormFactor( vtkPolyData* mesh, int cell1, int cell2, double * ff);

static void plotSurfaces(double (*p1)[3], int  n_points1, double (*p2)[3], int n_points2);

};

#endif

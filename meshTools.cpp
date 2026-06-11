#include "meshTools.h"



// auxiliary


void meshTools::printVec(double* vec, const char *name, int length)
{
	cout << name << " = " ;
	for(int i=0;i<length;i++)
		cout << vec[i] << " ";
	
	cout << endl;
}

void meshTools::printVec(int* vec, const char *name, int length)
{
	cout << name << " = " ;
	for(int i=0;i<length;i++)
		cout << vec[i] << " ";
	
	cout << endl;
}


void meshTools::extractMultiblockSubset(char *inputFilename, double distance, vtkStructuredGrid *grid, char *outputFilename)
{
	// Iterate through all sub-grids, and add only those with centroids within a certain distance from the origin
		
		//vtkMultiBlockDataSet *gridSet = vtkMultiBlockDataSet::New();
		vtkXMLMultiBlockDataReader *MBreader = vtkXMLMultiBlockDataReader::New();
		MBreader->SetFileName(inputFilename);
		//time_t start = std::time(0);
		//MBreader->Update();
		
		//time_t end = std::time(0);
		//double time = difftime(end, start);
		//cout << time/3600 << "hours to read multiblock" << endl;
		//gridSet->DeepCopy(MBreader->GetOutput());
		
		//cout << gridSet->GetNumberOfBlocks() << endl;

		//vtkStructuredGridAppend *gridAppend = vtkStructuredGridAppend::New();

		std::vector<vtkStructuredGrid*> inputs;
		vtkXMLStructuredGridReader* reader = vtkXMLStructuredGridReader::New();
		vtkXMLStructuredGridWriter* writer = vtkXMLStructuredGridWriter::New();

		DIR *pDIR;
		struct dirent *entry;
		int inputCount = 0;
		int maxCount = 100;
		int count = 0;

		if( (pDIR=opendir(inputFilename)) ){
			while((entry = readdir(pDIR)))
			{
				if( strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0 )
				{
					// read .vts file
					char *tempFilename = new char [100];
					memcpy(tempFilename, inputFilename, sizeof(char)*100);
					strcat(tempFilename, "/");
					strcat(tempFilename, entry->d_name);
					reader->SetFileName(tempFilename);
					reader->Update();
					
					// store in tempgrid
					vtkStructuredGrid *tempgrid = vtkStructuredGrid::New();
					tempgrid->DeepCopy(reader->GetOutput());

					// if grid centroid is close enought to origin, add it to the input list for the append filter
					if(vtkMath::Norm(tempgrid->GetCenter()) < distance)	
					{
						cout <<  entry->d_name << endl;
						//inputs.push_back(tempgrid);
						//inputCount++;
					}

					count ++;
				}
				
			}
			closedir(pDIR);
		}
		
		for(int i=0;i<inputCount; i++)
		{
			//gridAppend->AddInputData(inputs[i]);
		}

		//cout << gridAppend->GetNumberOfInputs() << endl;
		//gridAppend->Update();


		//grid->DeepCopy(gridAppend->GetOutput());
		double *bounds = new double [6];
		grid->GetBounds(bounds);
		meshTools::printVec(bounds, "Bounds", 6);
		
		// save subset as structured grid file
		strcat(inputFilename, "/");
		strcat(inputFilename, outputFilename);
		writer->SetInputData(grid);
		writer->SetFileName(inputFilename);
		writer->Update();
}

void meshTools::readList(vtkIdTypeArray * array, int floor)
{
	for(int i = 0; i<array->GetNumberOfTuples();i++)
	{
		double * thisTuple = array->GetTuple(i);
		if(*thisTuple > floor)
		{
			cout << *thisTuple << endl;
			cout << i<< endl;
		}
	}
}

void meshTools::displayAll(vtkPolyData *mesh, double chip_size, double imageLoc, int n_pixels)
{
	// extract points and cells from original mesh
	vtkSmartPointer<vtkCellArray> meshCells = 
	vtkSmartPointer<vtkCellArray>::New();
	meshCells = mesh->GetPolys();

	vtkSmartPointer<vtkPoints>	meshPoints =
	vtkSmartPointer<vtkPoints>::New();
	meshPoints = mesh->GetPoints();

	
	//Add new points
	
	vtkQuad *pixel = vtkQuad::New();

	double pixel_size = chip_size/n_pixels;

	int nOrigPoints = mesh->GetNumberOfPoints();

	for(int i=0;i<n_pixels+1;i++)
	{
		for(int j=0;j<n_pixels+1; j++)
			meshPoints->InsertNextPoint(-chip_size/2+j*pixel_size, imageLoc, -chip_size/2+i*pixel_size);			
	}
	
	int *points = new int;
	//Add new cells
	for(int i=0;i<n_pixels;i++)
	{
		for(int j=0;j<n_pixels; j++)
		{ 

			*points = j+i*(n_pixels+1);
			*(points+1) = j+1+i*(n_pixels+1);
			*(points+2) = j+1+(i+1)*(n_pixels+1);
			*(points+3) = j+(i+1)*(n_pixels+1);

			for(int k=0;k<4;k++) // add points to cells
			{ 
				pixel->GetPointIds()->SetId(k, *(points+k)+nOrigPoints);
			} 

			meshCells->InsertNextCell(pixel);
		}		
	}
	
	vtkPolyData *newMesh = vtkPolyData::New();
	// construct new mesh with chip cells added
	newMesh->SetPoints(meshPoints);
	newMesh->SetPolys(meshCells);

	meshTools::show(mesh);

	newMesh->Delete();
}

void meshTools::gridScalarContours(vtkStructuredGrid *grid, vtkPolyData *mesh, int nContours)
{
	// Assign refractive-index field as scalar to contour with. NOTE: read "rho_s" from the
	// input GRID (was mesh, the empty output polydata -> null). rho_s must hold n = 1+K_GD*rho.
		grid->GetCellData()->SetScalars(grid->GetCellData()->GetArray("rho_s"));

		// Transfer density from cells to points
		vtkCellDataToPointData *cell2point = vtkCellDataToPointData::New();
		cell2point->SetInputData(grid);
		cell2point->Update();
		grid->DeepCopy(cell2point->GetOutput());

		// have to transfer from field data to attribute data for the marching contour filter. 
		vtkFieldDataToAttributeDataFilter *field2attribute = vtkFieldDataToAttributeDataFilter::New();
		field2attribute->SetInputData(grid);
		field2attribute->SetInputFieldToPointDataField();
		field2attribute->SetScalarComponent(0, "rho_s", 0);
		field2attribute->SetOutputAttributeDataToPointData();
		field2attribute->Update();
		grid->DeepCopy(field2attribute->GetOutput());

		// Extract contour from data set
		double *range = new double [2];
		grid->GetPointData()->GetScalars()->GetRange(range);
		//printVec(range, "Min/Max density = ", 2);

		vtkAppendPolyData* append = vtkAppendPolyData::New();

		vtkMarchingContourFilter *filter = vtkMarchingContourFilter::New();
		filter->UseScalarTreeOn();
		filter->SetInputData(grid);
		filter->GenerateValues(nContours, range);
		double *contourValues = filter->GetValues();
		filter->SetNumberOfContours(1); 

		// Create the color map
		vtkScalarsToColors *colorLookupTable = vtkScalarsToColors::New();
		colorLookupTable->SetRange(range);
		colorLookupTable->Build();
 
		//isolate density contours one by one and assign refractive index ratio to them
		for(int k = 1; k<nContours-1;k++)
		{
			filter->SetValue(0, contourValues[k]);
			filter->Update();

			vtkPolyData *tempMesh = vtkPolyData::New();
			tempMesh->DeepCopy(filter->GetOutput());

			// PER-SHELL arrays (were allocated once outside the loop and resized each
			// iteration -> every shell ended up sharing the LAST shell's array, corrupting
			// the index ratios and overrunning on append). Allocate fresh each shell.
			vtkDoubleArray* index = vtkDoubleArray::New();
			index->SetNumberOfComponents(1); index->SetName("index ratio");
			index->SetNumberOfTuples(tempMesh->GetNumberOfCells());
			index->FillComponent(0, contourValues[k]/contourValues[k+1]);
			tempMesh->GetCellData()->AddArray(index);
			index->Delete();   // tempMesh keeps a reference

			// add colors for plotting (dcolor was `new double` (1 elem) used as a 3-vector)
			vtkIntArray *colors = vtkIntArray::New();
			colors->SetNumberOfComponents(3); colors->SetName("Colors");
			double dcolor[3] = {0,0,0};
			colorLookupTable->GetColor(contourValues[k], dcolor);
			colors->SetNumberOfTuples(tempMesh->GetNumberOfPoints());
			for(int l=0;l<3;l++)
			{
				colors->FillComponent(l, int(dcolor[l]*255));
			}

			tempMesh->GetCellData()->SetScalars(colors);
			colors->Delete();   // tempMesh keeps a reference
			// add to appendPolyData filter
			append->AddInputData(tempMesh);
			tempMesh->Delete();
		}
	
		append->Update();
		// extract all contours with associated density ratios
		mesh->DeepCopy(append->GetOutput());

	
		//free memory
		filter->Delete();
		field2attribute->Delete();
		cell2point->Delete();
}

//display

void meshTools::show(vtkSmartPointer<vtkRenderer> renderer,  vtkPolyData * mesh)
{
	vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(mesh);

	vtkSmartPointer<vtkActor> actor =
	vtkSmartPointer<vtkActor>::New();
	actor->GetProperty()->SetOpacity(0.25);
	actor->SetMapper(mapper);
 
 
	vtkSmartPointer<vtkRenderWindow> renderWindow =
	vtkSmartPointer<vtkRenderWindow>::New();

	renderWindow->AddRenderer(renderer);

	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
	vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	
	vtkSmartPointer<vtkAxesActor> axes = 
    vtkSmartPointer<vtkAxesActor>::New();
 
	vtkSmartPointer<vtkOrientationMarkerWidget> widget = 
	vtkSmartPointer<vtkOrientationMarkerWidget>::New();
	widget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
	widget->SetOrientationMarker( axes );
	widget->SetInteractor( renderWindowInteractor );
	widget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
	widget->SetEnabled( 1 );
	widget->InteractiveOn();
 
 
	renderer->AddActor(actor);
	renderer->SetBackground(1, 1, 1); // Background color white

	
	renderer->ResetCamera();
	renderWindow->Render();
	renderWindowInteractor->Start();
 }

void meshTools::show(vtkPolyData * mesh)
{
	vtkSmartPointer<vtkRenderer> renderer = 
		vtkSmartPointer<vtkRenderer>::New();

    vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
    
	mapper->SetInputData(mesh);

	vtkSmartPointer<vtkActor> actor =
	vtkSmartPointer<vtkActor>::New();
	//actor->GetProperty()->SetOpacity(0.25);
	actor->SetMapper(mapper);
 
 
	vtkSmartPointer<vtkRenderWindow> renderWindow =
	vtkSmartPointer<vtkRenderWindow>::New();

	renderWindow->AddRenderer(renderer);

	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
	vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	vtkSmartPointer<vtkAxesActor> axes = 
    vtkSmartPointer<vtkAxesActor>::New();
 
	vtkSmartPointer<vtkOrientationMarkerWidget> widget = 
	vtkSmartPointer<vtkOrientationMarkerWidget>::New();
	widget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
	widget->SetOrientationMarker( axes );
	widget->SetInteractor( renderWindowInteractor );
	widget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
	widget->SetEnabled( 1 );
	widget->InteractiveOn();
 
	renderer->AddActor(actor);
	renderer->SetBackground(1, 1, 1); // Background color black

	
	renderer->ResetCamera();
	renderWindow->Render();
	renderWindowInteractor->Start();
 }



void meshTools::show(vtkPolyData * mesh, vtkScalarsToColors*table)
{
	vtkSmartPointer<vtkRenderer> renderer = 
		vtkSmartPointer<vtkRenderer>::New();

	vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(mesh);
	mapper->SetLookupTable(table);

	vtkSmartPointer<vtkActor> actor =
	vtkSmartPointer<vtkActor>::New();
	//actor->GetProperty()->SetOpacity(0.25);
	actor->SetMapper(mapper);
 
 
	vtkSmartPointer<vtkRenderWindow> renderWindow =
	vtkSmartPointer<vtkRenderWindow>::New();

	renderWindow->AddRenderer(renderer);

	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
	vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	vtkSmartPointer<vtkAxesActor> axes = 
    vtkSmartPointer<vtkAxesActor>::New();
 
	vtkSmartPointer<vtkOrientationMarkerWidget> widget = 
	vtkSmartPointer<vtkOrientationMarkerWidget>::New();
	widget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
	widget->SetOrientationMarker( axes );
	widget->SetInteractor( renderWindowInteractor );
	widget->SetViewport( 0.0, 0.0, 0.4, 0.4 );
	widget->SetEnabled( 1 );
	widget->InteractiveOn();
 
	renderer->AddActor(actor);
	renderer->SetBackground(1, 1, 1); // Background color black

	
	renderer->ResetCamera();
	renderWindow->Render();
	renderWindowInteractor->Start();
 }

//property extraction

void meshTools::getCellPoints(vtkPolyData *mesh, int cellNum, double (*&p1)[3], int n_points)
{
	vtkIdList* pts = vtkIdList::New();
	mesh->GetCellPoints(cellNum, pts);
	n_points = pts->GetNumberOfIds();
				
	p1 = new double [n_points][3];
				
	for(int j = 0;j<n_points;j++) // go through the points in the target cell and 
	{
		double point[3];
		mesh->GetPoint(pts->GetId(j), point);

		for(int k=0;k<3;k++)
		{
			p1[j][k] = point[k];
		}
	}
}

void meshTools::getCellNormal(vtkPolyData *mesh, int cellNum, double * norm)
{
		 vtkSmartPointer<vtkDataArray> normals = mesh->GetCellData()->GetNormals();
		 normals->GetTuple(cellNum, norm);
		 vtkMath::Normalize(norm);
}

//form factors

void meshTools::cellFormFactor(vtkPolyData* mesh, int cell1, int cell2, double * ff)
{
				double (*p1)[3] = new double [1][3];
				int n_point1 = 0;
				meshTools::getCellPoints(mesh, cell1, p1, n_point1);

				double (*p2)[3] = new double [1][3];
				int n_point2 = 0;
				meshTools::getCellPoints(mesh, cell2, p2, n_point2);	

				//*ff = FormFactor(p1, n_point1, p2, n_point2);
}

void meshTools::plotSurfaces(double (*p1)[3], int  n_points1, double (*p2)[3], int n_points2)
{

	vtkSmartPointer<vtkPoints>	meshPoints =
	vtkSmartPointer<vtkPoints>::New();
	

	vtkSmartPointer<vtkTriangle> face1 =
	vtkSmartPointer<vtkTriangle>::New();
	if(n_points1 ==4)
	{
	vtkSmartPointer<vtkQuad> face1 =
	vtkSmartPointer<vtkQuad>::New();
	}

	for(int k=0;k<n_points1;k++)
	{
		face1->GetPointIds()->SetId( k, k);
		meshPoints->InsertNextPoint(p1[k][0], p1[k][1], p1[k][2]);
	}


	vtkSmartPointer<vtkTriangle> face2 =
	vtkSmartPointer<vtkTriangle>::New();
	if(n_points2 ==4)
	{
		vtkSmartPointer<vtkQuad> face2 =
		vtkSmartPointer<vtkQuad>::New();
	}

	for(int k=0;k<n_points2;k++)
	{
		face2->GetPointIds()->SetId( k, k+n_points1);
		meshPoints->InsertNextPoint(p2[k][0], p2[k][1], p2[k][2]);
	}

	vtkSmartPointer<vtkCellArray> meshCells = 
		vtkSmartPointer<vtkCellArray>::New();
	meshCells->InsertNextCell(face1);
	meshCells->InsertNextCell(face2);

	
	vtkPolyData * newMesh = 
	vtkPolyData::New();
	newMesh->SetPoints(meshPoints);
	newMesh->SetPolys(meshCells);

	vtkSmartPointer<vtkRenderer> renderer = 
		vtkSmartPointer<vtkRenderer>::New();
	
	meshTools::show(renderer, newMesh);


}






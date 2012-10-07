// Image types and other standard things
#include <iostream>
#include <cstdlib>
#include "vtkUnstructuredGrid.h"
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include "itkTimeVaryingBSplineVelocityFieldTransform.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkImageFileWriter.h"
#include <itkPointSet.h>
#include <itkImageRegionIterator.h>

#include "itkTransformToDisplacementFieldSource.h"

int main( int argc, char * argv[] )
{
	if (argc < 3)
	{
		std::cerr << "Usage : " << argv[0] << "InputMeshList.txt outputMeshPrefix" << std::endl;
		exit(EXIT_FAILURE);
	}
	const unsigned int MAX_PATH_LENGTH = 1024;
	const unsigned int ImageDimension = 3;

	int arg = 1;
	const char* InputMeshListFileName = argv[arg]; arg++;
	const char* outputMeshPrefixFileName = argv[arg]; arg++;

	// reading input list
	std::ifstream f (InputMeshListFileName, std::ifstream::in);
	if (f==0)
	{
		return EXIT_FAILURE;
	}

	std::string buffer;
	std::vector<std::string> inputUGFileNames;
	f>>buffer;

	while (!f.eof())
	{
		inputUGFileNames.push_back(buffer);
		f>>buffer;
	}

	// reading all meshes
	std::vector<vtkSmartPointer<vtkUnstructuredGrid> > ugs;
	for (unsigned int t=0; t<inputUGFileNames.size(); t++)
	{
		vtkSmartPointer<vtkUnstructuredGridReader> reader
			= vtkSmartPointer<vtkUnstructuredGridReader>::New();
		reader->SetFileName(inputUGFileNames[t].c_str());
		reader->Update();
		ugs.push_back(reader->GetOutput());
	}

	// check that we have at least to meshes
	if (ugs.size() <=1)
	{
		std::cerr << "Insufficient (" << ugs.size() << ") number of meshes" << std::endl;
		exit(EXIT_FAILURE);
	}

	// Instantiate velocity and adjust bspline field to mesh bounding box
	typedef itk::TimeVaryingBSplineVelocityFieldTransform< double, ImageDimension > TransformType;
	typedef TransformType::DisplacementVectorType DisplacementVectorType;
	typedef TransformType::TimeVaryingVelocityFieldControlPointLatticeType VelocityFieldType;
	typedef itk::Image<DisplacementVectorType, ImageDimension + 1> DisplacementFieldType;
	typedef itk::PointSet<DisplacementVectorType, ImageDimension + 1> PointSetType;
	typedef itk::BSplineScatteredDataPointSetToImageFilter<PointSetType, VelocityFieldType> BSplineFilterType;

	// Filling point set with velocity (incr. displ.) data
	PointSetType::Pointer velocityFieldPoints = PointSetType::New();
    velocityFieldPoints->Initialize();

	typedef BSplineFilterType::WeightsContainerType WeightsContainerType;
    WeightsContainerType::Pointer velocityFieldWeights = WeightsContainerType::New();

	PointSetType::PointType spatioTemporalPoint;
	DisplacementVectorType velocity;
	PointSetType::PointIdentifier id = 0;
	for (unsigned int t=0; t<ugs.size()-1; t++)
	{
		for (int p=0; p<ugs[t]->GetNumberOfPoints(); p+=50,id++)
		{
			for (int d=0; d<ImageDimension; d++)
			{
				spatioTemporalPoint[d] = ugs[t]->GetPoints()->GetPoint(p)[d];
				velocity[d] =
				(   ugs[t+1]->GetPoints()->GetPoint(p)[d]
				-	ugs[t  ]->GetPoints()->GetPoint(p)[d] )
				* (double) ugs.size();
			}
			spatioTemporalPoint[ImageDimension] = (double)t / (double) ugs.size();

			velocityFieldPoints->SetPoint(id,  spatioTemporalPoint);
			velocityFieldPoints->SetPointData(id, velocity );
			velocityFieldWeights->InsertElement( id, 1. );
		}
	}

	// Setting up Bspline regression of the velocity field
	BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
	unsigned int nocp[4] = {20,20,20,10};
	BSplineFilterType::SizeType domainSize;
	BSplineFilterType::ArrayType numberOfControlPoints(nocp);
	BSplineFilterType::PointType origin;
	BSplineFilterType::SpacingType spacing;

 double bb[6];
	ugs[0]->GetBounds(bb);

	spacing.Fill( 1.0 );
	domainSize[0] = (int)( bb[1] - bb[0] + 0.5 ) + 3;
	domainSize[1] = (int)( bb[3] - bb[2] + 0.5 ) + 3;
	domainSize[2] = (int)( bb[5] - bb[4] + 0.5 ) + 3;
	domainSize[3] = ugs.size() + 3;

	origin[0] = bb[0] - spacing[0];
	origin[1] = bb[2] - spacing[1];
	origin[2] = bb[4] - spacing[2];
	origin[3] = 0 - spacing[3];

	bspliner->SetOrigin( origin );
	bspliner->SetSpacing( spacing );
	bspliner->SetNumberOfControlPoints( numberOfControlPoints );
   	bspliner->SetNumberOfLevels( 1 );
	bspliner->SetSplineOrder( 3 );
	bspliner->SetSize( domainSize );
	bspliner->SetNumberOfLevels( 1 );
	bspliner->SetInput( velocityFieldPoints );
	bspliner->SetPointWeights( velocityFieldWeights );
	bspliner->GenerateOutputImageOn();

	try
	{
		bspliner->Update();
	}
	catch( itk::ExceptionObject & err )
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

// 	typedef itk::ImageFileWriter<VelocityFieldType> WriterType;
// 	WriterType::Pointer writer = WriterType::New();
// 	writer->SetFileName( "velocityField.nii.gz" );
// 	writer->SetInput( bspliner->GetOutput() );
// 	writer->Update();

	TransformType::TimeVaryingVelocityFieldControlPointLatticePointer lattice = bspliner->GetPhiLattice();
// 	lattice->FillBuffer(10.);// !!!!!!!!!!!!!!!!!!!!!!!!!!

	// Transformation that is supposed to approximate our velocity field
	TransformType::Pointer transform = TransformType::New();
	transform->SetVelocityFieldOrigin( origin );// Why is lattice->GetOrigin() giving an incorrect value?
	transform->SetVelocityFieldSize( domainSize );
	transform->SetVelocityFieldSpacing( spacing );
	transform->SetSplineOrder( 3 );
	transform->SetNumberOfIntegrationSteps( 4 );
	transform->SetTimeVaryingVelocityFieldControlPointLattice( lattice );

	// Allocating
	std::vector<vtkSmartPointer<vtkUnstructuredGrid> > out_ugs;
	for (unsigned int t=0; t<ugs.size(); t++)
	{
		vtkSmartPointer<vtkUnstructuredGrid> oug = vtkSmartPointer<vtkUnstructuredGrid>::New();
		oug->DeepCopy(ugs[0]);
		out_ugs.push_back(oug);
	}

	// Recomputing transform for all points
	TransformType::OutputPointType outp;
	for (unsigned int t=1; t<ugs.size(); t++)
	{
	 std::cout << "Transform " << t << std::endl;

		transform->SetUpperTimeBound((double)t / (double)ugs.size());
		transform->SetLowerTimeBound((double)(t-1)/ (double)ugs.size());
		transform->SetNumberOfIntegrationSteps( 4 );
		transform->Modified();
		transform->IntegrateVelocityField();

		for (int p=0; p<ugs[0]->GetNumberOfPoints(); p++)
		{
			for (int d=0; d<ImageDimension; d++)
			{
				spatioTemporalPoint[d] = ugs[t-1]->GetPoints()->GetPoint(p)[d];
			}
			spatioTemporalPoint[ImageDimension] = 0.;
			outp = transform->TransformPoint(ugs[0]->GetPoints()->GetPoint(p));
			out_ugs[t]->GetPoints()->SetPoint(p, outp.GetDataPointer());
		}
	}

	// Writing output
	char outputFileName[MAX_PATH_LENGTH];
	for (unsigned int t=0; t<ugs.size(); t++)
	{
		sprintf(outputFileName, "%s%03d.vtk", outputMeshPrefixFileName, t);
		vtkSmartPointer<vtkUnstructuredGridWriter> wr = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
		wr->SetFileName(outputFileName);
		wr->SetInput(out_ugs[t]);
		wr->Update();
	}


	exit(EXIT_SUCCESS);
}



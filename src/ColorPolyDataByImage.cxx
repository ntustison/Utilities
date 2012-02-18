#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkMarchingContourFilter.h"

#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkImageFileReader.h"

int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " inputPolyData inputImage outputPolyData" << std::endl;
    return EXIT_FAILURE;
    }

  vtkPolyDataReader *dataReader = vtkPolyDataReader::New();
  dataReader->SetFileName( argv[1] );
  dataReader->Update();

  vtkPolyData *polydata = dataReader->GetOutput();
  vtkPoints *points = polydata->GetPoints();
  vtkFloatArray *scalars = vtkFloatArray::New();
  scalars->SetNumberOfValues( polydata->GetNumberOfPoints() );

  typedef itk::Image<float, 3> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName( argv[2] );
  imageReader->Update();

  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, float> InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage( imageReader->GetOutput() );

  for( unsigned int n = 0; n < points->GetNumberOfPoints(); n++ )
    {
    double *coords = points->GetPoint( n );

    InterpolatorType::PointType point;
    for( unsigned int d = 0; d < 3; d++ )
      {
      point[d] = coords[d];
      }
    InterpolatorType::OutputType value = interpolator->Evaluate( point );

    std::cout << point << " -> " << value << std::endl;
    scalars->SetValue( n, value );
    }
  polydata->GetPointData()->SetScalars( scalars );

  vtkPolyDataWriter *writerF = vtkPolyDataWriter::New();
  writerF->SetInput( polydata );
//  writerF->SetFileTypeToBinary();
  writerF->SetFileName( argv[3] );
  writerF->Write();

  return 0;
}

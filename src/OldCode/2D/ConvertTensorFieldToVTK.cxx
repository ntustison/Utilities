#include "itkVectorImageFileReader.h"
#include "itkVector.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "vtkStructuredGrid.h"
#include "vtkStructuredGridWriter.h"
#include "vtkFloatArray.h"
#include "vtkPoints.h"
#include "vtkPointData.h"


#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " inputTensorField outputVTKFile [slice] [whichAxis]" << std::endl;
    exit( 1 );
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  
  typedef double RealType;
  typedef itk::Vector<RealType, 6> TensorType;
  typedef itk::Image<TensorType, ImageDimension> TensorFieldType;
  
  typedef itk::VectorImageFileReader<ImageType, TensorFieldType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  double origin[ImageDimension];
  double spacing[ImageDimension];
  int    size[ImageDimension];
  int totalsize = 1;

  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    origin[i] = reader->GetOutput()->GetOrigin()[i];
    spacing[i] = reader->GetOutput()->GetSpacing()[i];
    size[i] = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[i] - 1;
    if ( argc > 3 && atoi( argv[4] ) == static_cast<int>( i ) )
      {
      size[i] = 1;
      }              
    totalsize *= size[i];
    }


  vtkStructuredGrid *field = vtkStructuredGrid::New();
  field->SetDimensions( size );
  vtkFloatArray *tensors = vtkFloatArray::New();
  tensors->SetNumberOfComponents( 9 );
  tensors->SetNumberOfTuples( totalsize );
  vtkPoints *points = vtkPoints::New();
  points->Allocate( totalsize );

  float x[ImageDimension], t[9];
  int offset;
  
  offset = 0;

  itk::ImageRegionIteratorWithIndex<TensorFieldType> It
    ( reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    TensorFieldType::IndexType idx = It.GetIndex();
    
    if ( argc > 3 && idx[atoi( argv[4] )] != atoi( argv[3] ) )
      {
      continue;  
      }              
    TensorFieldType::PointType point;
    reader->GetOutput()->TransformIndexToPhysicalPoint( idx, point ); 
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      x[i] = point[i];
      }
    
    t[0] = It.Get()[0]; 
    t[1] = t[3] = It.Get()[1];
    t[4] = It.Get()[2]; 
    t[2] = t[6] = It.Get()[3]; 
    t[5] = t[7] = It.Get()[4]; 
    t[8] = It.Get()[5]; 

    points->InsertPoint( offset, x );
    tensors->InsertTuple( offset++, t );
    }
  field->SetPoints( points );
  points->Delete();
  field->GetPointData()->SetTensors( tensors );
  tensors->Delete();

  vtkStructuredGridWriter *writer = vtkStructuredGridWriter::New();
  writer->SetInput( field );
  writer->SetFileTypeToBinary();
  writer->SetFileName( argv[2] );
  writer->Write();

  return 0;
}

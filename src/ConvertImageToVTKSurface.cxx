#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"

#include "itkImage.h"
#include "itkVector.h"
#include "itkVectorImageFileReader.h"
#include "itkNeighborhoodIterator.h"

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " inputFile outputFile" << std::endl;
    exit( 1 );
    }

  const unsigned int ImageDimension = 2;

  typedef float RealType;
  typedef itk::Image<RealType, ImageDimension> ImageType;
  typedef itk::Vector<RealType, 3> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;

  typedef itk::VectorImageFileReader<ImageType, DeformationFieldType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->SetUseAvantsNamingConvention( true );
  reader->Update();

  int totalPoints = 1;
  int totalCells = 1;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    totalPoints *= reader->GetOutput()->GetLargestPossibleRegion().GetSize()[i];
    totalCells *= ( reader->GetOutput()->GetLargestPossibleRegion().GetSize()[i] - 1 );
    } 

  vtkPolyData *surface = vtkPolyData::New();
  surface->Allocate( totalCells );
  vtkPoints *points = vtkPoints::New();
  points->Allocate( totalPoints );
  vtkCellArray *polys = vtkCellArray::New();
  
  float x[3];
  vtkIdType pts[4];
  
  unsigned long count = 0;

  typedef itk::NeighborhoodIterator<DeformationFieldType> NeighborhoodIteratorType;
  NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It ( radius, 
    reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    VectorType pt = It.GetCenterPixel();
    x[0] = pt[0];  
    x[1] = pt[1];  
    x[2] = pt[2];  

    points->InsertPoint( count++, x );

    if ( !reader->GetOutput()->GetLargestPossibleRegion().IsInside( It.GetIndex( 2 ) ) )
      {
      continue;
      }  
    ImageType::IndexType index = It.GetIndex( 1 );  
    pts[0] = index[0] + 
      index[1]*reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0];  
    index = It.GetIndex( 2 );  
    pts[1] = index[0] + 
      index[1]*reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0];  
    index = It.GetIndex( 5 );  
    pts[2] = index[0] + 
      index[1]*reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0];  
    index = It.GetIndex( 4 );  
    pts[3] = index[0] + 
      index[1]*reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0];  
    
    polys->InsertNextCell( 4, pts );
    }
  surface->SetPoints( points );
  surface->SetPolys( polys );
  points->Delete();
  polys->Delete();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput( surface );
//  writer->SetFileTypeToBinary();
  writer->SetFileName( argv[2] );
  writer->Write();

  return 0;
}

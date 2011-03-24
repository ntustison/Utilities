#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"

#include "itkImage.h"
#include "itkVector.h"
#include "itkImageFileReader.h"
#include "itkNeighborhoodIterator.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " inputFile outputFile [scalar]" << std::endl;
    exit( 1 );
    }

  typedef float RealType;
  typedef itk::Image<PixelType, 2> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  int totalPoints = 1;
  int totalCells = 1;
  for ( unsigned int i = 0; i < 2; i++ )
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

  typedef itk::NeighborhoodIterator<ImageType> NeighborhoodIteratorType;
  NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It ( radius, 
    reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    x[0] = reader->GetOutput()->GetOrigin()[0] 
      + reader->GetOutput()->GetSpacing()[0] * It.GetIndex()[0];  
    x[1] = reader->GetOutput()->GetOrigin()[1] 
      + reader->GetOutput()->GetSpacing()[1] * It.GetIndex()[1];  
    x[2] = It.GetCenterPixel();
    if ( argc == 4 )
      {
      x[2] *= atof( argv[3] );
      }    

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

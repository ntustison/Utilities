#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"

#include "itkImage.h"
#include "itkVector.h"
#include "itkImageFileReader.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNeighborhoodIterator.h"

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " inputFile outputFile [scalar] [outputGridFile.txt] "
              << "[numberOfGridLinesX] [numberOfGridLinesY]" << std::endl;
    exit( 1 );
    }

  typedef float RealType;
  typedef itk::Image<RealType, 2> ImageType;

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
  float factor = ( ( argc < 4 ) ? 1.0 : atof( argv[3] ) );

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
    x[2] = It.GetCenterPixel() * factor;

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

  if( argc < 5 )
    {
    return 0;
    }

  /**
   * Now write the mesh grid
   */

  typedef itk::LinearInterpolateImageFunction<ImageType, RealType> InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage( reader->GetOutput() );

  ofstream str2( argv[4] );

  str2 << "0 0 0 0" << std::endl;

  ImageType::PointType maxBoundary;
  ImageType::SpacingType gridSpacing;
  for ( unsigned int d = 0; d < 2; d++ )
    {
    maxBoundary[d] = ( reader->GetOutput()->GetLargestPossibleRegion().GetSize()[d] - 1 )
      * reader->GetOutput()->GetSpacing()[d] + reader->GetOutput()->GetOrigin()[d];

    unsigned int numberOfGridLines = 10;
    if ( static_cast<unsigned>( argc ) > 5 + d )
      {
      numberOfGridLines = static_cast<unsigned int>( atoi( argv[d+5] ) );
      }

    gridSpacing[d] = ( maxBoundary[d] - reader->GetOutput()->GetOrigin()[d] )
      / static_cast<RealType>( numberOfGridLines - 1 ) - 0.0001;
    }

  unsigned int lineCount = 0;
  for ( unsigned int d = 0; d < 2; d++ )
    {
    RealType delta = reader->GetOutput()->GetSpacing()[d];

    InterpolatorType::PointType point;

    for ( unsigned int c = 0; c < 2; c++ )
      {
      point[c] = reader->GetOutput()->GetOrigin()[c];
      }

    while ( true )
      {
      InterpolatorType::OutputType scale = interpolator->Evaluate( point );

      str2 << point[0] << " " << point[1] << " " << scale*factor << " "
        << lineCount+1 << std::endl;

      point[d] += delta;

      unsigned int doneFlag = 0;
      if ( point[d] > maxBoundary[d] )
        {

        lineCount++;
        point[d] = reader->GetOutput()->GetOrigin()[d];
        doneFlag++;

        for ( unsigned int c = 0; c < 1; c++ )
          {
          point[(d+c+1)%2] += gridSpacing[(d+c+1)%2];
          if ( point[(d+c+1)%2] > maxBoundary[(d+c+1)%2] )
            {
            point[(d+c+1)%2] = reader->GetOutput()->GetOrigin()[(d+c+1)%2];
            doneFlag++;
            }
          else
            {
            break;
            }
          }
        }
      if ( doneFlag == 2 )
        {
        break;
        }
      }
    }
  str2 << "0 0 0 0" << std::endl;
  str2.close();



  return 0;
}

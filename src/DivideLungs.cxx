#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkBoundingBox.h"

int main( int argc, char *argv[] )
{
  // Given an image where the left and right lungs have label '2' and '3',
  // the output image has left lung labels "1...numberOfLeftLungDivisions"
  // and right lung labels "numberOfLeftLungDivisions+1+(1...numberOfRightLungDivisions).

  if ( argc < 5 )
    {
    std::cout << "Usage: " << argv[0] << " lobeLabelImage outputImage "
              << "numberOfLeftLungDivisions numberOfRightLungDivisions [axis] "
              << "[lobeLabel1] [lobeLabel2]" << std::endl;
    exit( 1 );
    }

  typedef int PixelType;
  const unsigned int ImageDimension = 3;
  typedef double RealType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::BoundingBox<unsigned long,
       ImageDimension, RealType> BoundingBoxType;
  BoundingBoxType::Pointer box1 = BoundingBoxType::New();
  BoundingBoxType::Pointer box2 = BoundingBoxType::New();
  BoundingBoxType::PointsContainerPointer Points1
       = BoundingBoxType::PointsContainer::New();
  BoundingBoxType::PointsContainerPointer Points2
       = BoundingBoxType::PointsContainer::New();
  itk::Point<double, ImageDimension> point;

  int lobeLabel1 = 2;
  if( argc > 6 )
    {
    lobeLabel1 = atoi( argv[6] );
    }
  int lobeLabel2 = 3;
  if( argc > 7 )
    {
    lobeLabel2 = atoi( argv[7] );
    }


  int idx1 = 0;
  int idx2 = 0;
  itk::ImageRegionIteratorWithIndex<ImageType> It(
    reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( It.Get() == lobeLabel1 )
      {
      for (unsigned int i = 0; i < ImageDimension; i++ )
        {
        point[i] = static_cast<RealType>( It.GetIndex()[i] );
        }
      Points1->InsertElement( idx1++, point );
      }
    else if ( It.Get() == lobeLabel2 )
      {
      for (unsigned int i = 0; i < ImageDimension; i++ )
        {
        point[i] = static_cast<RealType>( It.GetIndex()[i] );
        }
      Points2->InsertElement( idx2++, point );
      }
    }
  box1->SetPoints( Points1 );
  box1->ComputeBoundingBox();
  box2->SetPoints( Points2 );
  box2->ComputeBoundingBox();

  int axis = 2;
  if( argc > 5 )
    {
    axis = atoi( argv[5] );
    }

  RealType dz1 = ( box1->GetMaximum()[axis]-box1->GetMinimum()[axis] ) / atof( argv[3] );
  RealType dz2 = ( box2->GetMaximum()[axis]-box2->GetMinimum()[axis] ) / atof( argv[4] );


  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( It.Get() == lobeLabel1 )
      {
      for ( int i = 0; i < atoi( argv[3] ); i++ )
        {
        RealType minimum = box1->GetMinimum()[axis] + static_cast<RealType>( i ) * dz1;
        RealType maximum = minimum + dz1;
        if ( It.GetIndex()[axis] >= minimum && It.GetIndex()[axis] <= maximum )
          {
          It.Set( i + 1 );
          }
        }
      }
    else if ( It.Get() == lobeLabel2 )
      {
      for ( int i = 0; i < atoi( argv[4] ); i++ )
        {
        RealType minimum = box2->GetMinimum()[axis] + static_cast<RealType>( i ) * dz2;
        RealType maximum = minimum + dz2;
        if ( It.GetIndex()[axis] >= minimum && It.GetIndex()[axis] <= maximum )
          {
          It.Set( atoi( argv[3] ) + i + 1 );
          }
        }
      }
    else
      {
      It.Set( 0 );
      }
    }


  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( reader->GetOutput() );
  writer->Update();


  return 0;
}

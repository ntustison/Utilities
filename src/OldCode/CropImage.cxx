#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkBoundingBox.h"
#include "itkExtractImageFilter.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc != 6 )
    {
    std::cout << "Usage: " << argv[0] << "image_filename segmented_image_filename size_per_dimension output_image_filename output_segmented_image_filename" << std::endl;
    exit( 1 );
    }

  typedef itk::BoundingBox<unsigned long, 
       ImageDimension, double> BoundingBoxType;  
  BoundingBoxType::Pointer box = BoundingBoxType::New();
  BoundingBoxType::PointsContainerPointer Points 
       = BoundingBoxType::PointsContainer::New();
  itk::Point<double, ImageDimension> point;
  
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  int idx = 0;
  itk::ImageRegionIteratorWithIndex<ImageType> It( 
    reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( It.Get() > 0 )
      {
      for (unsigned int i = 0; i < ImageDimension; i++ )
        {
        point[i] = static_cast<double>( It.GetIndex()[i] );
        }
      Points->InsertElement( idx++, point ); 
      }  
    }
  box->SetPoints( Points );
  box->ComputeBoundingBox();
  //box->GetBounds();
  
  ImageType::RegionType region;
  ImageType::RegionType::SizeType size;
  ImageType::RegionType::IndexType index;
  size.Fill( atoi( argv[3] ) );

  ImageType::RegionType::SizeType sz =
          reader->GetOutput()->GetLargestPossibleRegion().GetSize();
  
  itk::Point<double, ImageDimension> minPoint;
  itk::Point<double, ImageDimension> maxPoint;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    minPoint[i] = box->GetCenter()[i] - 0.5*atoi( argv[3] );
    maxPoint[i] = box->GetCenter()[i] + 0.5*atoi( argv[3] );
    
    if ( minPoint[i] < 0 )
      {
      minPoint[i] = 0;
      maxPoint[i] = static_cast<double>( atoi( argv[3] ) - 1 );
      }
    if ( maxPoint[i] >= sz[i] )
      {
      maxPoint[i] = static_cast<double>( sz[i] - 1 );
      minPoint[i] = maxPoint[i] - static_cast<double>( atoi( argv[3] ) - 1 ) ;
      }
    
    
    index[i] = static_cast<int>( minPoint[i] );
    }        
    
  box->SetMinimum( minPoint );
  box->SetMaximum( maxPoint );
  
  typedef itk::ExtractImageFilter<ImageType, ImageType> CropperType;
  CropperType::Pointer cropper = CropperType::New();

  region.SetSize( size );
  region.SetIndex( index );

  cropper->SetInput( reader->GetOutput() );
  cropper->SetExtractionRegion( region );
  cropper->Update();  
                  
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  
  writer->SetInput( cropper->GetOutput() );
  writer->SetFileName( argv[5] );                                          
  writer->Update();
  
  reader->SetFileName( argv[1] );
  reader->Update();
  cropper->SetInput( reader->GetOutput() );
  cropper->SetExtractionRegion( region );
  cropper->Update();  
  
  writer->SetInput( cropper->GetOutput() );
  writer->SetFileName( argv[4] );                                          
  writer->Update();
  
    
  return 0;
}

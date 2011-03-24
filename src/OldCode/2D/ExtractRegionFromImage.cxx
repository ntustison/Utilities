#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkExtractImageFilter.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 7 )
    {
    std::cout << "Usage: " << argv[0] << "inputImage outputImage minIndex[0] maxIndex[0] minIndex[1] maxIndex[1] [minIndex[2]] [maxIndex[2]] " << std::endl;
    exit( 1 );
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  ImageType::RegionType region;
  ImageType::RegionType::SizeType size;
  ImageType::RegionType::IndexType index;

  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    index[i] = atoi( argv[3+2*i] );
    size[i] = atoi( argv[3+2*i+1] ) - atoi( argv[3+2*i] ) + 1;
    }
  region.SetSize( size );
  region.SetIndex( index );

  typedef itk::ExtractImageFilter<ImageType, ImageType> CropperType;
  CropperType::Pointer cropper = CropperType::New();
  cropper->SetInput( reader->GetOutput() );
  cropper->SetExtractionRegion( region );
  cropper->Update();  
                  
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( cropper->GetOutput() );
  writer->SetFileName( argv[2] );                                          
  writer->Update();
  
  
    
  return 0;
}

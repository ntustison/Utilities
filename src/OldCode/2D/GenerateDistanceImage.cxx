#include <stdio.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc != 3 )
    {
    std::cout << "Usage: " << argv[0] << " inputImage outputImage";
    exit( 1 );
    }
  typedef itk::Image<PixelType, ImageDimension> ImageType;   

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::SignedMaurerDistanceMapImageFilter<ImageType, ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->SetSquaredDistance( false );
  filter->SetUseImageSpacing( false );
  filter->SetInsideIsPositive( true );
  filter->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();
  
  return 0;
}

#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkRescaleIntensityImageFilter.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc != 5 )
    {
    std::cout << "Usage: " << argv[0] << " image_filename output_filename min max" << std::endl;
    exit( 1 );
    }
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();
  
  typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescalerType;
  RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetInput( reader->GetOutput() );
  rescaler->SetOutputMinimum( atof( argv[3] ) );
  rescaler->SetOutputMaximum( atof( argv[4] ) );
  rescaler->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( rescaler->GetOutput() );
  writer->Update();

 return 0;
}

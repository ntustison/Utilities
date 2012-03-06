#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkRescaleIntensityImageFilter.h"

template <unsigned int ImageDimension>
int RescaleImageIntensity( int argc, char *argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescalerType;
  typename RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetInput( reader->GetOutput() );
  rescaler->SetOutputMinimum( atof( argv[4] ) );
  rescaler->SetOutputMaximum( atof( argv[5] ) );
  rescaler->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( rescaler->GetOutput() );
  writer->Update();

 return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 6 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage outputImage "
     << " min max" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     RescaleImageIntensity<2>( argc, argv );
     break;
   case 3:
     RescaleImageIntensity<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


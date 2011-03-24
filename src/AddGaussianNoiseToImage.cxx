#include "itkAdditiveGaussianNoiseImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

template <unsigned int ImageDimension>
int AddGaussianNoiseToImage( int argc, char *argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::AdditiveGaussianNoiseImageFilter<ImageType> NoiseFilterType;
  typename NoiseFilterType::Pointer noise = NoiseFilterType::New();
  noise->SetInput( reader->GetOutput() );
  if ( argc <= 4 )
    {
    noise->SetMean( 0 );
    noise->SetStandardDeviation( 1 );
    }
  else
    {
    noise->SetMean( atof( argv[4] ) );
    noise->SetStandardDeviation( atof( argv[5] ) );
    }
  noise->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( noise->GetOutput() );
  writer->SetFileName( argv[3] );
  writer->Update();

  return 0;
}


int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage outputImage " << std::flush;
    std::cout << "[mean] [standard deviation]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     AddGaussianNoiseToImage<2>( argc, argv );
     break;
   case 3:
     AddGaussianNoiseToImage<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

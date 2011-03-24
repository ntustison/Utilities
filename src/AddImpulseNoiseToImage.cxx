#include "itkImpulseNoiseImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

template <unsigned int ImageDimension>
int AddImpulseNoiseToImage( int argc, char *argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );  
  reader->Update();

  typedef itk::ImpulseNoiseImageFilter<ImageType> NoiseFilterType;
  typename NoiseFilterType::Pointer noise = NoiseFilterType::New();
  noise->SetInput( reader->GetOutput() );
  if ( argc <= 4 )
    {
    noise->SetProbability( 0.5 );
    } 
  else  
    {
    noise->SetProbability( atof( argv[4] ) );
    } 
  if ( argc >= 6 )
    {
    noise->SetMinimum( static_cast<PixelType>( atof( argv[5] ) ) );
    } 
  if ( argc >= 7 )
    {
    noise->SetMinimum( static_cast<PixelType>( atof( argv[6] ) ) );
    } 

  noise->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( noise->GetOutput() );
  writer->SetFileName( argv[2] );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage outputImage " << std::flush;
    std::cout << "[probability] [min] [max]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) ) 
   {
   case 2:
     AddImpulseNoiseToImage<2>( argc, argv );
     break;
   case 3:
     AddImpulseNoiseToImage<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

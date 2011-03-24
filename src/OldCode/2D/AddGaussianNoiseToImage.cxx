#include "itkAdditiveGaussianNoiseImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " inputImage outputImage " << std::flush;
    std::cout << "[mean] [standard deviation]" << std::endl;
    exit( 1 );
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );  
  reader->Update();

  typedef itk::AdditiveGaussianNoiseImageFilter<ImageType> NoiseFilterType;
  NoiseFilterType::Pointer noise = NoiseFilterType::New();
  noise->SetInput( reader->GetOutput() );
  if ( argc <= 3 )
    {
    noise->SetMean( 0 );
    noise->SetStandardDeviation( 1 );
    } 
  else  
    {
    noise->SetMean( atof( argv[3] ) );
    noise->SetStandardDeviation( atof( argv[4] ) );
    } 
  noise->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( noise->GetOutput() );
  writer->SetFileName( argv[2] );
  writer->Update();

  return 0;
}

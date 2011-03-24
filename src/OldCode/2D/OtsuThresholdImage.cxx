#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkOtsuMultipleThresholdsImageFilter.h"
#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " inputImage outputImage numberOfBins [numberOfThresholds]" << std::endl;
    exit( 1 );
    }
  
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::OtsuMultipleThresholdsImageFilter<ImageType, ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  if ( argc == 4 )
    {
    filter->SetNumberOfThresholds( 1 );
    }
  else
    {
    filter->SetNumberOfThresholds( atoi( argv[4] ) );
    } 
  filter->SetLabelOffset( 0 );
  filter->SetNumberOfHistogramBins( atoi( argv[3] ) ); 
  
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

  return 0;
}

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "global.h"

int main(int argc, char *argv[])
{
  if ( argc < 5 )
    {
      std::cout << "Usage: "<< argv[0] << " inputImage outputImage lowerThreshold upperThreshold [insideValue] [outsideValue]" << std::endl;
      exit( 1 );
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->SetLowerThreshold( static_cast<PixelType>( atof( argv[3] ) ) );
  filter->SetUpperThreshold( static_cast<PixelType>( atof( argv[4] ) ) );

  if ( argc >= 6 )
    {
    filter->SetInsideValue( static_cast<PixelType>( atof( argv[5] ) ) );
    }
  else
    {
    filter->SetInsideValue( static_cast<PixelType>( 1 ) );
    }
  if ( argc >= 7 )
    {
    filter->SetOutsideValue( static_cast<PixelType>( atof( argv[6] ) ) );
    }
  else
    {
    filter->SetOutsideValue( static_cast<PixelType>( 0 ) );
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

 return 0;
}

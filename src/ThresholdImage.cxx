#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryThresholdImageFilter.h"

template <unsigned int ImageDimension>
int ThresholdImage(int argc, char *argv[])
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<int, ImageDimension> LabelImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->SetLowerThreshold( static_cast<PixelType>( atof( argv[4] ) ) );
  filter->SetUpperThreshold( static_cast<PixelType>( atof( argv[5] ) ) );

  if ( argc >= 7 )
    {
    filter->SetInsideValue( static_cast<PixelType>( atof( argv[6] ) ) );
    }
  else
    {
    filter->SetInsideValue( static_cast<PixelType>( 1 ) );
    }
  if ( argc >= 8 )
    {
    filter->SetOutsideValue( static_cast<PixelType>( atof( argv[7] ) ) );
    }
  else
    {
    filter->SetOutsideValue( static_cast<PixelType>( 0 ) );
    }
//   filter->Print( std::cout, 5 );


  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

 return 0;
}

int main( int argc, char *argv[] )
  {
  if ( argc < 6 )
    {
    std::cout << "Usage: "<< argv[0] << " imageDimension inputImage outputImage "
      "lowerThreshold upperThreshold [insideValue] [outsideValue]" << std::endl;
    exit( 1 );
    }
 
switch( atoi( argv[1] ) )
    {

   case 2:
						ThresholdImage<2>( argc, argv );
						break;
				case 3:
						ThresholdImage<3>( argc, argv );
						break;
				case 4:
						ThresholdImage<3>( argc, argv );
						break;
				default:
						std::cerr << "Unsupported dimension" << std::endl;
						exit( EXIT_FAILURE );
				}
  }

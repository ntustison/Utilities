#include "itkAdaptImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRGBPixel.h"
#include "itkRGBAPixel.h"
#include "itkRedPixelAccessor.h"
#include "itkGreenPixelAccessor.h"
#include "itkBluePixelAccessor.h"

template <unsigned int ImageDimension>
int ExtractRGBComponent( int argc, char *argv[] )
{
  typedef unsigned int PixelType;
  typedef itk::RGBPixel<unsigned char> RGBPixelType;
//  typedef itk::RGBAPixel<unsigned char> RGBPixelType;

  typedef float RealType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<float, ImageDimension> RealImageType;
  typedef itk::Image<RGBPixelType, ImageDimension> RGBImageType;

  typedef itk::ImageFileReader<RGBImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::RedPixelAccessor<float> myRedAccessorType;
  typedef itk::GreenPixelAccessor<float> myGreenAccessorType;
  typedef itk::BluePixelAccessor<float> myBlueAccessorType;

  char component = ( argv[4] )[0];

  switch( component )
    {
    case 'r':
      {
      typename itk::AdaptImageFilter<RGBImageType, ImageType, myRedAccessorType>::Pointer adaptImageToRed =
        itk::AdaptImageFilter<RGBImageType, ImageType, myRedAccessorType>::New();
      adaptImageToRed->SetInput( reader->GetOutput() );
      adaptImageToRed->UpdateLargestPossibleRegion();
      adaptImageToRed->SetFunctor( adaptImageToRed->GetFunctor() );

      typedef itk::ImageFileWriter<ImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( adaptImageToRed->GetOutput() );
      writer->SetFileName( argv[3] );
      writer->Update();

      break;
      }
    case 'g':
      {
      typename itk::AdaptImageFilter<RGBImageType, ImageType, myGreenAccessorType>::Pointer adaptImageToGreen =
        itk::AdaptImageFilter<RGBImageType, ImageType, myGreenAccessorType>::New();
      adaptImageToGreen->SetInput( reader->GetOutput() );
      adaptImageToGreen->UpdateLargestPossibleRegion();
      adaptImageToGreen->SetFunctor( adaptImageToGreen->GetFunctor() );

      typedef itk::ImageFileWriter<ImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( adaptImageToGreen->GetOutput() );
      writer->SetFileName( argv[3] );
      writer->Update();

      break;
      }
    case 'b':
      {
      typename itk::AdaptImageFilter<RGBImageType, ImageType, myBlueAccessorType>::Pointer adaptImageToBlue =
        itk::AdaptImageFilter<RGBImageType, ImageType, myBlueAccessorType>::New();
      adaptImageToBlue->SetInput( reader->GetOutput() );
      adaptImageToBlue->UpdateLargestPossibleRegion();
      adaptImageToBlue->SetFunctor( adaptImageToBlue->GetFunctor() );

      typedef itk::ImageFileWriter<ImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( adaptImageToBlue->GetOutput() );
      writer->SetFileName( argv[3] );
      writer->Update();
      break;
      }
    default:
      {
      std::cout << "Unrecognized option." << std::endl;
      return EXIT_FAILURE;
      }
    }




  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage outputImage component(r,g,b)" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     ExtractRGBComponent<2>( argc, argv );
     break;
   case 3:
     ExtractRGBComponent<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


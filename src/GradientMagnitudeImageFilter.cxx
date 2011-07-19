#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkBoundedReciprocalImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"

template<unsigned int ImageDimension>
int GradientMagnitudeImageFilter( int argc, char * argv[] )
{
  typedef float PixelType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;


  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );

  typedef itk::GradientMagnitudeRecursiveGaussianImageFilter
    <ImageType, ImageType>  FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->SetSigma( ( argc > 4 ) ? atof( argv[4] ) : 1.0 );
  filter->Update();

  typedef itk::ImageFileWriter<ImageType>  WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cerr << argv[0] << " imageDimension inputImage outputImage "
      << "[sigma] " << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     GradientMagnitudeImageFilter<2>( argc, argv );
     break;
   case 3:
     GradientMagnitudeImageFilter<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


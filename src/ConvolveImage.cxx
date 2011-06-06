#include "itkConvolutionImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

template <unsigned int ImageDimension>
int ConvolveImage( int argc, char *argv[] )
{
  typedef float PixelType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typename ReaderType::Pointer kernelReader = ReaderType::New();
  kernelReader->SetFileName( argv[3] );
  kernelReader->Update();

  typedef itk::ConvolutionImageFilter<ImageType> FilterType;
  typename FilterType::Pointer convoluter = FilterType::New();
  convoluter->SetInput( reader->GetOutput() );
  convoluter->SetImageKernelInput( kernelReader->GetOutput() );
  convoluter->NormalizeOn();
  convoluter->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( convoluter->GetOutput() );
  writer->SetFileName( argv[4] );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension "
      << "inputImage kernelImage outputImage" << std::endl;
    return 0;
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     ConvolveImage<2>( argc, argv );
     break;
   case 3:
     ConvolveImage<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


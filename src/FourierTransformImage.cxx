#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkFFTShiftImageFilter.h"
#include "itkFFTComplexToComplexImageFilter.h"

#include "itkComplexToRealImageFilter.h"
#include "itkComplexToImaginaryImageFilter.h"

//#if !defined(USE_FFTWF)
////#error "This example only works when single precision FFTW is used
////Changing WorkPixeltype to double and changing this conditional to USE_FFTWD
////will also work.
//#endif

template <unsigned int ImageDimension>
int FourierTransform( int argc, char * argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  typename ReaderType::Pointer inputRealReader = ReaderType::New();
  typename ReaderType::Pointer inputImagReader = ReaderType::New();

  inputRealReader->SetFileName( argv[2] );
  inputImagReader->SetFileName( argv[3] );

  typename ImageType::Pointer realImage = ImageType::New();
  realImage = NULL;
  typename ImageType::Pointer imagImage = ImageType::New();
  imagImage = NULL;

  realImage = inputRealReader->GetOutput();

  try
    {
    inputImagReader->Update();
    imagImage = inputImagReader->GetOutput();
    }
  catch(...)
    {
    imagImage->SetDirection( realImage->GetDirection() );
    imagImage->SetSpacing( realImage->GetSpacing() );
    imagImage->SetOrigin( realImage->GetOrigin() );
    imagImage->SetRegions( realImage->GetLargestPossibleRegion() );
    imagImage->Allocate();
    imagImage->FillBuffer( 0 );
    std::cout << "Setting imaginary image to all zeros.";
    }

  typedef std::complex<float> ComplexPixelType;
  typedef itk::Image<ComplexPixelType, ImageDimension> ComplexImageType;

  typename ComplexImageType::Pointer complexImage = ComplexImageType::New();
  complexImage->CopyInformation( realImage );
  complexImage->SetRegions( realImage->GetLargestPossibleRegion() );
  complexImage->Allocate();

  itk::ImageRegionIteratorWithIndex<ComplexImageType> It( complexImage,
    complexImage->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    PixelType real = realImage->GetPixel( It.GetIndex() );
    PixelType imag = imagImage->GetPixel( It.GetIndex() );

    ComplexPixelType complex( real, imag );

    It.Set( complex );
    }

  typedef itk::FFTComplexToComplexImageFilter<ComplexImageType> FFTFilterType;
  typename FFTFilterType::Pointer fftoutput = FFTFilterType::New();
  typedef itk::FFTComplexToComplexImageFilter
    <ComplexImageType> invFFTFilterType;
  typename invFFTFilterType::Pointer invfftoutput = invFFTFilterType::New();

  if( argc > 7 && atoi( argv[7] ) != 0 )
    {
    typedef itk::FFTShiftImageFilter<ComplexImageType, ComplexImageType>
      ShifterType;
    typename ShifterType::Pointer shifter = ShifterType::New();
    shifter->SetInput( complexImage );
    shifter->SetInverse( ( argc > 9 ) ? atoi( argv[9] ) : 1 );
    shifter->Update();
    complexImage = shifter->GetOutput();
    complexImage->DisconnectPipeline();
    }

  if( argc > 6 && atoi( argv[6] ) == 1 )
    {
    invfftoutput->SetTransformDirection( FFTFilterType::INVERSE );
    invfftoutput->SetInput( complexImage ); // compute inverse FFT
    invfftoutput->Update();
    complexImage = invfftoutput->GetOutput();
    complexImage->DisconnectPipeline();
    }
  else
    {
    fftoutput->SetTransformDirection( FFTFilterType::DIRECT );
    fftoutput->SetInput( complexImage ); // compute forward FFT
    fftoutput->Update();
    complexImage = fftoutput->GetOutput();
    complexImage->DisconnectPipeline();
    }

  if( argc > 8 && atoi( argv[8] ) != 0 )
    {
    typedef itk::FFTShiftImageFilter<ComplexImageType, ComplexImageType>
      ShifterType;
    typename ShifterType::Pointer shifter = ShifterType::New();
    shifter->SetInput( complexImage );
    shifter->SetInverse( ( argc > 9 ) ? atoi( argv[9] ) : 1 );
    shifter->Update();
    complexImage = shifter->GetOutput();
    complexImage->DisconnectPipeline();
    }

  typedef itk::ComplexToRealImageFilter<ComplexImageType, ImageType> RealerFilterType;
  typename RealerFilterType::Pointer realer = RealerFilterType::New();
  realer->SetInput( complexImage );
  realer->Update();

  typedef itk::ComplexToImaginaryImageFilter<ComplexImageType, ImageType> ImaginerFilterType;
  typename ImaginerFilterType::Pointer imaginer = ImaginerFilterType::New();
  imaginer->SetInput( complexImage );
  imaginer->Update();

// Write the output
  typedef itk::ImageFileWriter<ImageType> WriterType;

  typename WriterType::Pointer realWriter = WriterType::New();
  realWriter->SetFileName( argv[4] );
  realWriter->SetInput( realer->GetOutput() );
  realWriter->Update();
  typename WriterType::Pointer imagWriter = WriterType::New();
  imagWriter->SetFileName( argv[5] );
  imagWriter->SetInput( imaginer->GetOutput() );
  imagWriter->Update();

  return EXIT_SUCCESS;

}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " imageDimension inputRealImage inputImaginaryImage "
      << "outputRealImage outputImaginaryImage [direction] [prefftshift] "
      << " [postfftshift] [fftShiftInverse]" << std::endl;
    std::cout << "   direction=0 -> forward" << std::endl;
    std::cout << "   direction=1 -> inverse" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     FourierTransform<2>( argc, argv );
     break;
   case 3:
     FourierTransform<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}



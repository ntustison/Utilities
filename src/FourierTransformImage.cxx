#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkFFTShiftImageFilter.h"
#include "itkFFTWComplexToComplexImageFilter.h"
#include "itkRealAndImaginaryToComplexImageFilter.h"
#include "itkComplexToImaginaryImageFilter.h"
#include "itkComplexToRealImageFilter.h"

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

  typedef itk::RealAndImaginaryToComplexImageFilter<PixelType, PixelType,
    PixelType, ImageDimension> ComplexFilterType;
  typename ComplexFilterType::Pointer complexer = ComplexFilterType::New();
  complexer->SetInput1( realImage );
  complexer->SetInput2( imagImage );
  complexer->Update();

  typedef itk::FFTWComplexToComplexImageFilter
    <PixelType, ImageDimension> FFTFilterType;
  typename FFTFilterType::Pointer fftoutput = FFTFilterType::New();
  typedef itk::FFTWComplexToComplexImageFilter
    <PixelType, ImageDimension> invFFTFilterType;
  typename invFFTFilterType::Pointer invfftoutput = invFFTFilterType::New();

  typedef typename FFTFilterType::OutputImageType ComplexImageType;
  typename ComplexImageType::Pointer complexImage = ComplexImageType::New();
  complexImage = complexer->GetOutput();

  if( argc > 7 && atoi( argv[7] ) != 0 )
    {
    typedef itk::FFTShiftImageFilter<ComplexImageType, ComplexImageType>
      ShifterType;
    typename ShifterType::Pointer shifter = ShifterType::New();
    shifter->SetInput( complexImage );
    shifter->SetInverse( ( argc > 9 ) ? atoi( argv[9] ) : 1 );
    shifter->Update();
    complexImage = shifter->GetOutput();
    }

  if( argc > 6 && atoi( argv[6] ) == 1 )
    {
    invfftoutput->SetTransformDirection( FFTFilterType::INVERSE );
    invfftoutput->SetInput( complexImage ); // compute inverse FFT
    invfftoutput->Update();
    complexImage = invfftoutput->GetOutput();
    }
  else
    {
    fftoutput->SetTransformDirection( FFTFilterType::DIRECT );
    fftoutput->SetInput( complexImage ); // compute forward FFT
    fftoutput->Update();
    complexImage = fftoutput->GetOutput();
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









//#include "itkImage.h"
//#include "itkImageFileReader.h"
//#include "itkImageFileWriter.h"
//
//#include "itkImageRegionIteratorWithIndex.h"
//
//#include "itkFFTShiftImageFilter.h"
//#include "itkFFTWComplexConjugateToRealImageFilter.h"
//#include "itkFFTWRealToComplexConjugateImageFilter.h"
//
//#include "itkVnlFFTComplexConjugateToRealImageFilter.h"
//#include "itkVnlFFTRealToComplexConjugateImageFilter.h"
//
//#include "itkBinaryMagnitudeImageFilter.h"
//#include "itkComplexToImaginaryImageFilter.h"
//#include "itkComplexToRealImageFilter.h"
//
//#define USE_FFTW 1
//
//
//template <unsigned int ImageDimension>
//int FFT( int argc, char *argv[] )
//{
//  typedef float RealType;
//  typedef RealType PixelType;
//  typedef itk::Image<PixelType, ImageDimension> ImageType;
//  typedef itk::Image<float, ImageDimension> RealImageType;
//
//  typedef itk::ImageFileReader<ImageType> ReaderType;
//  typename ReaderType::Pointer reader = ReaderType::New();
//  reader->SetFileName( argv[2] );
//
//
//  typename RealImageType::SpacingType spacing;
//  for ( unsigned int i = 0; i < ImageDimension; i++ )
//    {
//    spacing[i] = 1.0 / reader->GetOutput()->GetSpacing()[i];
//    }
//
//
//#ifdef USE_FFTW
//  typedef itk::FFTWRealToComplexConjugateImageFilter<RealType, ImageDimension> FFTFilterType;
//#else
//  typedef itk::VnlFFTRealToComplexConjugateImageFilter<RealType, ImageDimension> FFTFilterType;
//#endif
//  typename FFTFilterType::Pointer fftFilter = FFTFilterType::New();
//  fftFilter->SetInput( reader->GetOutput() );
//  fftFilter->Update();
//  fftFilter->GetOutput()->SetSpacing( spacing );
//
//  typedef typename FFTFilterType::OutputImageType ComplexImageType;
//  typename ComplexImageType::Pointer fft = ComplexImageType::New();
//
//#ifdef USE_FFTW
//  fft->SetSpacing( spacing );
//  fft->SetOrigin( fftFilter->GetOutput()->GetOrigin() );
//  fft->SetRegions( reader->GetOutput()->GetRequestedRegion() );
//  fft->Allocate();
//
//  itk::ImageRegionIteratorWithIndex<ComplexImageType> ItF( fft,
//    fft->GetLargestPossibleRegion() );
//
//  for ( ItF.GoToBegin(); !ItF.IsAtEnd(); ++ItF )
//    {
//    if ( ItF.GetIndex()[0] <= 0.5 * fft->GetLargestPossibleRegion().GetSize()[0] )
//      {
//      ItF.Set( fftFilter->GetOutput()->GetPixel( ItF.GetIndex() ) );
//      }
//    else
//      {
//      typename ComplexImageType::IndexType index;
//      for ( unsigned int i = 0; i < ImageDimension; i++ )
//        {
//        index[i] = fft->GetLargestPossibleRegion().GetSize()[i] - ItF.GetIndex()[i];
//        }
//      typename ComplexImageType::PixelType pixel
//        = fftFilter->GetOutput()->GetPixel( index );
//      ItF.Set( std::complex<RealType>( pixel.real(), -pixel.imag() ) );
//      }
//    }
//
//#else
//  fft = fftFilter->GetOutput();
//#endif
//
//  typedef itk::FFTShiftImageFilter<ComplexImageType, ComplexImageType> ShifterType;
//  typename ShifterType::Pointer shifter = ShifterType::New();
//  shifter->SetInput( fft );
//  shifter->SetInverse( 1 );
//  shifter->Update();
//
//  /**
//   * Debug: Write magnitude image of imaginary and real parts.
//   */
//  typedef itk::ComplexToImaginaryImageFilter<ComplexImageType, RealImageType> ImaginerType;
//  typename ImaginerType::Pointer imaginer = ImaginerType::New();
//  imaginer->SetInput( shifter->GetOutput() );
//  imaginer->Update();
//
//  typedef itk::ComplexToRealImageFilter<ComplexImageType, RealImageType> RealerType;
//  typename RealerType::Pointer realer = RealerType::New();
//  realer->SetInput( shifter->GetOutput() );
//  realer->Update();
//
//  typedef itk::BinaryMagnitudeImageFilter<RealImageType, RealImageType, RealImageType> BinarierType;
//  typename BinarierType::Pointer binarier = BinarierType::New();
//  binarier->SetInput1( imaginer->GetOutput() );
//  binarier->SetInput2( realer->GetOutput() );
//  binarier->Update();
//
//  typedef itk::ImageFileWriter<ImageType> WriterType;
//  typename WriterType::Pointer writer = WriterType::New();
//  writer->SetFileName( argv[3] );
//  if ( atoi( argv[4] ) == 0 )
//    {
//    writer->SetInput( realer->GetOutput() );
//    }
//  else if ( atoi( argv[4] ) == 1 )
//    {
//    writer->SetInput( imaginer->GetOutput() );
//    }
//  else
//    {
//    writer->SetInput( binarier->GetOutput() );
//    }
//  writer->Update();
//
//  return 0;
//}
//
//int main( int argc, char *argv[] )
//{
//  if ( argc < 4 )
//    {
//    std::cout << argv[0] << " imageDimension inputImage outputImage which" << std::endl;
//    std::cout << "   0. Real part " << std::endl;
//    std::cout << "   1. Imaginary part " << std::endl;
//    std::cout << "   2. Magnitude " << std::endl;
//    exit( 1 );
//    }
//
//  switch( atoi( argv[1] ) )
//   {
//   case 2:
//     FFT<2>( argc, argv );
//     break;
//   case 3:
//     FFT<3>( argc, argv );
//     break;
//   default:
//      std::cerr << "Unsupported dimension" << std::endl;
//      exit( EXIT_FAILURE );
//   }
//}

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageRegionIteratorWithIndex.h"

#include "itkFFTShiftImageFilter.h"
#include "itkFFTWComplexConjugateToRealImageFilter.h"
#include "itkFFTWRealToComplexConjugateImageFilter.h"

#include "itkVnlFFTComplexConjugateToRealImageFilter.h"
#include "itkVnlFFTRealToComplexConjugateImageFilter.h"

#include "itkBinaryMagnitudeImageFilter.h"
#include "itkComplexToImaginaryImageFilter.h"
#include "itkComplexToRealImageFilter.h"

#include "global.h"

#define USE_FFTW

int main( int argc, char *argv[] )
{
  if ( argc != 4 )
    {
    std::cout << "Usage: " << argv[0] << " intputImage outputImage which" << std::endl;
    std::cout << "   0. Real part " << std::endl;
    std::cout << "   1. Imaginary part " << std::endl;
    std::cout << "   2. Magnitude " << std::endl;
    exit( 1 );
    }

  typedef float RealType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<float, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );  


  RealImageType::SpacingType spacing;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    { 
    spacing[i] = 1.0 / reader->GetOutput()->GetSpacing()[i];
    }


#ifdef USE_FFTW   
  typedef itk::FFTWRealToComplexConjugateImageFilter<RealType, ImageDimension> FFTFilterType;
#else
  typedef itk::VnlFFTRealToComplexConjugateImageFilter<RealType, ImageDimension> FFTFilterType;
#endif
  FFTFilterType::Pointer fftFilter = FFTFilterType::New();
  fftFilter->SetInput( reader->GetOutput() );
  fftFilter->Update(); 
  fftFilter->GetOutput()->SetSpacing( spacing );

  typedef FFTFilterType::OutputImageType ComplexImageType;
  ComplexImageType::Pointer fft = ComplexImageType::New();

#ifdef USE_FFTW   
  fft->SetSpacing( spacing );
  fft->SetOrigin( fftFilter->GetOutput()->GetOrigin() );
  fft->SetRegions( reader->GetOutput()->GetRequestedRegion() );  
  fft->Allocate();

  itk::ImageRegionIteratorWithIndex<ComplexImageType> ItF( fft,
    fft->GetLargestPossibleRegion() ); 

  for ( ItF.GoToBegin(); !ItF.IsAtEnd(); ++ItF )
    {
    if ( ItF.GetIndex()[0] <= 0.5 * fft->GetLargestPossibleRegion().GetSize()[0] )
      {
      ItF.Set( fftFilter->GetOutput()->GetPixel( ItF.GetIndex() ) );  
      } 
    else
      {
      ComplexImageType::IndexType index;
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        index[i] = fft->GetLargestPossibleRegion().GetSize()[i] - ItF.GetIndex()[i];
        } 
      ComplexImageType::PixelType pixel
        = fftFilter->GetOutput()->GetPixel( index );
      ItF.Set( std::complex<RealType>( pixel.real(), -pixel.imag() ) ); 
      }
    }

#else
  fft = fftFilter->GetOutput();
#endif 

  typedef itk::FFTShiftImageFilter<ComplexImageType, ComplexImageType> ShifterType;
  ShifterType::Pointer shifter = ShifterType::New();
  shifter->SetInput( fft );
  shifter->SetInverse( 1 );
  shifter->Update();
  
  /** 
   * Debug: Write magnitude image of imaginary and real parts.
   */ 
  typedef itk::ComplexToImaginaryImageFilter<ComplexImageType, RealImageType> ImaginerType;
  ImaginerType::Pointer imaginer = ImaginerType::New();
  imaginer->SetInput( shifter->GetOutput() );
  imaginer->Update();

  typedef itk::ComplexToRealImageFilter<ComplexImageType, RealImageType> RealerType;
  RealerType::Pointer realer = RealerType::New();
  realer->SetInput( shifter->GetOutput() );
  realer->Update();

  typedef itk::BinaryMagnitudeImageFilter<RealImageType, RealImageType, RealImageType> BinarierType;
  BinarierType::Pointer binarier = BinarierType::New();
  binarier->SetInput1( imaginer->GetOutput() );
  binarier->SetInput2( realer->GetOutput() );
  binarier->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  if ( atoi( argv[3] ) == 0 )
    {
    writer->SetInput( realer->GetOutput() );
    } 
  else if ( atoi( argv[3] ) == 1 )
    {
    writer->SetInput( imaginer->GetOutput() );
    } 
  else
    {
    writer->SetInput( binarier->GetOutput() );
    } 
  writer->Update();

  return 0;
}

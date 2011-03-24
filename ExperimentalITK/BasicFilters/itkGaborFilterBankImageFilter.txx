#ifndef __itkGaborFilterBankImageFilter_txx
#define __itkGaborFilterBankImageFilter_txx

#include "itkGaborFilterBankImageFilter.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkEuler3DTransform.h"
#include "itkExtractImageFilter.h"
#include "itkFFTShiftImageFilter.h"
#include "itkFFTWComplexConjugateToRealImageFilter.h"
#include "itkFFTWRealToComplexConjugateImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkVnlFFTComplexConjugateToRealImageFilter.h"
#include "itkVnlFFTRealToComplexConjugateImageFilter.h"

#include "itkBinaryMagnitudeImageFilter.h"
#include "itkComplexToImaginaryImageFilter.h"
#include "itkComplexToRealImageFilter.h"
#include "itkImageFileWriter.h"

#include "vnl/vnl_math.h"

#define USE_FFTW

namespace itk
{

template <class TInputImage, class TOutputImage>
GaborFilterBankImageFilter<TInputImage, TOutputImage>
::GaborFilterBankImageFilter()
{
}

template <class TInputImage, class TOutputImage>
GaborFilterBankImageFilter<TInputImage, TOutputImage>
::~GaborFilterBankImageFilter()
{  
}

template <class TInputImage, class TOutputImage>
void
GaborFilterBankImageFilter<TInputImage, TOutputImage>
::GenerateData()
{

  this->AllocateOutputs();

  this->GetOutput( NumericTraits<RealType>::NonpositiveMin() );

  /**
   * Note regarding tagging geometry:  Assume that the tagging planes are perpendicular
   * to the imaging planes.  We set the x-y plane of the coordinate system so that it
   * is parallel to the imaging planes.  (theta, psi, phi) are the Euler angles around
   * the x, y, and z axes, respectively.
   */

  typename InputImageType::SpacingType spacing;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    { 
    spacing[i] = 1.0 / this->GetInput()->GetSpacing()[i];
    }

  /**
   * Generate the fourier transform of the input image
   */
#ifdef USE_FFTW   
  typedef FFTWRealToComplexConjugateImageFilter
    <RealType, ImageDimension> FFTFilterType;
#else
  typedef VnlFFTRealToComplexConjugateImageFilter
    <RealType, ImageDimension> FFTFilterType;
#endif
  typename FFTFilterType::Pointer fftFilter = FFTFilterType::New();
  fftFilter->SetInput( this->GetInput() );
  fftFilter->Update(); 
  fftFilter->GetOutput()->SetSpacing( spacing );

  typedef typename FFTFilterType::OutputImageType ComplexImageType;
  typename ComplexImageType::Pointer fft = ComplexImageType::New();

#ifdef USE_FFTW   
  fft->SetSpacing( spacing );
  fft->SetOrigin( fftFilter->GetOutput()->GetOrigin() );
  fft->SetRegions( this->GetInput()->GetRequestedRegion() ); 
  fft->SetDirection( this->GetInput()->GetDirection() ); 
  fft->Allocate();

  ImageRegionIteratorWithIndex<ComplexImageType> ItF( fft,
    fft->GetRequestedRegion() ); 

  for ( ItF.GoToBegin(); !ItF.IsAtEnd(); ++ItF )
    {
    if ( ItF.GetIndex()[0] <= 0.5 * fft->GetRequestedRegion().GetSize()[0] )
      {
      ItF.Set( fftFilter->GetOutput()->GetPixel( ItF.GetIndex() ) );  
      } 
    else
      {
      typename ComplexImageType::IndexType index;
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        index[i] = fft->GetRequestedRegion().GetSize()[i] - ItF.GetIndex()[i];
        } 
      typename ComplexImageType::PixelType pixel
        = fftFilter->GetOutput()->GetPixel( index );
      ItF.Set( std::complex<RealType>( pixel.real(), -pixel.imag() ) ); 
      }
    }

#else
  fft = fftFilter->GetOutput();
#endif 

  typedef FFTShiftImageFilter<ComplexImageType, ComplexImageType> ShifterType;
  typename ShifterType::Pointer shifter = ShifterType::New();
  shifter->SetInput( fft );
  shifter->SetInverse( 1 );
  shifter->Update();

  
  /** 
   * Debug: Write magnitude image of imaginary and real parts.
   */ 
  /*
  typedef ComplexToImaginaryImageFilter<ComplexImageType, RealImageType> ImaginerType;
  typename ImaginerType::Pointer imaginer = ImaginerType::New();
  imaginer->SetInput( shifter->GetOutput() );
  imaginer->Update();

  typedef ComplexToImaginaryImageFilter<ComplexImageType, RealImageType> RealerType;
  typename RealerType::Pointer realer = RealerType::New();
  realer->SetInput( shifter->GetOutput() );
  realer->Update();

  typedef BinaryMagnitudeImageFilter<RealImageType, RealImageType, RealImageType> BinarierType;
  typename BinarierType::Pointer binarier = BinarierType::New();
  binarier->SetInput1( imaginer->GetOutput() );
  binarier->SetInput2( realer->GetOutput() );
  binarier->Update();

  typedef ImageFileWriter<RealImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( binarier->GetOutput() );
  writer->SetFileName( "FFTimage.hdr" );
  writer->Update();
  */  

  unsigned int iteration = 0;

  RealType deltaGaborSpacing = vnl_math_max( static_cast<RealType>( 1.0 ), 
    this->m_GaborSpacingMaximum - this->m_GaborSpacingMinimum );
  RealType deltaPhiSpacing = vnl_math_max( static_cast<RealType>( 1.0 ), 
    this->m_RotationAngleMaximum[2] - this->m_RotationAngleMinimum[2] );
  RealType deltaPsiSpacing = vnl_math_max( static_cast<RealType>( 1.0 ), 
    this->m_RotationAngleMaximum[1] - this->m_RotationAngleMinimum[1] );
  RealType deltaThetaSpacing = vnl_math_max( static_cast<RealType>( 1.0 ), 
    this->m_RotationAngleMaximum[0] - this->m_RotationAngleMinimum[0] );

  for ( RealType gaborSpacing = this->m_GaborSpacingMinimum; 
        gaborSpacing <= this->m_GaborSpacingMaximum; 
        gaborSpacing += deltaGaborSpacing 
          / static_cast<RealType>( this->m_NumberOfGaborSpacingSteps - 1 ) 
      )
    {
    for ( RealType phi = this->m_RotationAngleMinimum[2]; 
          phi <= this->m_RotationAngleMaximum[2]; 
          phi += deltaPhiSpacing 
            / static_cast<RealType>( this->m_NumberOfRotationAngleSteps[2] - 1 )  
        )
      {
      for ( RealType psi = this->m_RotationAngleMinimum[1]; 
            psi <= this->m_RotationAngleMaximum[1]; 
            psi += deltaPsiSpacing 
              / static_cast<RealType>( this->m_NumberOfRotationAngleSteps[1] - 1 )  
          )
        {
        for ( RealType theta = this->m_RotationAngleMinimum[0]; 
              theta <= this->m_RotationAngleMaximum[0]; 
              theta += deltaThetaSpacing 
                / static_cast<RealType>( this->m_NumberOfRotationAngleSteps[0] - 1 )  
            )
          {
          std::cout << "Iteration " << iteration << ": " << gaborSpacing << ", " << phi << ", " << psi << ", " << theta << std::endl; 
  
          /**
           * Calculate the initial parameters
           */
  
          ArrayType fundamentalFrequency;
          fundamentalFrequency[0] = gaborSpacing; 
          fundamentalFrequency[1] = 0; 
          fundamentalFrequency[2] = 0; 

          ArrayType sigma;
          sigma[0] = 1.0 / gaborSpacing;
          sigma[1] = sigma[2] = 2.0 * sigma[0]; 
        
          /**
           * Generate the fourier transform of the gabor filter
           */
          typename RealImageType::Pointer fourierGaborImage 
            = RealImageType::New();
          fourierGaborImage->SetOrigin( this->GetInput()->GetOrigin() );
          fourierGaborImage->SetSpacing( spacing );
          fourierGaborImage->SetRegions( this->GetInput()->GetRequestedRegion().GetSize() );
          fourierGaborImage->Allocate();
          fourierGaborImage->FillBuffer( 0.0 );
             
          for ( RealType d = -1.0; d <= 1.01; d+=2.0 )
            {  
            typename FourierTransformGaborImageSourceType::Pointer 
              fourierGabor = FourierTransformGaborImageSourceType::New();
            fourierGabor->SetOrigin( this->GetInput()->GetOrigin() );
            fourierGabor->SetSpacing( spacing );
            fourierGabor->SetSize( this->GetInput()->GetRequestedRegion().GetSize() );
          
            typename FourierTransformGaborImageSourceType::ArrayType gaussianSigma;
            typename FourierTransformGaborImageSourceType::ArrayType mean;
          
            for ( unsigned int i = 0; i < ImageDimension; i++ )
              {
              gaussianSigma[i] = 1.0 / ( 2.0 * sigma[i] * vnl_math::pi );
              mean[i] = this->GetInput()->GetOrigin()[i] 
                + 0.5 * ( this->GetInput()->GetRequestedRegion().GetSize()[i] - 1 ) 
                      * ( spacing[i] ) + d * fundamentalFrequency[i];
              } 
    
            fourierGabor->SetSigma( gaussianSigma );
            fourierGabor->SetMean( mean );
            fourierGabor->SetScale( 1.0 );
            fourierGabor->SetNormalized( false );
            fourierGabor->Update();
            
            ImageRegionIterator<RealImageType> It( fourierGabor->GetOutput(),
              fourierGabor->GetOutput()->GetRequestedRegion() );
            ImageRegionIterator<RealImageType> ItF( fourierGaborImage,
              fourierGaborImage->GetRequestedRegion() );
            for ( It.GoToBegin(), ItF.GoToBegin(); !It.IsAtEnd(); ++It, ++ItF )  
              {
              ItF.Set( ItF.Get() + It.Get() );
              }  
            }       
       
          /**
           * Rotate the fourier transformed gabor filter
           */
          typedef Euler3DTransform<RealType> TransformType;
          typename TransformType::Pointer transform = TransformType::New();
        
          transform->SetRotation( theta, psi, phi );
          typename TransformType::OutputVectorType translation;
          translation.Fill( 0 );
          transform->SetTranslation( translation );
          typename TransformType::InputPointType center;
          for ( unsigned int i = 0; i < ImageDimension; i++ )
            {
            center[i] = this->GetInput()->GetOrigin()[i] + 0.5
              * ( this->GetInput()->GetRequestedRegion().GetSize()[i] - 1 ) 
              * ( spacing[i] );
            }
          transform->SetCenter( center );
          
          typedef LinearInterpolateImageFunction
            <RealImageType, RealType> LinearInterpolatorType;
          typename LinearInterpolatorType::Pointer interpolator 
            = LinearInterpolatorType::New();
          interpolator->SetInputImage( fourierGaborImage );
        
          typedef ResampleImageFilter
            <RealImageType, RealImageType, RealType> ResamplerType;
          typename ResamplerType::Pointer resampler = ResamplerType::New();
          resampler->SetDefaultPixelValue( NumericTraits<RealType>::Zero );
          resampler->SetTransform( transform );
          resampler->SetInterpolator( interpolator );
          resampler->SetInput( fourierGaborImage );
          resampler->SetOutputSpacing( fourierGaborImage->GetSpacing() );
          resampler->SetOutputOrigin( fourierGaborImage->GetOrigin() );
          resampler->SetSize( 
            fourierGaborImage->GetRequestedRegion().GetSize() );
          resampler->Update();
  
          /**
           * Multiply in frequency space
           */
          typename ComplexImageType::Pointer filteredImage 
            = ComplexImageType::New();
          filteredImage->SetRegions( 
            resampler->GetOutput()->GetRequestedRegion() );
          filteredImage->SetOrigin( resampler->GetOutput()->GetOrigin() );
          filteredImage->SetSpacing( resampler->GetOutput()->GetSpacing() );
          filteredImage->Allocate();
        
          ImageRegionIterator<ComplexImageType> ItS( shifter->GetOutput(), 
            shifter->GetOutput()->GetRequestedRegion() );
          ImageRegionIterator<RealImageType> ItG( resampler->GetOutput(), 
            resampler->GetOutput()->GetRequestedRegion() );
          ImageRegionIterator<ComplexImageType> It( filteredImage, 
            filteredImage->GetRequestedRegion() );
        
          ItS.GoToBegin();
          ItG.GoToBegin();
          It.GoToBegin();
          while ( !It.IsAtEnd()  )
            {
            It.Set( ItS.Get() * std::complex<RealType>( 
              -ItG.Get(), NumericTraits<RealType>::Zero ) );     
        
            ++ItS;
            ++ItG;  
            ++It;
            }  
        
          /**
           * Generate the fourier transform of the filtered image
           */
        
          typename ShifterType::Pointer inverseShifter = ShifterType::New();
          inverseShifter->SetInput( filteredImage );
          inverseShifter->SetInverse( 1 );
          inverseShifter->Update();

          typename ComplexImageType::Pointer ifft = ComplexImageType::New();

#ifdef USE_FFTW   
          typedef FFTWComplexConjugateToRealImageFilter
            <RealType, ImageDimension> InverseFFTFilterType;

          typedef ExtractImageFilter
            <ComplexImageType, ComplexImageType> ExtracterType;
          typename ExtracterType::Pointer extracter = ExtracterType::New();
          typename ComplexImageType::RegionType region;
          typename ComplexImageType::RegionType::SizeType size;
          size = inverseShifter->GetOutput()->GetRequestedRegion().GetSize();
          size[0] = static_cast<unsigned int>( 0.5 * size[0] ) + 1;
          region.SetIndex( 
            inverseShifter->GetOutput()->GetRequestedRegion().GetIndex() );
          region.SetSize( size );
          extracter->SetInput( inverseShifter->GetOutput() );
          extracter->SetExtractionRegion( region );  
          extracter->Update();

          ifft = extracter->GetOutput();
#else
          typedef VnlFFTComplexConjugateToRealImageFilter
            <RealType, ImageDimension> InverseFFTFilterType;
          ifft = inverseShifter->GetOutput();
#endif
          typename InverseFFTFilterType::Pointer ifftFilter 
            = InverseFFTFilterType::New();
          ifftFilter->SetInput( ifft );
          ifftFilter->Update(); 
          ifftFilter->GetOutput()->SetSpacing( this->GetInput()->GetSpacing() );

          ImageRegionIterator<RealImageType> ItI( ifftFilter->GetOutput(), 
            ifftFilter->GetOutput()->GetRequestedRegion() );
          ImageRegionIterator<RealImageType> ItR( this->GetOutput(), 
            this->GetOutput()->GetRequestedRegion() );
     
          for ( ItI.GoToBegin(), ItR.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItR )
            {
            ItR.Set( vnl_math_max( ItI.Get(), ItR.Get() ) );
            }  
  
          iteration++;
          } 
        }
      }
    }      

}

/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TOutputImage>
void
GaborFilterBankImageFilter<TInputImage, TOutputImage>
::PrintSelf(
  std::ostream& os, 
  Indent indent) const
{
}



}  //end namespace itk

#endif

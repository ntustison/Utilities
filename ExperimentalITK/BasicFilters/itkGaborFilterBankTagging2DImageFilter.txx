#ifndef __itkGaborFilterBankTagging2DImageFilter_txx
#define __itkGaborFilterBankTagging2DImageFilter_txx

#include "itkGaborFilterBankTagging2DImageFilter.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkEuler2DTransform.h"
#include "itkExtractImageFilter.h"
#include "itkFFTShiftImageFilter.h"
#include "itkFFTWComplexConjugateToRealImageFilter.h"
#include "itkFFTWRealToComplexConjugateImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkOtsuMultipleThresholdsImageFilter.h"
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
GaborFilterBankTagging2DImageFilter<TInputImage, TOutputImage>
::GaborFilterBankTagging2DImageFilter()
{
  if ( ImageDimension != 2 )
    {
    itkExceptionMacro( << "Image dimension must be equal to 2." );
    }
  this->m_MaskImage = NULL;
  this->m_ThresholdPercentage = 0.9;
  this->m_UseAutomaticThresholding = true;
}

template <class TInputImage, class TOutputImage>
GaborFilterBankTagging2DImageFilter<TInputImage, TOutputImage>
::~GaborFilterBankTagging2DImageFilter()
{  
}

template <class TInputImage, class TOutputImage>
void
GaborFilterBankTagging2DImageFilter<TInputImage, TOutputImage>
::GenerateData()
{

  this->m_MaximumResponseImage = RealImageType::New();
  this->m_MaximumResponseImage->SetOrigin( this->GetInput()->GetOrigin() );
  this->m_MaximumResponseImage->SetSpacing( this->GetInput()->GetSpacing() );
  this->m_MaximumResponseImage->SetRegions( this->GetInput()->GetRequestedRegion() );
  this->m_MaximumResponseImage->Allocate();
  this->m_MaximumResponseImage->FillBuffer( NumericTraits<RealType>::NonpositiveMin() );
 
  if ( !this->m_MaskImage )
    { 
    this->m_MaskImage = LabelImageType::New();
    this->m_MaskImage->SetOrigin( this->GetInput()->GetOrigin() );
    this->m_MaskImage->SetSpacing( this->GetInput()->GetSpacing() );
    this->m_MaskImage->SetRegions( this->GetInput()->GetRequestedRegion() );
    this->m_MaskImage->Allocate();
    this->m_MaskImage->FillBuffer( 1 );
    } 

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
  typedef FFTWRealToComplexConjugateImageFilter<RealType, ImageDimension> FFTFilterType;
#else
  typedef VnlFFTRealToComplexConjugateImageFilter<RealType, ImageDimension> FFTFilterType;
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
  fft->Allocate();

  ImageRegionIteratorWithIndex<ComplexImageType> ItF( fft,
    fft->GetLargestPossibleRegion() ); 

  for ( ItF.GoToBegin(); !ItF.IsAtEnd(); ++ItF )
    {
    if ( ItF.GetIndex()[0] <= 0.5 * fft->GetLargestPossibleRegion().GetSize()[0] )
      {
      ItF.Set( fftFilter->GetOutput()->GetPixel( ItF.GetIndex() ) );  
      } 
    else
      {
      typename ComplexImageType::IndexType index;
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        index[i] = fft->GetLargestPossibleRegion().GetSize()[i] - ItF.GetIndex()[i];
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

  typedef ComplexToRealImageFilter<ComplexImageType, RealImageType> RealerType;
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
  writer->SetFileName( "FFTimage.nii" );
  writer->Update();
  */

  unsigned int iteration = 0;

  for ( RealType tagSpacing = this->m_TagSpacingMinimum; 
        tagSpacing <= this->m_TagSpacingMaximum; 
        tagSpacing += ( this->m_TagSpacingMaximum - this->m_TagSpacingMinimum ) 
                      / static_cast<RealType>( this->m_NumberOfTagSpacingSteps - 1 ) 
      )
    {
    for ( RealType psi = this->m_RotationAngleMinimum[1]; 
          psi <= this->m_RotationAngleMaximum[1]; 
          psi += ( this->m_RotationAngleMaximum[1] - this->m_RotationAngleMinimum[1] ) 
                 / static_cast<RealType>( this->m_NumberOfRotationAngleSteps[1] - 1 )  
        )
      {
      for ( RealType theta = this->m_RotationAngleMinimum[0]; 
            theta <= this->m_RotationAngleMaximum[0]; 
            theta += ( this->m_RotationAngleMaximum[0] - this->m_RotationAngleMinimum[0] ) 
                   / static_cast<RealType>( this->m_NumberOfRotationAngleSteps[0] - 1 )  
          )
        {
//        std::cout << "Iteration " << iteration << ": " << tagSpacing  << ", " << psi << ", " << theta << std::endl; 

        /**
         * Calculate the initial parameters
         */

        ArrayType fundamentalFrequency;
        fundamentalFrequency[0] = tagSpacing; 
        fundamentalFrequency[1] = 0; 
/*
        fundamentalFrequency[1] = tagSpacing / 
          vcl_sqrt( ( 1 + vnl_math_sqr( vcl_tan( phi ) ) )
                  * ( 1 + vnl_math_sqr( vcl_tan( psi ) ) ) );
        fundamentalFrequency[0] = fundamentalFrequency[1] * vcl_tan( phi );
        fundamentalFrequency[2] = vcl_tan( psi )
          * vcl_sqrt( vnl_math_sqr( fundamentalFrequency[0] ) 
                    + vnl_math_sqr( fundamentalFrequency[1] ) );
*/
        ArrayType sigma;
        sigma[0] = 0.4 / tagSpacing;
        sigma[1] = 7.5 * sigma[0]; 
            
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
              + 0.5 * ( this->GetInput()->GetLargestPossibleRegion().GetSize()[i] - 1 ) 
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
        typedef Euler2DTransform<RealType> TransformType;
        typename TransformType::Pointer transform = TransformType::New();
      
        transform->SetRotation( theta );
        typename TransformType::OutputVectorType translation;
        translation.Fill( 0 );
        transform->SetTranslation( translation );
        typename TransformType::InputPointType center;
        for ( unsigned int i = 0; i < ImageDimension; i++ )
          {
          center[i] = this->GetInput()->GetOrigin()[i] 
            + 0.5 * ( this->GetInput()->GetLargestPossibleRegion().GetSize()[i] - 1 ) 
                  * ( spacing[i] );
          }
        transform->SetCenter( center );
        
        typedef LinearInterpolateImageFunction<RealImageType, RealType> LinearInterpolatorType;
        typename LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();
        interpolator->SetInputImage( fourierGaborImage );
      
        typedef ResampleImageFilter<RealImageType, RealImageType, RealType> ResamplerType;
        typename ResamplerType::Pointer resampler = ResamplerType::New();
        resampler->SetDefaultPixelValue( 0 );
        resampler->SetTransform( transform );
        resampler->SetInterpolator( interpolator );
        resampler->SetInput( fourierGaborImage );
        resampler->SetOutputSpacing( fourierGaborImage->GetSpacing() );
        resampler->SetOutputOrigin( fourierGaborImage->GetOrigin() );
        resampler->SetSize( fourierGaborImage->GetLargestPossibleRegion().GetSize() );
        resampler->Update();

        /**
         * Multiply in frequency space
         */
        typename ComplexImageType::Pointer filteredImage = ComplexImageType::New();
        filteredImage->SetRegions( resampler->GetOutput()->GetLargestPossibleRegion() );
        filteredImage->SetOrigin( resampler->GetOutput()->GetOrigin() );
        filteredImage->SetSpacing( resampler->GetOutput()->GetSpacing() );
        filteredImage->Allocate();
      
        ImageRegionIterator<ComplexImageType> ItS( shifter->GetOutput(), 
          shifter->GetOutput()->GetLargestPossibleRegion() );
        ImageRegionIterator<RealImageType> ItG( resampler->GetOutput(), 
          resampler->GetOutput()->GetLargestPossibleRegion() );
        ImageRegionIterator<ComplexImageType> It( filteredImage, 
          filteredImage->GetLargestPossibleRegion() );
      
        ItS.GoToBegin();
        ItG.GoToBegin();
        It.GoToBegin();
        while ( !It.IsAtEnd()  )
          {
          It.Set( ItS.Get() * std::complex<RealType>( -ItG.Get(), 0 ) );     
      
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
        typedef FFTWComplexConjugateToRealImageFilter<RealType, ImageDimension> InverseFFTFilterType;

        typedef ExtractImageFilter<ComplexImageType, ComplexImageType> ExtracterType;
        typename ExtracterType::Pointer extracter = ExtracterType::New();
        typename ComplexImageType::RegionType region;
        typename ComplexImageType::RegionType::SizeType size;
        size = inverseShifter->GetOutput()->GetLargestPossibleRegion().GetSize();
        size[0] = static_cast<unsigned int>( 0.5 * size[0] ) + 1;
        region.SetIndex( inverseShifter->GetOutput()->GetLargestPossibleRegion().GetIndex() );
        region.SetSize( size );
        extracter->SetInput( inverseShifter->GetOutput() );
        extracter->SetExtractionRegion( region );  
        extracter->Update();

        ifft = extracter->GetOutput();
#else
        typedef VnlFFTComplexConjugateToRealImageFilter<RealType, ImageDimension> InverseFFTFilterType;
        ifft = inverseShifter->GetOutput();
#endif
        typename InverseFFTFilterType::Pointer ifftFilter = InverseFFTFilterType::New();
        ifftFilter->SetInput( ifft );
        ifftFilter->Update(); 
        ifftFilter->GetOutput()->SetSpacing( this->GetInput()->GetSpacing() );
   
//        {
//        itk::OStringStream buf;
//        buf << iteration;
//        std::string filename = std::string( "ifft" )
//          + buf.str() + std::string( ".nii.gz" );       
//
//        typedef ImageFileWriter<RealImageType> WriterType;
//        typename WriterType::Pointer writer = WriterType::New();
//        writer->SetInput( ifftFilter->GetOutput() );
//        writer->SetFileName( filename.c_str() );
//        writer->Update();
//        } 

        ImageRegionIterator<RealImageType> ItI( ifftFilter->GetOutput(), 
          ifftFilter->GetOutput()->GetLargestPossibleRegion() );
        ImageRegionIterator<RealImageType> ItR( this->m_MaximumResponseImage, 
          this->m_MaximumResponseImage->GetLargestPossibleRegion() );
   
        for ( ItI.GoToBegin(), ItR.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItR )
          {
          ItR.Set( vnl_math_max( ItI.Get(), ItR.Get() ) );
          }  

        iteration++;
        } 
      }
    }

  typename OutputImageType::Pointer output = OutputImageType::New();
  output->SetOrigin( this->GetInput()->GetOrigin() );
  output->SetSpacing( this->GetInput()->GetSpacing() );
  output->SetRegions( this->GetInput()->GetRequestedRegion() );
  output->Allocate();
  output->FillBuffer( 0 );

  if ( this->m_UseAutomaticThresholding )
    {
    typedef OtsuMultipleThresholdsImageFilter<RealImageType, OutputImageType> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput( this->m_MaximumResponseImage );
    filter->SetNumberOfThresholds( 1 );
    filter->SetLabelOffset( 0 );
    filter->SetNumberOfHistogramBins( 100 );
    filter->Update();

    ImageRegionIterator<LabelImageType> ItM( this->m_MaskImage, 
      this->m_MaskImage->GetLargestPossibleRegion() );
    ImageRegionIterator<OutputImageType> ItO( output, 
      output->GetLargestPossibleRegion() );
    ImageRegionIterator<OutputImageType> ItF( filter->GetOutput(), 
      filter->GetOutput()->GetLargestPossibleRegion() );

    for ( ItF.GoToBegin(), ItM.GoToBegin(), ItO.GoToBegin(); !ItM.IsAtEnd(); ++ItM, ++ItF, ++ItO )
      {
      if ( ItF.Get() > 0 )
        {
        ItO.Set( ItM.Get() );
        }   
      }   

    } 
  else
    {
  
    /**
     * Threshold filtered image and incorporate into final output
     */  
    typedef LabelStatisticsImageFilter<RealImageType, LabelImageType> HistogramGeneratorType;
    typename HistogramGeneratorType::Pointer stats = HistogramGeneratorType::New();
    stats->SetInput( this->m_MaximumResponseImage );
    stats->SetLabelInput( this->m_MaskImage );
    stats->UseHistogramsOff();
    stats->Update();
  
    unsigned int numberOfBins = 100;
    
    for ( unsigned n = 1; n < stats->GetNumberOfLabels(); n++ )
      {
      RealType delta = ( static_cast<RealType>( stats->GetMaximum( n ) ) 
        - static_cast<RealType>( stats->GetMinimum( n ) ) )
        / static_cast<RealType>( numberOfBins );
      RealType lowerBound = static_cast<RealType>( stats->GetMinimum( n ) ) - delta;
      RealType upperBound = static_cast<RealType>( stats->GetMaximum( n ) ) + delta;
    
      stats->UseHistogramsOn();
      stats->SetHistogramParameters( numberOfBins, lowerBound, upperBound );
      stats->Update();
  
      typedef typename HistogramGeneratorType::HistogramType  HistogramType;
      const HistogramType *histogram = stats->GetHistogram( n );
      
  /*
      unsigned long cumulativeSum[numberOfBins];
      cumulativeSum[0] = histogram->GetFrequency( 0, 0 );
      RealType percentage = static_cast<RealType>( cumulativeSum[0] ) 
        / static_cast<RealType>( histogram->GetTotalFrequency() );        
    
      unsigned int i = 1;
    
      while ( percentage < this->m_ThresholdPercentage && i < numberOfBins )
        {
        cumulativeSum[i] = cumulativeSum[i-1] + histogram->GetFrequency( i, 0 );
        percentage = static_cast<RealType>( cumulativeSum[i] ) 
          / static_cast<RealType>( histogram->GetTotalFrequency() );
        i++;        
        }
    
      RealType m2 = static_cast<RealType>( histogram->GetMeasurement( i, 0 ) );
      RealType m1 = static_cast<RealType>( histogram->GetMeasurement( i-1, 0 ) );
      RealType p2 = static_cast<RealType>( cumulativeSum[i] ) 
        / static_cast<RealType>( histogram->GetTotalFrequency() );        
      RealType p1 = static_cast<RealType>( cumulativeSum[i-1] ) 
        / static_cast<RealType>( histogram->GetTotalFrequency() );        
      RealType lowerThreshold = m2 
        - ( p2 - this->m_ThresholdPercentage ) * ( m2 - m1 ) / ( p2 - p1 );            
  */
      RealType lowerThreshold = histogram->Quantile( 0, this->m_ThresholdPercentage );
      RealType upperThreshold = stats->GetMaximum( n );
  
      if ( lowerThreshold > upperThreshold )
        {
        lowerThreshold = upperThreshold;
        }   
  
      typedef BinaryThresholdImageFilter<RealImageType, RealImageType> ThresholderType;
      typename ThresholderType::Pointer thresholder = ThresholderType::New();
      thresholder->SetInput( this->m_MaximumResponseImage );
      thresholder->SetInsideValue( 1 );
      thresholder->SetOutsideValue( 0 );
      thresholder->SetLowerThreshold( lowerThreshold ); 
      thresholder->SetUpperThreshold( upperThreshold ); 
      thresholder->Update();
  
      ImageRegionIterator<LabelImageType> ItM( this->m_MaskImage, 
        this->m_MaskImage->GetLargestPossibleRegion() );
      ImageRegionIterator<OutputImageType> ItO( output, 
        output->GetLargestPossibleRegion() );
      ImageRegionIterator<RealImageType> ItR( this->m_MaximumResponseImage, 
        this->m_MaximumResponseImage->GetLargestPossibleRegion() );
      ImageRegionIterator<RealImageType> ItT( thresholder->GetOutput(), 
        thresholder->GetOutput()->GetLargestPossibleRegion() );
       
      ItM.GoToBegin();
      ItO.GoToBegin();
      ItT.GoToBegin();
      while ( !ItM.IsAtEnd() )
        {
        if ( ItT.Get() > 0 && ItM.Get() == n )
          {
          ItO.Set( static_cast<typename OutputImageType::PixelType>( n ) ); 
          }  
        ++ItM; 
        ++ItO;
        ++ItT;
        }
      }  
    }     
  this->SetNthOutput( 0, output );
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TOutputImage>
void
GaborFilterBankTagging2DImageFilter<TInputImage, TOutputImage>
::PrintSelf(
  std::ostream& os, 
  Indent indent) const
{
}



}  //end namespace itk

#endif

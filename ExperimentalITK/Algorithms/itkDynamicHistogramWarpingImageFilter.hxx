/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkDynamicHistogramWarpingImageFilter_hxx
#define __itkDynamicHistogramWarpingImageFilter_hxx

#include "itkDynamicHistogramWarpingImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkNumericTraits.h"
#include "itkVariableSizeMatrix.h"

namespace itk
{
/**
 *
 */
template<class TInputImage, class TOutputImage, class THistogramMeasurement>
DynamicHistogramWarpingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>
::DynamicHistogramWarpingImageFilter() :
  m_NumberOfHistogramLevels( 256 ),
  m_MaximumSourceBinCompressionSize( 30 ),
  m_MaximumReferenceBinCompressionSize( 30 ),
  m_ThresholdAtMeanIntensity( true ),
  m_SourceIntensityThreshold( 0 ),
  m_ReferenceIntensityThreshold( 0 ),
  m_OutputIntensityThreshold( 0 ),
  m_SourceMinValue( 0 ),
  m_SourceMaxValue( 0 ),
  m_SourceMeanValue( 0 ),
  m_ReferenceMinValue( 0 ),
  m_ReferenceMaxValue( 0 ),
  m_ReferenceMeanValue( 0 ),
  m_OutputMinValue( 0 ),
  m_OutputMaxValue( 0 ),
  m_OutputMeanValue( 0 )
{
  this->SetNumberOfRequiredInputs( 2 );

  this->m_SourceHistogram = HistogramType::New();
  this->m_ReferenceHistogram = HistogramType::New();
  this->m_OutputHistogram = HistogramType::New();
}

/*
 *
 */
template<class TInputImage, class TOutputImage, class THistogramMeasurement>
void
DynamicHistogramWarpingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Number of histogram levels: ";
  os << this->m_NumberOfHistogramLevels << std::endl;
  os << indent << "Maximum reference bin compression size: ";
  os << this->m_MaximumReferenceBinCompressionSize << std::endl;
  os << indent << "Maximum source bin compression size: ";
  os << this->m_MaximumSourceBinCompressionSize << std::endl;
  os << indent << "Source histogram: ";
  os << this->m_SourceHistogram.GetPointer() << std::endl;
  os << indent << "Reference histogram: ";
  os << this->m_ReferenceHistogram.GetPointer() << std::endl;

  os << indent << "Threshold at mean intensity: ";
  os << this->m_ThresholdAtMeanIntensity << std::endl;

  os << indent << "Source intensity threshold: ";
  os << this->m_SourceIntensityThreshold << std::endl;
  os << indent << "Reference intensity threshold: ";
  os << this->m_ReferenceIntensityThreshold << std::endl;
  os << indent << "Output intensity threshold: ";
  os << this->m_OutputIntensityThreshold << std::endl;

  os << indent << "Source mean value: ";
  os << this->m_SourceMeanValue << std::endl;
  os << indent << "Reference mean value: ";
  os << this->m_ReferenceMeanValue << std::endl;
  os << indent << "Output mean value: ";
  os << this->m_OutputMeanValue << std::endl;

  os << indent << "Source min value: ";
  os << this->m_SourceMinValue << std::endl;
  os << indent << "Reference min value: ";
  os << this->m_ReferenceMinValue << std::endl;
  os << indent << "Output min value: ";
  os << this->m_OutputMinValue << std::endl;

  os << indent << "Source max value: ";
  os << this->m_SourceMaxValue << std::endl;
  os << indent << "Reference max value: ";
  os << this->m_ReferenceMaxValue << std::endl;
  os << indent << "Output max value: ";
  os << this->m_OutputMaxValue << std::endl;
}

/*
 *
 */
template<class TInputImage, class TOutputImage, class THistogramMeasurement>
void
DynamicHistogramWarpingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>
::SetReferenceImage( const InputImageType *reference )
{
  this->ProcessObject::SetNthInput( 1,
    const_cast<InputImageType *>( reference ) );
}

/*
 *
 */
template<class TInputImage, class TOutputImage, class THistogramMeasurement>
const typename DynamicHistogramWarpingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>
::InputImageType *
DynamicHistogramWarpingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>
::GetReferenceImage()
{
  if( this->GetNumberOfInputs() < 2 )
    {
    return NULL;
    }

  return dynamic_cast<TInputImage *>(
    this->ProcessObject::GetInput(1) );
}

/*
 * This filter requires all of the input images to be
 * in the buffer.
 */
template<class TInputImage, class TOutputImage, class THistogramMeasurement>
void
DynamicHistogramWarpingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>
::GenerateInputRequestedRegion()
{
  this->Superclass::GenerateInputRequestedRegion();

  for( unsigned int idx = 0; idx < this->GetNumberOfInputs(); ++idx )
    {
    if( this->GetInput( idx ) )
      {
      InputImagePointer image =
        const_cast<InputImageType *>( this->GetInput( idx ) );
      image->SetRequestedRegionToLargestPossibleRegion();
      }
    }
}

/**
 *
 */
template<class TInputImage, class TOutputImage, class THistogramMeasurement>
void
DynamicHistogramWarpingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>
::BeforeThreadedGenerateData()
{
  InputImageConstPointer source = this->GetSourceImage();
  InputImageConstPointer reference = this->GetReferenceImage();

  this->ComputeMinMaxMean( source, this->m_SourceMinValue,
    this->m_SourceMaxValue, this->m_SourceMeanValue );
  this->ComputeMinMaxMean( reference, this->m_ReferenceMinValue,
    this->m_ReferenceMaxValue, this->m_ReferenceMeanValue );


  if( this->m_ThresholdAtMeanIntensity )
    {
    this->m_SourceIntensityThreshold =
      static_cast<InputPixelType>( this->m_SourceMeanValue );
    this->m_ReferenceIntensityThreshold =
      static_cast<InputPixelType>( this->m_ReferenceMeanValue );
    }
  else
    {
    this->m_SourceIntensityThreshold =
      static_cast<InputPixelType>( this->m_SourceMinValue );
    this->m_ReferenceIntensityThreshold =
      static_cast<InputPixelType>( this->m_ReferenceMinValue );
    }

  this->ConstructHistogram( source, this->m_SourceHistogram,
    this->m_SourceIntensityThreshold, this->m_SourceMaxValue );

  this->ConstructHistogram( reference, this->m_ReferenceHistogram,
    this->m_ReferenceIntensityThreshold, this->m_ReferenceMaxValue );

  this->CalculateOptimalHistogramMapping( this->m_ReferenceHistogram,
    this->m_SourceHistogram );
}

template<class TInputImage, class TOutputImage, class THistogramMeasurement>
void
DynamicHistogramWarpingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>
::AfterThreadedGenerateData()
{
  OutputImagePointer output = this->GetOutput();

  this->ComputeMinMaxMean( output, this->m_OutputMinValue,
    this->m_OutputMaxValue, this->m_OutputMeanValue );

  if ( this->m_ThresholdAtMeanIntensity )
    {
    this->m_OutputIntensityThreshold =
      static_cast<OutputPixelType>( this->m_OutputMeanValue );
    }
  else
    {
    this->m_OutputIntensityThreshold =
      static_cast<OutputPixelType>( this->m_OutputMinValue );
    }

  this->ConstructHistogram( output, this->m_OutputHistogram,
    this->m_OutputIntensityThreshold, this->m_OutputMaxValue );
}

template<class TInputImage, class TOutputImage, class THistogramMeasurement>
void
DynamicHistogramWarpingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>
::ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread,
  int threadId )
{
  // Transform the source image and write to output.
  ImageRegionConstIterator<InputImageType> inIter(
    this->GetInput(), outputRegionForThread );

  ImageRegionIterator<OutputImageType> outIter(
    this->GetOutput(), outputRegionForThread );

  // support progress methods/callbacks
  SizeValueType updateVisits = 0;
  SizeValueType totalPixels = 0;
  if( threadId == 0 )
    {
    totalPixels = outputRegionForThread.GetNumberOfPixels();
    updateVisits = totalPixels / 10;
    if( updateVisits < 1 )
      {
      updateVisits = 1;
      }
    }

  RealType srcValue = 0.0;
  RealType mappedValue = 0.0;

  for( unsigned int i = 0; !outIter.IsAtEnd(); ++inIter, ++outIter, i++ )
    {
    if( threadId == 0 && !( i % updateVisits ) )
      {
      this->UpdateProgress( static_cast<float>( i ) /
        static_cast<float>( totalPixels ) );
      }

    srcValue = static_cast<RealType>( inIter.Get() );

    if( srcValue <= this->m_WarpedSourceHistogramProfile[0] )
      {
      mappedValue = this->m_WarpedReferenceHistogramProfile[0];
      }
    else if( srcValue >= this->m_WarpedSourceHistogramProfile[
      this->m_WarpedSourceHistogramProfile.size()-1] )
      {
      mappedValue = this->m_WarpedSourceHistogramProfile[
        this->m_WarpedSourceHistogramProfile.size()-1];
      }
    else
      {
      typename std::deque<RealType>::const_iterator it;
      for( it = this->m_WarpedSourceHistogramProfile.begin() + 1; it !=
        this->m_WarpedSourceHistogramProfile.end(); ++it )
        {
        if( *it > srcValue )
          {
          SizeValueType index = it -
            this->m_WarpedSourceHistogramProfile.begin();

          RealType proportion = 1.0;
          if( this->m_WarpedSourceHistogramProfile[index] !=
            this->m_WarpedSourceHistogramProfile[index+1] )
            {
            proportion = ( srcValue -
              this->m_WarpedSourceHistogramProfile[index] ) /
              ( this->m_WarpedSourceHistogramProfile[index+1] -
              this->m_WarpedSourceHistogramProfile[index] );
            }
          mappedValue = this->m_WarpedReferenceHistogramProfile[index] +
            proportion * ( this->m_WarpedReferenceHistogramProfile[index+1] -
            this->m_WarpedReferenceHistogramProfile[index] );
          break;
          }
        }
      }

    outIter.Set( static_cast< OutputPixelType >( mappedValue ) );
    }
}

/**
 * Compute min, max and mean of an image.
 */
template<class TInputImage, class TOutputImage, class THistogramMeasurement>
void
DynamicHistogramWarpingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>
::ComputeMinMaxMean( const InputImageType *image,
  THistogramMeasurement & minValue, THistogramMeasurement & maxValue,
  THistogramMeasurement & meanValue )
{
  typedef ImageRegionConstIterator< InputImageType > ConstIterator;
  ConstIterator iter( image, image->GetBufferedRegion() );

  RealType sum = 0.0;

  SizeValueType count = 0;

  minValue = static_cast< THistogramMeasurement >( iter.Get() );
  maxValue = minValue;

  while( !iter.IsAtEnd() )
    {
    const THistogramMeasurement value =
      static_cast< THistogramMeasurement >( iter.Get() );
    sum += static_cast< double >( value );

    if( value < minValue )
      {
      minValue = value;
      }
    if( value > maxValue )
      {
      maxValue = value;
      }
    ++iter;
    ++count;
    }

  meanValue = static_cast<THistogramMeasurement>( sum
    / static_cast<RealType>( count ) );
}

/**
 * Construct a histogram from an image.
 */
template<class TInputImage, class TOutputImage, class THistogramMeasurement>
void
DynamicHistogramWarpingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>
::ConstructHistogram( const InputImageType *image, HistogramType  *histogram,
  const THistogramMeasurement minValue, const THistogramMeasurement maxValue )
{
  // allocate memory for the histogram
  typename HistogramType::SizeType size;
  typename HistogramType::MeasurementVectorType lowerBound;
  typename HistogramType::MeasurementVectorType upperBound;

  size.SetSize( 1 );
  lowerBound.SetSize( 1 );
  upperBound.SetSize( 1 );
  histogram->SetMeasurementVectorSize(1);

  size[0] = this->m_NumberOfHistogramLevels;
  lowerBound.Fill( minValue );
  upperBound.Fill( maxValue );

  //Initialize with equally spaced bins.
  histogram->Initialize( size, lowerBound, upperBound );
  histogram->SetToZero();

  typename HistogramType::MeasurementVectorType measurement;
  measurement.SetSize( 1 );

  typedef typename HistogramType::MeasurementType MeasurementType;
  measurement[0] = NumericTraits< MeasurementType >::Zero;

  // put each image pixel into the histogram
  ImageRegionConstIterator<InputImageType>
    iter( image, image->GetBufferedRegion() );

  iter.GoToBegin();
  while( !iter.IsAtEnd() )
    {
    InputPixelType value = iter.Get();

    if( static_cast<RealType>( value ) >= minValue &&
      static_cast<RealType>( value ) <= maxValue )
      {
      // add sample to histogram
      measurement[0] = value;
      histogram->IncreaseFrequencyOfMeasurement( measurement, 1 );
      }
    ++iter;
    }
}

/**
 * Construct cumulative histogram.
 */
template<class TInputImage, class TOutputImage, class THistogramMeasurement>
void
DynamicHistogramWarpingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>
::ConstructCumulativeHistogram( const HistogramType *histogram,
  HistogramType *cumulativeHistogram )
{
  typename HistogramType::SizeType size = histogram->GetSize();

  typename HistogramType::MeasurementVectorType lowerBound;
  typename HistogramType::MeasurementVectorType upperBound;

  lowerBound.SetSize( 1 );
  upperBound.SetSize( 1 );

  lowerBound.Fill( histogram->GetBinMin( 0, 0 ) );
  upperBound.Fill( histogram->GetBinMax( 0, size[0]-1 ) );

  cumulativeHistogram->SetMeasurementVectorSize( 1 );
  cumulativeHistogram->Initialize( size, lowerBound, upperBound );
  cumulativeHistogram->SetToZero();

  typename HistogramType::ConstIterator It = histogram->Begin();
  typename HistogramType::Iterator ItC = cumulativeHistogram->Begin();

  typename HistogramType::AbsoluteFrequencyType sum = 0.0;
  while( It != histogram->End() )
    {
    sum += It.GetFrequency();
    ItC.SetFrequency( sum );
    ++ItC;
    ++It;
    }
}

/**
 * Generate cost table.
 */
template<class TInputImage, class TOutputImage, class THistogramMeasurement>
void
DynamicHistogramWarpingImageFilter<TInputImage, TOutputImage, THistogramMeasurement>
::CalculateOptimalHistogramMapping( const HistogramType *referenceHistogram,
  const HistogramType *sourceHistogram )
{
  // Create cumulative histograms

  typename HistogramType::Pointer referenceCumulativeHistogram =
    HistogramType::New();
  this->ConstructCumulativeHistogram( referenceHistogram,
    referenceCumulativeHistogram );

  typename HistogramType::Pointer sourceCumulativeHistogram =
    HistogramType::New();
  this->ConstructCumulativeHistogram( sourceHistogram,
    sourceCumulativeHistogram );

  VariableSizeMatrix<RealType> D;
  D.SetSize( referenceCumulativeHistogram->GetSize( 0 ),
    sourceCumulativeHistogram->GetSize( 0 ) );

  // Formulate the cost table
  for( unsigned int m = 0; m < D.Rows(); m++ )
    {
    for( unsigned int n = 0; n < D.Cols(); n++ )
      {
      if( m == 0 && n == 0 )
        {
        D(m, n) = 0.0;
        }
      else if( m == 0 || n == 0 )
        {
        D(m, n) = NumericTraits<RealType>::max();
        }
      else
        {
        typename HistogramType::IndexType referenceIndex;
        referenceIndex.SetSize( 1 );
        referenceIndex[0] = m;
        typename HistogramType::IndexType sourceIndex;
        sourceIndex.SetSize( 1 );
        sourceIndex[0] = n;

        RealType case1 = D(m - 1, n - 1) +
          vnl_math_abs( referenceHistogram->GetFrequency( referenceIndex ) -
          sourceHistogram->GetFrequency( sourceIndex ) );
        RealType case2 = NumericTraits<RealType>::max();
        for( unsigned int k = 2;
          k <= this->m_MaximumReferenceBinCompressionSize; k++ )
          {
          typename HistogramType::IndexType kIndex;
          kIndex.SetSize( 1 );
          kIndex[0] = referenceIndex[0] - k;
          if( kIndex[0] < 0 )
            {
            break;
            }
          case2 = vnl_math_min( case2, D(m - k, n - 1) +
            vnl_math_abs(
              ( referenceCumulativeHistogram->GetFrequency( referenceIndex ) -
              referenceCumulativeHistogram->GetFrequency( kIndex ) ) -
              sourceHistogram->GetFrequency( sourceIndex )
              )
            );
          }
        RealType case3 = NumericTraits<RealType>::max();
        for( unsigned int l = 2;
          l <= this->m_MaximumSourceBinCompressionSize; l++ )
          {
          typename HistogramType::IndexType lIndex;
          lIndex.SetSize( 1 );
          lIndex[0] = sourceIndex[0] - l;
          if( lIndex[0] < 0 )
            {
            break;
            }
          case3 = vnl_math_min( case3, D(m - 1, n - l) +
            vnl_math_abs(
              referenceHistogram->GetFrequency( referenceIndex ) -
              ( sourceCumulativeHistogram->GetFrequency( sourceIndex ) -
              sourceCumulativeHistogram->GetFrequency( lIndex ) )
              )
            );
          }
        D(m, n) = vnl_math_min( case1, vnl_math_min( case2, case3 ) );
        }
      }
    }

  // Find the optimal histogram mapping by tracing back through the
  // cost table.
  unsigned int m = D.Rows() - 1;
  unsigned int n = D.Cols() - 1;
  while( m > 0 && n > 0 )
    {
    RealType referenceBinCenter = 0.5 * ( referenceHistogram->GetBinMin( 0, m ) +
      referenceHistogram->GetBinMax( 0, m ) );
    this->m_WarpedReferenceHistogramProfile.push_front( referenceBinCenter );

    RealType sourceBinCenter = 0.5 * ( sourceHistogram->GetBinMin( 0, n ) +
      sourceHistogram->GetBinMax( 0, n ) );
    this->m_WarpedSourceHistogramProfile.push_front( sourceBinCenter );

    unsigned int minM = m;
    unsigned int minN = n - 1;
    if( D(m - 1, n - 1) < D(minM, minN) )
      {
      minM = m - 1;
      minN = n - 1;
      }
    else if( D(m - 1, n) < D(minM, minN) )
      {
      minM = m - 1;
      minN = n;
      }
    m = minM;
    n = minN;
    }
}

} // end namespace itk

#endif

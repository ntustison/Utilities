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
#ifndef __itkTextureFeaturesImageFilter_hxx
#define __itkTextureFeaturesImageFilter_hxx

#include "itkTextureFeaturesImageFilter.h"

#include "itkHistogramToTextureFeaturesFilter.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkProgressReporter.h"

#include <vector>

namespace itk
{
namespace Statistics
{
template<class TInputImage, class TOutputImage>
TextureFeaturesImageFilter<TInputImage, TOutputImage>
::TextureFeaturesImageFilter()
{
  // Set the offset directions to their defaults: half of all the possible
  // directions 1 pixel away. (The other half is included by symmetry.)
  // We use a neighborhood iterator to calculate the appropriate offsets.
  typedef Neighborhood<typename InputImageType::PixelType, ImageDimension> NeighborhoodType;
  NeighborhoodType neighborhood;
  neighborhood.SetRadius( 1 );

  // select all "previous" neighbors that are face+edge+vertex
  // connected to the current pixel. do not include the center pixel.
  unsigned int        centerIndex = neighborhood.GetCenterNeighborhoodIndex();
  OffsetType          offset;

  this->m_Offsets.clear();
  for ( unsigned int d = 0; d < centerIndex; d++ )
    {
    offset = neighborhood.GetOffset( d );
    this->m_Offsets.push_back( offset );
    }

  this->m_NumberOfBinsPerAxis = 64;

  this->m_InsidePixelValue = 1;

  this->m_NeighborhoodRadius.Fill( 10 );
}

template<class TInputImage, class TOutputImage>
void
TextureFeaturesImageFilter<TInputImage, TOutputImage>
::SetMaskImage( const MaskImageType *mask )
{
  this->SetNthInput( 1, const_cast<MaskImageType *>( mask ) );
}

template<class TInputImage, class TOutputImage>
const typename TextureFeaturesImageFilter<TInputImage, TOutputImage>::MaskImageType *
TextureFeaturesImageFilter<TInputImage, TOutputImage>
::GetMaskImage() const
{
  const MaskImageType *maskImage = dynamic_cast<const MaskImageType *>( this->ProcessObject::GetInput( 1 ) );

  return maskImage;
}

template< class TInputImage, class TOutputImage>
void
TextureFeaturesImageFilter<TInputImage, TOutputImage>
::SetOffset( const OffsetType &offset )
{
  OffsetVectorType offsetVector;

  offsetVector.push_back( offset );
  this->SetOffsets( offsetVector );
}

template< class TInputImage, class TOutputImage>
void
TextureFeaturesImageFilter<TInputImage, TOutputImage>
::SetPixelValueMinMax( InputPixelType min, InputPixelType max )
{
  if( this->m_Min != min || this->m_Max != max )
    {
    itkDebugMacro( "setting Min to " << min << "and Max to " << max );
    this->m_Min = min;
    this->m_Max = max;
    this->Modified();
    }
}
template< class TInputImage, class TOutputImage>
void
TextureFeaturesImageFilter< TInputImage, TOutputImage>
::GenerateOutputInformation()
{
  // this methods is overloaded so that if the output image is a
  // VectorImage then the correct number of components are set.
  Superclass::GenerateOutputInformation();
  OutputImageType* output = this->GetOutput();

  if ( !output )
    {
    return;
    }
  if ( output->GetNumberOfComponentsPerPixel() != this->GetNumberOfOutputComponents() )
    {
    output->SetNumberOfComponentsPerPixel( this->GetNumberOfOutputComponents() );
    }
}


template<class TInputImage, class TOutputImage>
void
TextureFeaturesImageFilter<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
  this->m_CooccurenceOffsetVector.clear();

  OffsetType o1;
  OffsetType o2;

  // calculate all offsets pairs
  for( unsigned int d = 0; d < ImageDimension; ++d )
    {
    o1[d] = static_cast<OffsetValueType>( -this->m_NeighborhoodRadius[d] );
    }

  // iterate over all point in the window
  do
    {
    typename OffsetVectorType::const_iterator it;
    for( it = this->m_Offsets.begin(); it != this->m_Offsets.end(); ++it )
      {
      o2 = o1 + *it;

      bool ok = true;
      for( unsigned int d = 0; d < ImageDimension; ++d )
        {
        if( vcl_abs( o2[d] ) > this->m_NeighborhoodRadius[d] )
          {
          ok = false;
          break;
          }
        }

      if( ok )
        {
        this->m_CooccurenceOffsetVector.push_back( std::make_pair( o1, o2 ) );
        }
      }

    // increment the offset 1
    ++o1[0];
    for( unsigned int d = 1; d < ImageDimension; ++d )
      {
      if( o1[d-1] > static_cast<OffsetValueType>( this->m_NeighborhoodRadius[d-1] ) )
        {
        o1[d-1] = -this->m_NeighborhoodRadius[d-1];
        ++o1[d];
        }
      }
    } while ( o1[ImageDimension-1] <= static_cast<OffsetValueType>( this->m_NeighborhoodRadius[ImageDimension-1] ) );
}


template<class TInputImage, class TOutputImage>
void
TextureFeaturesImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateData( const RegionType & region, ThreadIdType threadId )
{
  const InputImageType *inputImage = this->GetInput();
  OutputImageType      *outputImage = this->GetOutput();
  const MaskImageType  *maskImage = this->GetMaskImage();

  // constant for a coocurrence matrix.
  const unsigned int measurementVectorSize = 2;

  //
  // Create Histogram
  //
  typename HistogramType::Pointer histogram = HistogramType::New();
  histogram->SetMeasurementVectorSize( measurementVectorSize );

  //initialize parameters
  MeasurementVectorType lowerBound( measurementVectorSize );
  MeasurementVectorType upperBound( measurementVectorSize );
  lowerBound.Fill( this->GetMin() );
  upperBound.Fill( this->GetMax() +1 );

  // First, create an appropriate histogram with the right number of bins
  // and mins and maxs correct for the image type.
  typename HistogramType::SizeType size( histogram->GetMeasurementVectorSize() );
  size.Fill( this->m_NumberOfBinsPerAxis );
  histogram->Initialize( size, lowerBound, upperBound );

  typedef Statistics::HistogramToTextureFeaturesFilter<HistogramType> FeatureFilterType;
  typename FeatureFilterType::Pointer featureFilter = FeatureFilterType::New();

  ProgressReporter progress( this, threadId, region.GetNumberOfPixels() );

  typedef typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType> FaceCalculatorType;
  FaceCalculatorType faceCalculator;

  typename FaceCalculatorType::FaceListType faceList = faceCalculator( inputImage, region, this->m_NeighborhoodRadius );

  typename FaceCalculatorType::FaceListType::iterator fit;

  for( fit = faceList.begin(); fit != faceList.end(); ++fit )
    {
    ConstNeighborhoodIteratorType It( this->m_NeighborhoodRadius, inputImage, *fit );
    NeighborhoodIterator<OutputImageType> ItO( this->m_NeighborhoodRadius, outputImage, *fit );

    OutputPixelType out;
    NumericTraits<OutputPixelType>::SetLength( out, this->GetNumberOfOutputComponents() );

    typename FeatureFilterType::Pointer featureFilter = FeatureFilterType::New();

    typename InputImageType::IndexType centerIndex = It.GetIndex();

    for( It.GoToBegin(), ItO.GoToBegin(); !It.IsAtEnd(); ++It, ++ItO )
      {
      if( maskImage && ( maskImage->GetPixel( centerIndex ) == this->m_InsidePixelValue ) )
        {
        continue;
        }
      MeasurementVectorType cooccur( histogram->GetMeasurementVectorSize() );

      typename std::vector<std::pair<OffsetType, OffsetType > >::const_iterator it;
      for( it = this->m_CooccurenceOffsetVector.begin(); it != this->m_CooccurenceOffsetVector.end(); ++it )
        {
        const InputPixelType p1 = It.GetPixel( it->first );
        const InputPixelType p2 = It.GetPixel( it->second );

        if( maskImage &&
          ( maskImage->GetPixel( centerIndex + it->first ) == this->m_InsidePixelValue ||
            maskImage->GetPixel( centerIndex + it->second ) == this->m_InsidePixelValue ) )
          {
          continue;
          }

        if( p1 >= this->m_Min && p2 >= this->m_Min && p1 <= this->m_Max && p2 <= this->m_Max )
          {
          cooccur[0] = p1;
          cooccur[1] = p2;
          histogram->IncreaseFrequencyOfMeasurement( cooccur, 1.0 );

          cooccur[1] = p1;
          cooccur[0] = p2;
          histogram->IncreaseFrequencyOfMeasurement( cooccur, 1.0 );
          }

        ++it;
        }

      featureFilter->SetInput( histogram );
      featureFilter->Modified();
      featureFilter->Update();

      out[0] = featureFilter->GetEnergy();
      out[1] = featureFilter->GetEntropy();
      out[2] = featureFilter->GetCorrelation();
      out[3] = featureFilter->GetInverseDifferenceMoment();
      out[4] = featureFilter->GetInertia();
      out[5] = featureFilter->GetClusterShade();
      out[6] = featureFilter->GetClusterProminence();
      out[7] = featureFilter->GetHaralickCorrelation();

      ItO.SetCenterPixel( out );

      histogram->SetToZero();
      progress.CompletedPixel();
      }
    }
}

template<class TInputImage, class TOutputImage>
void
TextureFeaturesImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

//   os << indent << "Offsets: " << this->m_Offsets << std::endl;
  os << indent << "Min: " << this->GetMin() << std::endl;
  os << indent << "Max: " << this->GetMax() << std::endl;
  os << indent << "NumberOfBinsPerAxis: " << this->GetNumberOfBinsPerAxis() << std::endl;
  os << indent << "Normalize: " << this->GetNormalize() << std::endl;
  }

} // end namespace Statistics
} // end namespace itk

#endif

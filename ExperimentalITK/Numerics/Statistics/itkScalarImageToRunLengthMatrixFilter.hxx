/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkScalarImageToRunLengthMatrixFilter.hxx,v $
  Language:  C++
  Date:      $Date: 2009/05/02 03:00:26 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkScalarImageToRunLengthMatrixFilter_hxx
#define _itkScalarImageToRunLengthMatrixFilter_hxx

#include "itkScalarImageToRunLengthMatrixFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhood.h"
#include "vnl/vnl_math.h"


namespace itk {
namespace Statistics {

template<class TImageType, class THistogramFrequencyContainer>
ScalarImageToRunLengthMatrixFilter<TImageType, THistogramFrequencyContainer>
::ScalarImageToRunLengthMatrixFilter()
{
  this->SetNumberOfRequiredInputs( 1 );
  this->SetNumberOfRequiredOutputs( 1 );

  this->ProcessObject::SetNthOutput( 0, this->MakeOutput( 0 ) );

  const unsigned int measurementVectorSize = 2;

  HistogramType * output = const_cast< HistogramType *>( this->GetOutput() );

  output->SetMeasurementVectorSize( measurementVectorSize );

  //initialize parameters
  this->m_LowerBound.SetSize( measurementVectorSize );
  this->m_UpperBound.SetSize( measurementVectorSize );

  this->m_LowerBound.Fill(NumericTraits<PixelType>::NonpositiveMin());
  this->m_UpperBound.Fill(NumericTraits<PixelType>::max() + 1);

  this->m_Min = NumericTraits<PixelType>::NonpositiveMin();
  this->m_Max = NumericTraits<PixelType>::max();

  this->m_MinDistance = NumericTraits<PixelType>::NonpositiveMin();
  this->m_MaxDistance = NumericTraits<PixelType>::max();

  //mask inside pixel value
  this->m_InsidePixelValue = NumericTraits<PixelType>::One;

  this->m_NumberOfBinsPerAxis = DefaultBinsPerAxis;

  // Get a set of default offset values.
  typedef Neighborhood<PixelType, ImageDimension> NeighborhoodType;
  NeighborhoodType neighborhood;
  typename NeighborhoodType::RadiusType radius;
  radius.Fill( 1 );
  neighborhood.SetRadius( radius );

  this->m_Offsets = OffsetVector::New();
  unsigned int numberOfOffsets = static_cast<unsigned int>(
    0.5 * neighborhood.Size() );
  for( unsigned int n = 0; n < numberOfOffsets; n++ )
    {
    this->AddOffset( neighborhood.GetOffset( n ) );
    }
  }

template<class TImageType, class THistogramFrequencyContainer>
void
ScalarImageToRunLengthMatrixFilter<TImageType, THistogramFrequencyContainer>
::SetOffset( const OffsetType offset )
{
  OffsetVectorPointer offsetVector = OffsetVector::New();
  offsetVector->push_back( offset );
  this->SetOffsets( offsetVector );
}

template<class TImageType, class THistogramFrequencyContainer>
void
ScalarImageToRunLengthMatrixFilter<TImageType, THistogramFrequencyContainer>
::AddOffset( const OffsetType offset )
{
  if( this->m_Offsets )
    {
    this->m_Offsets->push_back( offset );
    }
  else
    {
    this->SetOffset( offset );
    }
}

template<class TImageType, class THistogramFrequencyContainer>
void
ScalarImageToRunLengthMatrixFilter<TImageType,
THistogramFrequencyContainer>
::SetInput(const ImageType* image)
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( 0,
    const_cast< ImageType*>( image ) );
}
template<class TImageType, class THistogramFrequencyContainer>
void
ScalarImageToRunLengthMatrixFilter<TImageType,
THistogramFrequencyContainer>
::SetMaskImage(const ImageType* image)
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( 1,
    const_cast< ImageType*>( image ) );
}


template<class TImageType, class THistogramFrequencyContainer>
const TImageType*
ScalarImageToRunLengthMatrixFilter<TImageType, THistogramFrequencyContainer>
::GetInput() const
{
  if (this->GetNumberOfInputs() < 1)
    {
    return 0;
    }
  return static_cast<const ImageType *>( this->ProcessObject::GetInput( 0 ) );
}

template<class TImageType, class THistogramFrequencyContainer>
const TImageType*
ScalarImageToRunLengthMatrixFilter<TImageType, THistogramFrequencyContainer>
::GetMaskImage() const
{
  if ( this->GetNumberOfInputs() < 2 )
    {
    return 0;
    }
  return static_cast<const ImageType *>( this->ProcessObject::GetInput( 1 ) );
}

template<class TImageType, class THistogramFrequencyContainer>
const typename ScalarImageToRunLengthMatrixFilter<TImageType,
THistogramFrequencyContainer>::HistogramType *
ScalarImageToRunLengthMatrixFilter<TImageType, THistogramFrequencyContainer>
::GetOutput() const
{
  const HistogramType *output =
    static_cast<const HistogramType *>( this->ProcessObject::GetOutput( 0 ) );
  return output;
}

template<class TImageType, class THistogramFrequencyContainer>
typename ScalarImageToRunLengthMatrixFilter<TImageType,
THistogramFrequencyContainer>::DataObjectPointer
ScalarImageToRunLengthMatrixFilter<TImageType,
THistogramFrequencyContainer>
::MakeOutput( unsigned int itkNotUsed( idx ) )
{
  typename HistogramType::Pointer output = HistogramType::New();
  return static_cast< DataObject *>( output );
}

template<class TImageType, class THistogramFrequencyContainer>
void
ScalarImageToRunLengthMatrixFilter<TImageType,
THistogramFrequencyContainer>::
GenerateData( void )
{
  HistogramType * output =
   static_cast< HistogramType * >( this->ProcessObject::GetOutput( 0 ) );

  const ImageType *input = this->GetInput();
  // At this point input must be non-NULL because the ProcessObject
  // checks the number of required input to be non-NULL pointers before
  // calling this GenerateData() method.

  // First, create an appropriate histogram with the right number of bins
  // and mins and maxes correct for the image type.
  typename HistogramType::SizeType size( output->GetMeasurementVectorSize() );
  size.Fill( this->m_NumberOfBinsPerAxis );
  output->Initialize( size, this->m_LowerBound, this->m_UpperBound );

  // Next, find the minimum radius that encloses all the offsets.
  unsigned int minRadius = 0;
  typename OffsetVector::ConstIterator offsets;
  for( offsets = this->m_Offsets->Begin(); offsets != this->m_Offsets->End();
    offsets++ )
    {
    for( unsigned int i = 0; i < offsets.Value().GetOffsetDimension(); i++ )
      {
      unsigned int distance = vnl_math_abs( offsets.Value()[i] );
      if( distance > minRadius )
        {
        minRadius = distance;
        }
      }
    }

  RadiusType radius;
  radius.Fill( minRadius );

  const ImageType *maskImage = NULL;

  // Check if a mask image has been provided
  //
  if ( this->GetNumberOfInputs() > 1 )
    {
    maskImage = this->GetMaskImage();
    }

  // Now fill in the histogram
  if ( maskImage != NULL )
    {
    this->FillHistogramWithMask( radius, input->GetRequestedRegion(), maskImage );
    }
  else
    {
    this->FillHistogram( radius, input->GetRequestedRegion() );
    }
}

template<class TImageType, class THistogramFrequencyContainer>
void
ScalarImageToRunLengthMatrixFilter<TImageType,
THistogramFrequencyContainer>::
FillHistogramWithMask( RadiusType radius, RegionType region, const ImageType *maskImage )
{
  // Iterate over all of those pixels and offsets, adding each
  // co-occurrence pair to the histogram

  const ImageType *input = this->GetInput();

  HistogramType * output =
   static_cast< HistogramType * >( this->ProcessObject::GetOutput( 0 ) );

  typedef ConstNeighborhoodIterator<ImageType> NeighborhoodIteratorType;
  NeighborhoodIteratorType neighborIt, maskNeighborIt;
  neighborIt = NeighborhoodIteratorType( radius, input, region );
  maskNeighborIt = NeighborhoodIteratorType( radius, maskImage, region );

  typename OffsetVector::ConstIterator offsets;
  for( offsets = this->GetOffsets()->Begin();
    offsets != this->GetOffsets()->End(); offsets++ )
    {
    typedef Image<bool, ImageDimension> BoolImageType;
    typename BoolImageType::Pointer alreadyVisitedImage = BoolImageType::New();
    alreadyVisitedImage->SetRegions( input->GetRequestedRegion() );
    alreadyVisitedImage->SetOrigin( input->GetOrigin() );
    alreadyVisitedImage->SetSpacing( input->GetSpacing() );
    alreadyVisitedImage->Allocate();
    alreadyVisitedImage->FillBuffer( false );

    neighborIt.GoToBegin();
    OffsetType offset = offsets.Value();

    for(neighborIt.GoToBegin(), maskNeighborIt.GoToBegin();
      !neighborIt.IsAtEnd(); ++neighborIt, ++maskNeighborIt )
      {
      if ( maskNeighborIt.GetCenterPixel() != this->m_InsidePixelValue )
        {
        continue; // Go to the next loop if we're not in the mask
        }
      const PixelType centerPixelIntensity = neighborIt.GetCenterPixel();
      IndexType centerIndex = neighborIt.GetIndex();
      if ( centerPixelIntensity < this->m_Min ||
           centerPixelIntensity > this->m_Max ||
           alreadyVisitedImage->GetPixel( centerIndex ) )
        {
        continue; // don't put a pixel in the histogram if the value
                  // is out-of-bounds.
        }

      MeasurementType centerBinMin = output->
        GetBinMinFromValue( 0, centerPixelIntensity );
      MeasurementType centerBinMax = output->
        GetBinMaxFromValue( 0, centerPixelIntensity );

      IndexType index = centerIndex;
      PixelType pixelIntensity = input->GetPixel( index );
      while ( pixelIntensity >= centerBinMin &&
              pixelIntensity <= centerBinMax &&
              !alreadyVisitedImage->GetPixel( index ) )
        {
        alreadyVisitedImage->SetPixel( index, true );
        index += offset;
        if ( input->GetRequestedRegion().IsInside( index ) )
          {
          pixelIntensity = input->GetPixel( index );
          }
        else
          {
          break;
          }
        }

      PointType centerPoint;
      input->TransformIndexToPhysicalPoint( centerIndex, centerPoint );
      PointType point;
      input->TransformIndexToPhysicalPoint( index, point );

      MeasurementVectorType run( output->GetMeasurementVectorSize() );
      run[0] = centerPixelIntensity;
      run[1] = centerPoint.EuclideanDistanceTo( point );

      if ( run[1] >= this->m_MinDistance &&
           run[1] <= this->m_MaxDistance )
        {
        output->IncreaseFrequencyOfMeasurement( run, 1 );
        }
      }
    }
}

template<class TImageType, class THistogramFrequencyContainer>
void
ScalarImageToRunLengthMatrixFilter<TImageType,
THistogramFrequencyContainer>::
FillHistogram( RadiusType radius, RegionType region )
{
   // Iterate over all of those pixels and offsets, adding each
  // co-occurrence pair to the histogram

  const ImageType *input = this->GetInput();

  HistogramType * output =
   static_cast< HistogramType * >( this->ProcessObject::GetOutput( 0 ) );

  typedef ConstNeighborhoodIterator<ImageType> NeighborhoodIteratorType;
  NeighborhoodIteratorType neighborIt( radius,
    input, input->GetRequestedRegion() );

  typename OffsetVector::ConstIterator offsets;
  for( offsets = this->GetOffsets()->Begin();
    offsets != this->GetOffsets()->End(); offsets++ )
    {
    typedef Image<bool, ImageDimension> BoolImageType;
    typename BoolImageType::Pointer alreadyVisitedImage = BoolImageType::New();
    alreadyVisitedImage->SetRegions( input->GetRequestedRegion() );
    alreadyVisitedImage->SetOrigin( input->GetOrigin() );
    alreadyVisitedImage->SetSpacing( input->GetSpacing() );
    alreadyVisitedImage->Allocate();
    alreadyVisitedImage->FillBuffer( false );

    neighborIt.GoToBegin();
    OffsetType offset = offsets.Value();

    for ( neighborIt.GoToBegin(); !neighborIt.IsAtEnd(); ++neighborIt )
      {
      const PixelType centerPixelIntensity = neighborIt.GetCenterPixel();
      IndexType centerIndex = neighborIt.GetIndex();
      if ( centerPixelIntensity < this->m_Min ||
           centerPixelIntensity > this->m_Max ||
           alreadyVisitedImage->GetPixel( centerIndex ) )
        {
        continue; // don't put a pixel in the histogram if the value
                  // is out-of-bounds.
        }

      MeasurementType centerBinMin = output->
        GetBinMinFromValue( 0, centerPixelIntensity );
      MeasurementType centerBinMax = output->
        GetBinMaxFromValue( 0, centerPixelIntensity );

      IndexType index = centerIndex;
      PixelType pixelIntensity = input->GetPixel( index );
      while ( pixelIntensity >= centerBinMin &&
              pixelIntensity <= centerBinMax &&
              !alreadyVisitedImage->GetPixel( index ) )
        {
        alreadyVisitedImage->SetPixel( index, true );
        index += offset;
        if ( input->GetRequestedRegion().IsInside( index ) )
          {
          pixelIntensity = input->GetPixel( index );
          }
        else
          {
          break;
          }
        }

      PointType centerPoint;
      input->TransformIndexToPhysicalPoint( centerIndex, centerPoint );
      PointType point;
      input->TransformIndexToPhysicalPoint( index, point );

      MeasurementVectorType run( output->GetMeasurementVectorSize() );
      run[0] = centerPixelIntensity;
      run[1] = centerPoint.EuclideanDistanceTo( point );

      if ( run[1] >= this->m_MinDistance &&
           run[1] <= this->m_MaxDistance )
        {
        output->IncreaseFrequencyOfMeasurement( run, 1 );
        }
      }
    }
}

template<class TImageType, class THistogramFrequencyContainer>
void
ScalarImageToRunLengthMatrixFilter<TImageType,
THistogramFrequencyContainer>::
SetPixelValueMinMax( PixelType min, PixelType max )
{
  itkDebugMacro("setting Min to " << min << "and Max to " << max);
  this->m_Min = min;
  this->m_Max = max;
  this->m_LowerBound[0] = min;
  this->m_UpperBound[0] = max;
  this->Modified();
}

template<class TImageType, class THistogramFrequencyContainer>
void
ScalarImageToRunLengthMatrixFilter<TImageType,
THistogramFrequencyContainer>::
SetDistanceValueMinMax( RealType min, RealType max )
  {
  itkDebugMacro( "setting MinDistance to " << min <<
                 "and MaxDistance to " << max );
  this->m_MinDistance = min;
  this->m_MaxDistance = max;
  this->m_LowerBound[1] = min;
  this->m_UpperBound[1] = max;
  this->Modified();
  }

template<class TImageType, class THistogramFrequencyContainer>
void
ScalarImageToRunLengthMatrixFilter<TImageType,
THistogramFrequencyContainer>::
PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "Offsets: " << this->GetOffsets() << std::endl;
  os << indent << "Min: " << this->GetMin() << std::endl;
  os << indent << "Max: " << this->GetMax() << std::endl;
  os << indent << "MinDistance: " << this->GetMinDistance() << std::endl;
  os << indent << "MaxDistance: " << this->GetMaxDistance() << std::endl;
  os << indent << "NumberOfBinsPerAxis: "
     << this->GetNumberOfBinsPerAxis() << std::endl;
  os << indent << "InsidePixelValue: "
     << this->GetInsidePixelValue() << std::endl;
}

} // end of namespace Statistics
} // end of namespace itk


#endif

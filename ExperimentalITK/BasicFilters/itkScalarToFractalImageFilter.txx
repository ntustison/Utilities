/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkScalarToFractalImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2009/02/12 22:59:54 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkScalarToFractalImageFilter_txx
#define __itkScalarToFractalImageFilter_txx

#include "itkScalarToFractalImageFilter.h"

#include "itkNeighborhoodAlgorithm.h"
#include "itkProgressReporter.h"

#include <vector>

namespace itk {

template <class TInputImage, class TOutputImage>
ScalarToFractalImageFilter<TInputImage, TOutputImage>
::ScalarToFractalImageFilter()
{
  this->m_NeighborhoodRadius.Fill( 2 );

  this->m_MaskImage = NULL;
}

template<class TInputImage, class TOutputImage>
void
ScalarToFractalImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  this->AllocateOutputs();

  ProgressReporter progress( this, 0, 
    this->GetInput()->GetRequestedRegion().GetNumberOfPixels(), 100 );

  typedef typename NeighborhoodAlgorithm
    ::ImageBoundaryFacesCalculator<InputImageType> FaceCalculatorType;
  FaceCalculatorType faceCalculator;

  typename FaceCalculatorType::FaceListType faceList
    = faceCalculator( this->GetInput(),
    this->GetInput()->GetRequestedRegion(), this->m_NeighborhoodRadius );
  typename FaceCalculatorType::FaceListType::iterator fit;

  RealType minSpacing = this->GetInput()->GetSpacing()[0];
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    if( this->GetInput()->GetSpacing()[d] < minSpacing )
      {
      minSpacing = this->GetInput()->GetSpacing()[d];
      }
    }

  std::vector<RealType> distances;
  std::vector<RealType> distancesFrequency;
  std::vector<RealType> averageAbsoluteIntensityDifference;

  for( fit = faceList.begin(); fit != faceList.end(); ++fit )
    {
    ConstNeighborhoodIteratorType It(
      this->m_NeighborhoodRadius, this->GetInput(), *fit );
    NeighborhoodIterator<OutputImageType> ItO(
      this->m_NeighborhoodRadius, this->GetOutput(), *fit );

    for( It.GoToBegin(), ItO.GoToBegin(); !It.IsAtEnd(); ++It, ++ItO )
      {
      if( this->m_MaskImage &&
        !this->m_MaskImage->GetPixel( It.GetIndex() ) )
        {
        ItO.SetCenterPixel(
          NumericTraits<typename OutputImageType::PixelType>::Zero );
        progress.CompletedPixel();
        continue;
        }
      
      distances.clear();
      distancesFrequency.clear();
      averageAbsoluteIntensityDifference.clear();

      for( unsigned int i = 0; i < It.GetNeighborhood().Size(); i++ )
        {
        bool IsInBounds1;
        typename InputImageType::PixelType pixel1
          = It.GetPixel( i, IsInBounds1 );

        if( IsInBounds1 && ( !this->m_MaskImage || ( this->m_MaskImage &&
          this->m_MaskImage->GetPixel( It.GetIndex( i ) ) ) ) )
          {
          typename InputImageType::PointType point1;
          this->GetInput()->TransformIndexToPhysicalPoint(
            It.GetIndex( i ), point1 );

          for( unsigned int j = 0; j < It.GetNeighborhood().Size(); j++ )
            {
            if( i != j )
              {
              bool IsInBounds2;
              typename InputImageType::PixelType pixel2
                = It.GetPixel( j, IsInBounds2 );

              if( IsInBounds2 && ( !this->m_MaskImage || ( this->m_MaskImage &&
                this->m_MaskImage->GetPixel( It.GetIndex( j ) ) ) ) )
                {
                typename InputImageType::PointType point2;
                this->GetInput()->TransformIndexToPhysicalPoint(
                  It.GetIndex( j ), point2 );

                RealType distance
                  = point1.SquaredEuclideanDistanceTo( point2 );

                bool distanceFound = false;
                for( unsigned int k = 0; k < distances.size(); k++ )
                  {
                  if( vnl_math_abs( distances[k] - distance )
                    < 0.5 * minSpacing )
                    {
                    distancesFrequency[k]++;
                    averageAbsoluteIntensityDifference[k]
                      += vnl_math_abs( pixel1 - pixel2 );
                    distanceFound = true;
                    break;
                    }
                  }
                if( !distanceFound )
                  {
                  distances.push_back( distance );
                  distancesFrequency.push_back( 1 );
                  averageAbsoluteIntensityDifference.push_back(
                    vnl_math_abs( pixel1 - pixel2 ) );
                  }
                }
              }
            }
          }
        }
      RealType sumY = 0.0;
      RealType sumX = 0.0;
      RealType sumXY = 0.0;
      RealType sumXX = 0.0;

      for( unsigned int k = 0; k < distances.size(); k++ )
        {
        if( distancesFrequency[k] == 0 )
          {
          continue;
          }

        averageAbsoluteIntensityDifference[k]
          /= static_cast<RealType>( distancesFrequency[k] );
        averageAbsoluteIntensityDifference[k]
          = vcl_log( averageAbsoluteIntensityDifference[k] );

        RealType distance = vcl_log( vcl_sqrt( distances[k] ) );

        sumY += averageAbsoluteIntensityDifference[k];
        sumX += distance;
        sumXX += ( distance * distance );
        sumXY += ( averageAbsoluteIntensityDifference[k] * distance );
        }
      RealType N = static_cast<RealType>( distances.size() );

      RealType slope = ( N * sumXY - sumX * sumY )
        / ( N * sumXX - sumX * sumX );

      ItO.SetCenterPixel(
        static_cast<typename OutputImageType::PixelType>( 3.0 - slope ) );
      progress.CompletedPixel();
      }
    }
}

template<class TInputImage, class TOutputImage>
void
ScalarToFractalImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Neighborhood radius: "
    << this->m_NeighborhoodRadius << std::endl;
}

}// end namespace itk
#endif

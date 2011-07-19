/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRobustDemonsRegistrationFunction.hxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:13:43 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkRobustDemonsRegistrationFunction_hxx_
#define _itkRobustDemonsRegistrationFunction_hxx_

#include "itkRobustDemonsRegistrationFunction.h"
#include "itkExceptionObject.h"
#include "vnl/vnl_math.h"

namespace itk {

/*
 * Default constructor
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
RobustDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::RobustDemonsRegistrationFunction()
{

  RadiusType r;
  unsigned int j;
  for( j = 0; j < ImageDimension; j++ )
    {
    r[j] = 0;
    }
  this->SetRadius(r);

  this->m_TimeStep = 1.0;
  this->m_DenominatorThreshold = 1e-9;
  this->m_IntensityDifferenceThreshold = 0.001;
  this->m_MovingImage = NULL;
  this->m_FixedImage = NULL;
  this->m_FixedImageSpacing.Fill( 1.0 );
  this->m_FixedImageOrigin.Fill( 0.0 );
  this->m_Normalizer = 1.0;
  this->m_FixedImageGradientCalculator = GradientCalculatorType::New();
  
  this->m_Robustness=1.0;

  typename DefaultInterpolatorType::Pointer interp =
    DefaultInterpolatorType::New();

  this->m_MovingImageInterpolator = static_cast<InterpolatorType*>(
    interp.GetPointer() );

  this->m_Metric = NumericTraits<double>::max();
  this->m_SumOfSquaredDifference = 0.0;
  this->m_NumberOfPixelsProcessed = 0L;
  this->m_RMSChange = NumericTraits<double>::max();
  this->m_SumOfSquaredChange = 0.0;

}


/*
 * Standard "PrintSelf" method.
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
RobustDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "MovingImageIterpolator: ";
  os << this->m_MovingImageInterpolator.GetPointer() << std::endl;
  os << indent << "FixedImageGradientCalculator: ";
  os << this->m_FixedImageGradientCalculator.GetPointer() << std::endl;
  os << indent << "DenominatorThreshold: ";
  os << this->m_DenominatorThreshold << std::endl;
  os << indent << "IntensityDifferenceThreshold: ";
  os << this->m_IntensityDifferenceThreshold << std::endl;

  os << indent << "Metric: ";
  os << this->m_Metric << std::endl;
  os << indent << "SumOfSquaredDifference: ";
  os << this->m_SumOfSquaredDifference << std::endl;
  os << indent << "NumberOfPixelsProcessed: ";
  os << this->m_NumberOfPixelsProcessed << std::endl;
  os << indent << "RMSChange: ";
  os << this->m_RMSChange << std::endl;
  os << indent << "SumOfSquaredChange: ";
  os << this->m_SumOfSquaredChange << std::endl;

}


/*
 * Set the function state values before each iteration
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
RobustDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::InitializeIteration()
{
  if( !this->m_MovingImage || !this->m_FixedImage || !this->m_MovingImageInterpolator )
    {
    itkExceptionMacro( << "MovingImage, FixedImage and/or Interpolator not set" );
    }

  // cache fixed image information
  this->m_FixedImageSpacing    = this->m_FixedImage->GetSpacing();
  this->m_FixedImageOrigin     = this->m_FixedImage->GetOrigin();

  // compute the normalizer
  this->m_Normalizer      = 0.0;
  for( unsigned int k = 0; k < ImageDimension; k++ )
    {
    this->m_Normalizer += this->m_FixedImageSpacing[k] * this->m_FixedImageSpacing[k];
    }
  this->m_Normalizer /= static_cast<double>( ImageDimension );


  // setup gradient calculator
  this->m_FixedImageGradientCalculator->SetInputImage( this->m_FixedImage );

  // setup moving image interpolator
  this->m_MovingImageInterpolator->SetInputImage( this->m_MovingImage );

  // initialize metric computation variables
  this->m_SumOfSquaredDifference  = 0.0;
  this->m_NumberOfPixelsProcessed = 0L;
  this->m_SumOfSquaredChange      = 0.0;

}


/*
 * Compute update at a specify neighbourhood
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
typename RobustDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::PixelType
RobustDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::ComputeUpdate(const NeighborhoodType &it, void * gd,
                const FloatOffsetType& itkNotUsed(offset))
{

  PixelType update;
  unsigned int j;

  IndexType index = it.GetIndex();

  // Get fixed image related information
  double fixedValue;
  CovariantVectorType fixedGradient;
  double fixedGradientSquaredMagnitude = 0;

  // Note: no need to check the index is within
  // fixed image buffer. This is done by the external filter.
  fixedValue = (double) this->m_FixedImage->GetPixel( index );
  fixedGradient = this->m_FixedImageGradientCalculator->EvaluateAtIndex( index );
  for( unsigned int j = 0; j < ImageDimension; j++ )
    {
    fixedGradientSquaredMagnitude += vnl_math_sqr( fixedGradient[j] );
    } 

  // Get moving image related information
  double movingValue;
  PointType mappedPoint;

  for( j = 0; j < ImageDimension; j++ )
    {
    mappedPoint[j] = double( index[j] ) * this->m_FixedImageSpacing[j] + 
      this->m_FixedImageOrigin[j];
    mappedPoint[j] += it.GetCenterPixel()[j];
    }
  if( this->m_MovingImageInterpolator->IsInsideBuffer( mappedPoint ) )
    {
    movingValue = this->m_MovingImageInterpolator->Evaluate( mappedPoint );
    }
  else
    {
    for( j = 0; j < ImageDimension; j++ )
      {
      update[j] = 0.0;
      }
    return update;
    }

  /**
   * Compute Update.
   * In the original equation the denominator is defined as (g-f)^2 + grad_mag^2.
   * However there is a mismatch in units between the two terms. 
   * The units for the second term is intensity^2/mm^2 while the
   * units for the first term is intensity^2. This mismatch is particularly
   * problematic when the fixed image does not have unit spacing.
   * In this implemenation, we normalize the first term by a factor K,
   * such that denominator = (g-f)^2/K + grad_mag^2
   * where K = mean square spacing to compensate for the mismatch in units.
   */
  double speedValue = fixedValue - movingValue;
  
  if (  speedValue > this->m_Robustness ) speedValue=this->m_Robustness;
  if (  (speedValue*(-1.)) > this->m_Robustness ) speedValue=(-1.0)*this->m_Robustness;
  
  // update the metric
  GlobalDataStruct *globalData = (GlobalDataStruct *)gd;
  if ( globalData )
    {
    globalData->m_SumOfSquaredDifference += vnl_math_sqr( speedValue );
    this->m_Energy = globalData->m_SumOfSquaredDifference;
    globalData->m_NumberOfPixelsProcessed += 1;
    }

  double denominator = vnl_math_sqr( speedValue ) / this->m_Normalizer + 
    fixedGradientSquaredMagnitude;

  if ( vnl_math_abs(speedValue) < this->m_IntensityDifferenceThreshold || 
       denominator < this->m_DenominatorThreshold  )
    {
    for( j = 0; j < ImageDimension; j++ )
      {
      update[j] = 0.0;
      }
    return update;
    }

  for( j = 0; j < ImageDimension; j++ )
    {
    update[j] = speedValue * fixedGradient[j] / denominator;
    if ( globalData )
      {
      globalData->m_SumOfSquaredChange += vnl_math_sqr( update[j] );
      }
    }

  return update;

}

/*
 * Update the metric and release the per-thread-global data.
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
RobustDemonsRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::ReleaseGlobalDataPointer( void *gd ) const
{
  GlobalDataStruct * globalData = (GlobalDataStruct *) gd;

  this->m_MetricCalculationLock.Lock();
  this->m_SumOfSquaredDifference  += globalData->m_SumOfSquaredDifference;
  this->m_NumberOfPixelsProcessed += globalData->m_NumberOfPixelsProcessed;
  this->m_SumOfSquaredChange += globalData->m_SumOfSquaredChange;
  if ( this->m_NumberOfPixelsProcessed )
    {
    this->m_Metric = this->m_SumOfSquaredDifference / 
               static_cast<double>( this->m_NumberOfPixelsProcessed ); 
    this->m_RMSChange = vcl_sqrt( this->m_SumOfSquaredChange / 
               static_cast<double>( this->m_NumberOfPixelsProcessed ) ); 
    }
  this->m_MetricCalculationLock.Unlock();

  delete globalData;
}



} // end namespace itk

#endif

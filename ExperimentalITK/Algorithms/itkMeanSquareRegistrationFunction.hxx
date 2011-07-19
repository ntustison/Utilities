/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMeanSquareRegistrationFunction.hxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:13:43 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkMeanSquareRegistrationFunction_hxx_
#define _itkMeanSquareRegistrationFunction_hxx_

#include "itkMeanSquareRegistrationFunction.h"
#include "itkExceptionObject.h"
#include "vnl/vnl_math.h"

namespace itk {

/*
 * Default constructor
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
MeanSquareRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::MeanSquareRegistrationFunction()
{

  RadiusType r;
  unsigned int j;
  for( j = 0; j < ImageDimension; j++ )
    {
    r[j] = 0;
    }
  this->SetRadius(r);

  m_Normalizer = 1.0;
   Superclass::m_Energy = 0.0;
  m_TimeStep = 1.0;
  m_DenominatorThreshold = 1e-9;
  m_IntensityDifferenceThreshold = 0.001;
  Superclass::m_MovingImage = NULL;
  Superclass::m_FixedImage = NULL;
  m_FixedImageSpacing.Fill( 1.0 );
  m_FixedImageOrigin.Fill( 0.0 );
  m_FixedImageGradientCalculator = GradientCalculatorType::New();
  m_MovingImageGradientCalculator = GradientCalculatorType::New();


  typename DefaultInterpolatorType::Pointer interp =
    DefaultInterpolatorType::New();

  m_MovingImageInterpolator = static_cast<InterpolatorType*>(
    interp.GetPointer() );
    
   Superclass::m_GradientStep=1.e3;


  bool m_Symmetric=false;
  bool m_Robust=true;
  bool m_MovingGradient=false;
  
}


/*
 * Standard "PrintSelf" method.
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
MeanSquareRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  std::cout << " mean squar " << std::endl;
/*
  os << indent << "MovingImageIterpolator: ";
  os << m_MovingImageInterpolator.GetPointer() << std::endl;
  os << indent << "FixedImageGradientCalculator: ";
  os << m_FixedImageGradientCalculator.GetPointer() << std::endl;
  os << indent << "DenominatorThreshold: ";
  os << m_DenominatorThreshold << std::endl;
  os << indent << "IntensityDifferenceThreshold: ";
  os << m_IntensityDifferenceThreshold << std::endl;
*/
}


/*
 * Set the function state values before each iteration
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
MeanSquareRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::InitializeIteration()
{
  if( !Superclass::m_MovingImage || !Superclass::m_FixedImage || !m_MovingImageInterpolator )
    {
    itkExceptionMacro( << "MovingImage, FixedImage and/or Interpolator not set" );
    throw ExceptionObject(__FILE__,__LINE__);
    }

  // cache fixed image information
  m_FixedImageSpacing    = Superclass::m_FixedImage->GetSpacing();
  m_FixedImageOrigin     = Superclass::m_FixedImage->GetOrigin();

  // setup gradient calculator
  m_FixedImageGradientCalculator->SetInputImage( Superclass::m_FixedImage );
  m_MovingImageGradientCalculator->SetInputImage( Superclass::m_MovingImage );
  // compute the normalizer
  m_Normalizer      = 0.0;
  for( unsigned int k = 0; k < ImageDimension; k++ )
    {
    m_Normalizer += m_FixedImageSpacing[k] * m_FixedImageSpacing[k];
    }
  m_Normalizer /= static_cast<double>( ImageDimension );


  float scl=1.0;  
//  for (int i=0; i<ImageDimension; i++) scl*=(float)Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[i];
  // setup moving image interpolator
  m_MovingImageInterpolator->SetInputImage( Superclass::m_MovingImage );
//  std::cout << " MSQ metric " <<  Superclass::m_Energy/scl << std::endl;
   Superclass::m_Energy=0.0;
}


/*
 * Compute update at a non boundary neighbourhood
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
typename MeanSquareRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::PixelType
MeanSquareRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::ComputeUpdate(const NeighborhoodType &it, void * itkNotUsed(globalData),
                const FloatOffsetType& itkNotUsed(offset)) 
{

  PixelType update;
  unsigned int j;

  IndexType index = it.GetIndex();
  typedef typename TFixedImage::SizeType SizeType;
  SizeType size = Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize();
  for (j=0; j<ImageDimension;j++) 
    if (index[j] >= size[j]-1 || index[j] < 1 ) return update;

  // Get fixed image related information
  double fixedValue;
  CovariantVectorType fixedGradient;
  CovariantVectorType movingGradient;
  double fixedGradientSquaredMagnitude = 0;
  double movingGradientSquaredMagnitude = 0;

  // Note: no need to check the index is within
  // fixed image buffer. This is done by the external filter.
  fixedValue = (double) Superclass::m_FixedImage->GetPixel( index );
  fixedGradient = m_FixedImageGradientCalculator->EvaluateAtIndex( index );
  for( unsigned int j = 0; j < ImageDimension; j++ )
    {
    fixedGradientSquaredMagnitude += vnl_math_sqr( fixedGradient[j] ) ;
    } 

  // Get moving image related information
  double movingValue;
  PointType mappedPoint;
		typename TDeformationField::PixelType itvec;
		if ( Superclass::m_DeformationField )
				{
				itvec = Superclass::m_DeformationField->GetPixel(index);
				}
		else
				{
				itvec.Fill( 0 );
				} 
  
  IndexType movindex;
  for( j = 0; j < ImageDimension; j++ )
    {
     mappedPoint[j] = double( index[j] ) * m_FixedImageSpacing[j] + 
      m_FixedImageOrigin[j];
//     mappedPoint[j] += it.GetCenterPixel()[j];
      mappedPoint[j] += itvec[j];
      movindex[j] = (long)( double( index[j] ) +  itvec[j]/m_FixedImageSpacing[j] +0.5);
    }
  if( m_MovingImageInterpolator->IsInsideBuffer( mappedPoint ) )
    {
    movingValue = m_MovingImageInterpolator->Evaluate( mappedPoint );
//    movingGradient = m_MovingImageGradientCalculator->Evaluate( mappedPoint );
    movingGradient = m_MovingImageGradientCalculator->EvaluateAtIndex( movindex );
    for( unsigned int jj = 0; jj < ImageDimension; jj++ )
    {
    movingGradient[jj]*=(1.0);
    movingGradientSquaredMagnitude += vnl_math_sqr( movingGradient[jj] ) * m_FixedImageSpacing[j];
    } 
    }
  else
    {
    movingValue = 0.0;
    update.Fill(0.0);
    return update;
    }

  // Compute update
  
  //float ttt=1.e-2;
  //if (fixedValue > ttt) fixedValue=1.0-fixedValue; //else fixedValue=movingValue;
  
  double speedValue = fixedValue - movingValue;


  if (m_Symmetric) fixedGradient=fixedGradient+movingGradient*0.9;
  else if (m_MovingGradient) fixedGradient=movingGradient;

  if (m_Robust)
  {
    float m_RobustnessFactor=10.0*m_IntensityDifferenceThreshold;
    float frac=vnl_math_abs(speedValue)/m_RobustnessFactor;
    if ( frac > 1) speedValue/=frac;
//    if ( vnl_math_abs(speedValue) > 0.34) { update.Fill(0.0); return update; }
  }
  
  bool normalizemetric= Superclass::m_NormalizeGradient;  
  double denominator = 1.0;
  if (normalizemetric) 
  {   
//    denominator = (speedValue*speedValue/m_Normalizer*(fixedGradientSquaredMagnitude));
    denominator = fixedGradientSquaredMagnitude;
    denominator = sqrt(denominator);
  }

  denominator = vnl_math_sqr( speedValue ) /  m_Normalizer + 
    fixedGradientSquaredMagnitude;

  if (denominator == 0) denominator=1.0;


/* 
    denominator = (fixedGradientSquaredMagnitude);
    m_DenominatorThreshold= 1.e-11;
     
 denominator = 
    sqrt(fixedGradientSquaredMagnitude);*/

  if ( vnl_math_abs(fixedValue - movingValue) < m_IntensityDifferenceThreshold 
  ||   denominator < m_DenominatorThreshold )
    {
    for( j = 0; j < ImageDimension; j++ )
      {
      update[j] = 0.0;
      }
    return update;
    }

   Superclass::m_Energy+=speedValue*speedValue;
    
  for( j = 0; j < ImageDimension; j++ )
    {
    update[j] = (1.0)*speedValue * fixedGradient[j]* m_FixedImageSpacing[j] / 
      denominator* Superclass::m_GradientStep;//*exp(-0.2/fabs(denominator));
    }
//std::cout << " eno " << update  << " ind " << index << std::endl;

  return update;

}

 



} // end namespace itk

#endif

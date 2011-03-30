/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: itkLabelOverlapMeasuresImageFilter.txx,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDiReCTImageFilter_txx
#define __itkDiReCTImageFilter_txx

#include "itkDiReCTImageFilter.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkMaximumImageFilter.h"

namespace itk
{


template<class TInputImage, class TOutputImage>
DiReCTImageFilter<TInputImage, TOutputImage>
::DiReCTImageFilter() :
  m_MaximumNumberOfIterations( 50 ),
  m_ThicknessPriorEstimate( 6.0 ),
  m_SmoothingSigma( 1.0 ),
  m_GradientStep( 0.5 )
{
  this->SetNumberOfRequiredInputs( 3 );
}

template<class TInputImage, class TOutputImage>
DiReCTImageFilter<TInputImage, TOutputImage>
::~DiReCTImageFilter()
{
}

template<class TInputImage, class TOutputImage>
void
DiReCTImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  typedef BinaryThresholdImageFilter<InputImageType, InputImageType>
    ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( this->GetSegmentationImage() );
  thresholder->SetLowerThreshold( 3 );
  thresholder->SetUpperThreshold( 3 );
  thresholder->SetInsideValue( 1 );
  thresholder->SetOutsideValue( 0 );
  thresholder->Update();

  typename InputImageType::Pointer whiteMatter = InputImageType::New();
  whiteMatter = thresholder->GetOutput();
  whiteMatter->DisconnectPipeline();

  typedef CastImageFilter<InputImageType, RealImageType> CasterType;
  typename CasterType::Pointer caster = CasterType::New();
  caster->SetInput( whiteMatter );
  caster->Update();
  RealImagePointer laplacian = caster->GetOutput();
  laplacian->DisconnectPipeline();

  for( unsigned int n = 0; n < 100; n++ )
    {
    typedef DiscreteGaussianImageFilter<RealImageType, RealImageType> SmootherType;
    typename SmootherType::Pointer smoother = SmootherType::New();
    smoother->SetVariance( vnl_math_sqr( this->m_SmoothingSigma ) );
    smoother->SetUseImageSpacingOn();
    smoother->SetMaximumError( 0.01 );
    smoother->SetInput( laplacian );
    smoother->Update();

    typedef MaximumImageFilter<RealImageType, InputImageType, RealImageType>
      MaxFilterType;
    typename MaxFilterType::Pointer maxFilter = MaxFilterType::New();
    maxFilter->SetInput1( smoother->GetOutput() );
    maxFilter->SetInput2( whiteMatter );
    maxFilter->Update();
    laplacian = maxFilter->GetOutput();
    laplacian->DisconnectPipeline();
    }

  this->SetNthOutput( 0, laplacian );
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TOutputImage>
void
DiReCTImageFilter<TInputImage, TOutputImage>
::PrintSelf( std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  std::cout << indent << "Maximum number of iterations = "
    << this->m_MaximumNumberOfIterations << std::endl;
  std::cout << indent << "Thickness prior estimate = "
    << this->m_ThicknessPriorEstimate << std::endl;
  std::cout << indent << "Smoothing sigma = "
    << this->m_SmoothingSigma << std::endl;
  std::cout << indent << "Gradient step = "
    << this->m_GradientStep << std::endl;
}

} // end namespace itk

#endif

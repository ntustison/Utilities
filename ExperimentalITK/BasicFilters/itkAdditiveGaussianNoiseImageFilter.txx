/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAdditiveGaussianNoiseImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:49 $
  Version:   $Revision: 1.1.1.1 $
  Author:    Gavin Baker <gavinb@cs.mu.oz.au>

  Copyright (c) 2004 Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itkAdditiveGaussianNoiseImageFilter.h"

namespace itk
{

template <class TInputImage>
AdditiveGaussianNoiseImageFilter<TInputImage>
::AdditiveGaussianNoiseImageFilter()
{
  m_NoiseFilter = NoiseFilterType::New();
}


template <class TInputImage>
void
AdditiveGaussianNoiseImageFilter<TInputImage>
::GenerateData()
{
  // Setup grafted pipeline for composite filter

  m_NoiseFilter->SetInput( this->GetInput() );
  m_NoiseFilter->SetNumberOfThreads( 1 );
  m_NoiseFilter->Update();
  this->GraftOutput( m_NoiseFilter->GetOutput() );
}

template <class TInputImage>
void
AdditiveGaussianNoiseImageFilter<TInputImage>
::PrintSelf(std::ostream& os,
            Indent indent) const
{
  os
    << indent << "AdditiveGaussianNoiseImageFilter"
    << "\n  Mean: " << this->GetMean()
    << "\n  StandardDeviation: " << this->GetStandardDeviation()
    << std::endl;
}

} /* namespace itk */

/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImpulseNoiseImageFilter.hxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:52 $
  Version:   $Revision: 1.1.1.1 $
  Author:    Gavin Baker <gavinb@cs.mu.oz.au>

  Copyright (c) 2004 Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itkImpulseNoiseImageFilter.h"

namespace itk
{

template <class TInputImage>
ImpulseNoiseImageFilter<TInputImage>
::ImpulseNoiseImageFilter()
{
  m_NoiseFilter = NoiseFilterType::New();
}


template <class TInputImage>
void
ImpulseNoiseImageFilter<TInputImage>
::GenerateData()
{
  this->AllocateOutputs();

  // Setup grafted pipeline for composite filter

  m_NoiseFilter->SetInput( this->GetInput() );
  m_NoiseFilter->GraftOutput( this->GetOutput() );
  m_NoiseFilter->Update();
  this->GraftOutput( m_NoiseFilter->GetOutput() );
}

template <class TInputImage>
void
ImpulseNoiseImageFilter<TInputImage>
::PrintSelf(std::ostream& os,
            Indent indent) const
{
  os
    << indent << "ImpulseNoiseImageFilter"
    << "\n  Prob: " << this->GetProbability()
    << std::endl;
}

} /* namespace itk */

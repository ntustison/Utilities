/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGaborKernelFunction.cxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:20:03 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkGaborKernelFunction.h"

namespace itk
{

GaborKernelFunction::GaborKernelFunction()
{
  this->m_CalculateImaginaryPart = false;
  this->m_Sigma = 1.0;  
  this->m_Frequency = 0.4;
  this->m_PhaseOffset = 0.0;
}

GaborKernelFunction::~GaborKernelFunction()
{
}

} // namespace itk


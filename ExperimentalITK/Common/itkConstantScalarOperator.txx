/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkConstantScalarOperator.txx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:20:03 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkConstantScalarOperator_txx
#define _itkConstantScalarOperator_txx
#include "itkConstantScalarOperator.h"

namespace itk
{

template<class TPixel,unsigned int VDimension, class TAllocator>
typename ConstantScalarOperator<TPixel,VDimension, TAllocator>
::CoefficientVector
ConstantScalarOperator<TPixel,VDimension, TAllocator>
::GenerateCoefficients()
{
  unsigned long numberOfCoefficients = 1;
  for ( unsigned int i = 0; i < VDimension; i++ )
    {
    numberOfCoefficients *= this->GetRadius()[i];
    }

  CoefficientVector coeff( numberOfCoefficients, this->m_Scalar );
  return coeff;
}

}// end namespace itk

#endif

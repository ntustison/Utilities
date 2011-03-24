/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBinaryBoxStructuringElement.txx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:20:03 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBinaryBoxStructuringElement_txx
#define __itkBinaryBoxStructuringElement_txx
#include "itkBinaryBoxStructuringElement.h"

#include "itkNumericTraits.h"

namespace itk
{

// Create the structuring element
template <class TPixel, unsigned int VDimension, class TAllocator>
void
BinaryBoxStructuringElement<TPixel, VDimension, TAllocator>
::CreateStructuringElement()
{
  // 
  // Zero out the neighborhood
  //
  Iterator kernel_it;
  for (kernel_it=this->Begin(); kernel_it != this->End(); ++kernel_it)
    {
    *kernel_it = NumericTraits<TPixel>::One;
    }
}

} // namespace itk

#endif

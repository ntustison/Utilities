/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBinaryDiamondStructuringElement.hxx,v $
  Language:  C++
  Date:      $Date: 2008/12/03 18:37:17 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBinaryDiamondStructuringElement_hxx
#define __itkBinaryDiamondStructuringElement_hxx
#include "itkBinaryDiamondStructuringElement.h"

#include "itkNumericTraits.h"

#include "vnl/vnl_math.h"

namespace itk
{

// Create the structuring element
template <class TPixel, unsigned int VDimension, class TAllocator>
void
BinaryDiamondStructuringElement<TPixel, VDimension, TAllocator>
::CreateStructuringElement()
{

  unsigned int minRadius = this->GetRadius( 0 );
  for ( unsigned d = 1; d < NeighborhoodDimension; d++ )
    {
    if ( minRadius > this->GetRadius( d ) )
      {
      minRadius = this->GetRadius( d );
      }
    }

  for ( unsigned int n = 0; n < this->Size(); n++ )
    {
    OffsetType offset = this->GetOffset( n );
    unsigned int manhattanDistance = 0;
    for ( unsigned int d = 0; d < NeighborhoodDimension; d++ )
      {
      manhattanDistance += std::abs( offset[d] );
      }
    if ( manhattanDistance <= minRadius )
      {
      this->operator[]( n ) = NumericTraits<TPixel>::One;
      }
    else
      {
      this->operator[]( n ) = NumericTraits<TPixel>::Zero;
      }
    }
}

} // namespace itk

#endif

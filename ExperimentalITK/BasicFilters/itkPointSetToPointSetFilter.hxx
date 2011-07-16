/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPointSetToPointSetFilter.hxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:53 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkPointSetToPointSetFilter_hxx
#define _itkPointSetToPointSetFilter_hxx

#include "itkPointSetToPointSetFilter.h"


namespace itk
{
  
/**
 *
 */
template <class TInputPointSet, class TOutputPointSet>
PointSetToPointSetFilter<TInputPointSet,TOutputPointSet>
::PointSetToPointSetFilter()
{
  // Modify superclass default values, can be overridden by subclasses
  this->SetNumberOfRequiredInputs(1);

  this->ProcessObject::SetNumberOfRequiredOutputs(1);
}


/**
 *
 */
template <class TInputPointSet, class TOutputPointSet>
void 
PointSetToPointSetFilter<TInputPointSet,TOutputPointSet>
::SetInput(const TInputPointSet * pointSet)
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput(0, 
                                   const_cast< InputPointSetType * >( pointSet ) );
}

/**
 * Connect one of the operands for pixel-wise addition
 */
template <class TInputPointSet, class TOutputPointSet>
void
PointSetToPointSetFilter<TInputPointSet,TOutputPointSet>
::SetInput( unsigned int index, const TInputPointSet * pointSet ) 
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput(index, 
                                   const_cast< TInputPointSet *>( pointSet ) );
}


/**
 *
 */
template <class TInputPointSet, class TOutputPointSet>
const typename PointSetToPointSetFilter<TInputPointSet,TOutputPointSet>::InputPointSetType *
PointSetToPointSetFilter<TInputPointSet,TOutputPointSet>
::GetInput()
{
  if (this->GetNumberOfInputs() < 1)
    {
    return 0;
    }
  
  return static_cast<TInputPointSet*>
    (this->ProcessObject::GetInput(0));
}

  
/**
 *
 */
template <class TInputPointSet, class TOutputPointSet>
const typename PointSetToPointSetFilter<TInputPointSet,TOutputPointSet>::InputPointSetType *
PointSetToPointSetFilter<TInputPointSet,TOutputPointSet>
::GetInput(unsigned int idx)
{
  return static_cast<TInputPointSet*>
    (this->ProcessObject::GetInput(idx));
}


} // end namespace itk

#endif

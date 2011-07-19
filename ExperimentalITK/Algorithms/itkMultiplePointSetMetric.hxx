/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMultiplePointSetMetric.hxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:13:43 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkMultiplePointSetMetric_hxx
#define _itkMultiplePointSetMetric_hxx

#include "itkMultiplePointSetMetric.h"

namespace itk
{

/** Constructor */
template <class TPointSet>
MultiplePointSetMetric<TPointSet>
::MultiplePointSetMetric()
{
}

template <class TPointSet>
void
MultiplePointSetMetric<TPointSet>
::SetInput( const TPointSet *input )
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( 0,
                                    const_cast<PointSetType *>( input ) );
}

template <class TPointSet>
void
MultiplePointSetMetric<TPointSet>
::SetInput( unsigned int index, const TPointSet *points )
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( index,
                                    const_cast<TPointSet *>( points ) );
}

template <class TPointSet>
const typename MultiplePointSetMetric<TPointSet>::PointSetType*
MultiplePointSetMetric<TPointSet>
::GetInput( void )
{
  if( this->GetNumberOfInputs() < 1 )
    {
    return 0;
    }
  return static_cast<const TPointSet *>
    ( this->ProcessObject::GetInput( 0 ) );
}

template <class TPointSet>
const typename MultiplePointSetMetric<TPointSet>::PointSetType*
MultiplePointSetMetric<TPointSet>
::GetInput( unsigned int idx )
{
  return static_cast< const TPointSet * >
    ( this->ProcessObject::GetInput( idx ) );
}

/** Initialize the metric */
template <class TPointSet>
void
MultiplePointSetMetric<TPointSet>
::Initialize( void ) throw ( ExceptionObject )
{
  // If the PointSet is provided by a source, update the source.
  for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
    {
    if( this->GetInput( i )->GetSource() )
      {
      this->GetInput( i )->GetSource()->Update();
      }
    }
}


/** PrintSelf */
template <class TPointSet>
void
MultiplePointSetMetric<TPointSet>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}


} // end namespace itk

#endif


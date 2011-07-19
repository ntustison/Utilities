/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkInPlaceGraphFilter.hxx,v $
  Language:  C++
  Date:      $Date: 2008/11/11 03:08:24 $
  Version:   $Revision: 1.2 $

  Copyright ( c ) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkInPlaceGraphFilter_hxx
#define _itkInPlaceGraphFilter_hxx
#include "itkInPlaceGraphFilter.h"


namespace itk
{

/**
 *
 */
template <class TInputGraph, class TOutputGraph>
InPlaceGraphFilter<TInputGraph, TOutputGraph>
::InPlaceGraphFilter()
  : m_InPlace( true )
{
}

/**
 *
 */
template <class TInputGraph, class TOutputGraph>
InPlaceGraphFilter<TInputGraph, TOutputGraph>
::~InPlaceGraphFilter()
{
}



template<class TInputGraph, class TOutputGraph>
void
InPlaceGraphFilter<TInputGraph, TOutputGraph>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "InPlace: " 
     << ( this->m_InPlace ? "On" : "Off" ) << std::endl;
  if ( this->CanRunInPlace() )
    {
    os << indent << "The input and output to this filter are the same type. "
       << "The filter can be run in place." << std::endl;
    }
  else
    {
    os << indent << "The input and output to this filter are different types. "
       << "The filter cannot be run in place." << std::endl;
    }
}

template<class TInputGraph, class TOutputGraph>
void
InPlaceGraphFilter<TInputGraph, TOutputGraph>
::AllocateOutputs()
{
  // if told to run in place and the types support it,
  if ( this->m_InPlace && ( typeid( TInputGraph ) == typeid( TOutputGraph ) ) )
    {
    // Graft this first input to the output.  Later, we'll need to
    // remove the input's hold on the bulk data.
    //
    OutputGraphPointer inputAsOutput = dynamic_cast<TOutputGraph *>(
      const_cast<TInputGraph *>( this->GetInput() ) );
    if ( inputAsOutput )
      {
      this->GraftOutput( inputAsOutput );
      }
    else
      {
      // if we cannot cast the input to an output type, then allocate
      // an output usual.

      OutputGraphPointer outputPtr;

      outputPtr = this->GetOutput( 0 );
      outputPtr = TOutputGraph::New();
      }

    // If there are more than one outputs, allocate the remaining outputs
    for ( unsigned int i = 1; i < this->GetNumberOfOutputs(); i++ )
      {
      OutputGraphPointer outputPtr;

      outputPtr = this->GetOutput( i );
      outputPtr = TOutputGraph::New();
      }
    }
//  else
//    {
//    Superclass::AllocateOutputs();
//    }
}

template<class TInputGraph, class TOutputGraph>
void
InPlaceGraphFilter<TInputGraph, TOutputGraph>
::ReleaseInputs()
{
  // if told to run in place and the types support it,
  if ( this->m_InPlace && ( typeid( TInputGraph ) == typeid( TOutputGraph ) ) )
    {
    // Release any input where the ReleaseData flag has been set
    ProcessObject::ReleaseInputs();

    // Release input 0 by default since we overwrote it
    TInputGraph * ptr = const_cast<TInputGraph*>( this->GetInput() );
    if( ptr )
      {
      ptr->ReleaseData();
      }
    }
  else
    {
    Superclass::ReleaseInputs();
    }
}


} // end namespace itk

#endif

/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: itkGraphDuplicator.txx,v $
Language:  C++
Date:      $Date: 2008/11/11 03:08:24 $
Version:   $Revision: 1.2 $

Copyright ( c ) Insight Software Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkGraphDuplicator_txx
#define _itkGraphDuplicator_txx

#include "itkGraphDuplicator.h"

namespace itk
{

/** Constructor */
template<class TInputGraph>
GraphDuplicator<TInputGraph>
::GraphDuplicator()
{
  this->m_Input = NULL;
  this->m_Output = NULL;
  this->m_InternalGraphTime = 0;
}

/** Update function */
template<class TInputGraph>
void
GraphDuplicator<TInputGraph>
::Update( void )
{
  if( !this->m_Input )
    {
    itkExceptionMacro( "Input Graph has not been connected" );
    return;
    }

  // Update only if the input Graph has been modified
  unsigned long t, t1, t2;
  t1 = this->m_Input->GetPipelineMTime();
  t2 = this->m_Input->GetMTime();
  t = ( t1 > t2 ? t1 : t2 );

  if( t == this->m_InternalGraphTime )
    {
    return; // No need to update
    }

  // Cache the timestamp
  this->m_InternalGraphTime = t;

  // Allocate the graph
  this->m_Output = GraphType::New();
  this->m_Output->Graft( this->m_Input );
}

template<class TInputGraph>
void
GraphDuplicator<TInputGraph>
::PrintSelf(  std::ostream& os, Indent indent  ) const
{
  Superclass::PrintSelf( os,indent );
  os << indent << "Input Graph: " << this->m_Input << std::endl;
  os << indent << "Output Graph: " << this->m_Output << std::endl;
  os << indent << "Internal Graph Time: "
     << this->m_InternalGraphTime << std::endl;
}

} // end namespace itk

#endif

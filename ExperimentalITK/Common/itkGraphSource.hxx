/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGraphSource.hxx,v $
  Language:  C++
  Date:      $Date: 2008/11/11 03:08:24 $
  Version:   $Revision: 1.2 $

  Copyright ( c ) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkGraphSource_hxx
#define _itkGraphSource_hxx

#include "itkGraphSource.h"

namespace itk
{

/**
 *
 */
template<class TOutputGraph>
GraphSource<TOutputGraph>
::GraphSource()
{
  // Create the output. We use static_cast<> here because we know the default
  // output must be of type TOutputGraph
  OutputGraphPointer output
    = static_cast<TOutputGraph*>( this->MakeOutput( 0 ).GetPointer() );

  this->ProcessObject::SetNumberOfRequiredOutputs( 1 );
  this->ProcessObject::SetNthOutput( 0, output.GetPointer() );

  m_GenerateDataRegion = 0;
  m_GenerateDataNumberOfRegions = 0;
}

/**
 *
 */
template<class TOutputGraph>
typename GraphSource<TOutputGraph>::DataObjectPointer
GraphSource<TOutputGraph>
::MakeOutput( unsigned int )
{
  return static_cast<DataObject*>( TOutputGraph::New().GetPointer() );
}

/**
 *
 */
template<class TOutputGraph>
typename GraphSource<TOutputGraph>::OutputGraphType *
GraphSource<TOutputGraph>
::GetOutput( void )
{
  if ( this->GetNumberOfOutputs() < 1 )
    {
    return 0;
    }

  return static_cast<TOutputGraph*>
    ( this->ProcessObject::GetOutput( 0 ) );
}


/**
 *
 */
template<class TOutputGraph>
typename GraphSource<TOutputGraph>::OutputGraphType *
GraphSource<TOutputGraph>
::GetOutput( unsigned int idx )
{
  return static_cast<TOutputGraph*>
    ( this->ProcessObject::GetOutput( idx ) );
}


/**
 *
 */
template<class TOutputGraph>
void
GraphSource<TOutputGraph>
::SetOutput( OutputGraphType *output )
{
  this->ProcessObject::SetNthOutput( 0, output );
}


/**
 *
 */
template<class TOutputGraph>
void
GraphSource<TOutputGraph>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();
}

/**
 *
 */
template<class TOutputGraph>
void
GraphSource<TOutputGraph>
::GraftOutput( OutputGraphType *graft )
{
  OutputGraphType * output = this->GetOutput();

  if ( output && graft )
    {
    // grab a handle to the bulk data of the specified data object

    // copy the meta-information

    output->Graft( graft );
    }
}

/**
 *
 */
template<class TOutputGraph>
void
GraphSource<TOutputGraph>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os,indent );
}

} // end namespace itk

#endif

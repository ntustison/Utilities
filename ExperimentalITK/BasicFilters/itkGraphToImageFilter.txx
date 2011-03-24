/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGraphToImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2008/11/11 03:08:13 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGraphToImageFilter_txx
#define __itkGraphToImageFilter_txx

#include "itkGraphToImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkNumericTraits.h"

namespace itk
{

/** Constructor */
template <class TInputGraph, class TOutputImage>
GraphToImageFilter<TInputGraph,TOutputImage>
::GraphToImageFilter()
{
  this->SetNumberOfRequiredInputs( 1 );
  this->m_Size.Fill( 0 );

  for( unsigned int i = 0; i < OutputImageDimension; i++ )
    {
    // Set an image spacing for the user
    m_Spacing[i] = 1.0;
    m_Origin[i] = 0;
    }

  m_BackgroundValue = NumericTraits<ValueType>::Zero;
}

/** Destructor */
template <class TInputGraph, class TOutputImage>
GraphToImageFilter<TInputGraph,TOutputImage>
::~GraphToImageFilter()
{
}


/** Set the Input SpatialObject */
template <class TInputGraph, class TOutputImage>
void
GraphToImageFilter<TInputGraph,TOutputImage>
::SetInput( const InputGraphType *input )
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( 0,
    const_cast< InputGraphType * >( input ) );
}


/** Connect one of the operands  */
template <class TInputGraph, class TOutputImage>
void
GraphToImageFilter<TInputGraph,TOutputImage>
::SetInput( unsigned int index, const InputGraphType * graph )
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( index,
    const_cast< InputGraphType *>( graph ) );
}



/** Get the input Graph */
template <class TInputGraph, class TOutputImage>
const typename GraphToImageFilter<TInputGraph,TOutputImage>::InputGraphType *
GraphToImageFilter<TInputGraph,TOutputImage>
::GetInput( void )
{
  if( this->GetNumberOfInputs() < 1 )
    {
    return 0;
    }

  return static_cast<const TInputGraph * >
    ( this->ProcessObject::GetInput( 0 ) );
}

/** Get the input Graph */
template <class TInputGraph, class TOutputImage>
const typename GraphToImageFilter<TInputGraph,TOutputImage>::InputGraphType *
GraphToImageFilter<TInputGraph,TOutputImage>
::GetInput( unsigned int idx )
{
  return static_cast< const TInputGraph * >
    ( this->ProcessObject::GetInput( idx ) );
}

//----------------------------------------------------------------------------
template <class TInputGraph, class TOutputImage>
void
GraphToImageFilter<TInputGraph,TOutputImage>
::SetSpacing( const double spacing[Self::OutputImageDimension] )
{
  unsigned int i;
  for( i = 0; i < OutputImageDimension; i++ )
    {
    if( spacing[i] != this->m_Spacing[i] )
      {
      break;
      }
    }
  if( i < OutputImageDimension )
    {
    for( i = 0; i < OutputImageDimension; i++ )
      {
      this->m_Spacing[i] = spacing[i];
      }
    }
}

template <class TInputGraph, class TOutputImage>
void
GraphToImageFilter<TInputGraph,TOutputImage>
::SetSpacing( const float spacing[Self::OutputImageDimension] )
{
  unsigned int i;
  for( i = 0; i < OutputImageDimension; i++ )
    {
    if( static_cast<double>( spacing[i] ) != this->m_Spacing[i] )
      {
      break;
      }
    }
  if( i < OutputImageDimension )
    {
    for( i = 0; i < OutputImageDimension; i++ )
      {
      this->m_Spacing[i] = spacing[i];
      }
    }
}

template <class TInputGraph, class TOutputImage>
const double *
GraphToImageFilter<TInputGraph,TOutputImage>
::GetSpacing() const
{
  return this->m_Spacing;
}

//----------------------------------------------------------------------------
template <class TInputGraph, class TOutputImage>
void
GraphToImageFilter<TInputGraph,TOutputImage>
::SetOrigin( const double origin[Self::OutputImageDimension] )
{
  unsigned int i;
  for( i = 0; i < OutputImageDimension; i++ )
    {
    if( origin[i] != this->m_Origin[i] )
      {
      break;
      }
    }
  if( i < OutputImageDimension )
    {
    for( i = 0; i < OutputImageDimension; i++ )
      {
      this->m_Origin[i] = origin[i];
      }
    }
}

template <class TInputGraph, class TOutputImage>
void
GraphToImageFilter<TInputGraph,TOutputImage>
::SetOrigin( const float origin[Self::OutputImageDimension] )
{
  unsigned int i;
  for( i = 0; i < OutputImageDimension; i++ )
    {
    if( static_cast<double>( origin[i] ) != this->m_Origin[i] )
      {
      break;
      }
    }
  if( i < OutputImageDimension )
    {
    for( i = 0; i < OutputImageDimension; i++ )
      {
      this->m_Origin[i] = origin[i];
      }
    }
}

template <class TInputGraph, class TOutputImage>
void
GraphToImageFilter<TInputGraph,TOutputImage>
::GenerateData( void )
{
  unsigned int i;

  OutputImagePointer OutputImage = this->GetOutput();

  // Generate the image
  double origin[OutputImageDimension];
  SizeType size;

  size.Fill(0);
  for( i = 0; i < OutputImageDimension; i++ )
    {
    origin[i] = 0;
    }

  typename OutputImageType::IndexType index;
  index.Fill( 0 );
  typename OutputImageType::RegionType region;

  bool specified = false;
  for( i = 0; i < OutputImageDimension; i++ )
    {
    if( this->m_Size[i] != 0 )
      {
      specified = true;
      break;
      }
    }

  if( specified )
    {
    region.SetSize( this->m_Size );
    }
  else
    {
    itkExceptionMacro( "Currently, the user MUST specify an image size" );
    //region.SetSize( size );
    }
  region.SetIndex( index );

  OutputImage->SetLargestPossibleRegion( region );    //
  OutputImage->SetBufferedRegion( region );           // set the region
  OutputImage->SetRequestedRegion( region );          //

  // If the spacing has been explicitly specified, the filter
  // will set the output spacing to that explicit spacing,
  // otherwise the spacing from the spatial object is used as default.

  specified = false;
  for( i = 0; i < OutputImageDimension; i++ )
    {
    if( this->m_Spacing[i] != 0 )
      {
      specified = true;
      break;
      }
    }

  if( specified )
    {
    OutputImage->SetSpacing( this->m_Spacing );         // set spacing
    }
  else
    {
    itkExceptionMacro( "Currently, the user MUST specify an image spacing" )
    //OutputImage->SetSpacing(
    //  InputObject->GetIndexToObjectTransform()->GetScaleComponent() );
    }
  OutputImage->SetOrigin( origin );   //   and origin
  OutputImage->Allocate();   // allocate the image

  ImageRegionIterator<OutputImageType> imageIt( OutputImage, region );
  for( imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt )
    {
    imageIt.Set( this->m_BackgroundValue );
    }

  itkDebugMacro( "GraphToImageFilter::GenerateData() finished" );

}


template <class TInputGraph, class TOutputImage>
void
GraphToImageFilter<TInputGraph,TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Size : " << m_Size << std::endl;
  os << indent << "Background Value : " << m_BackgroundValue << std::endl;
}



} // end namespace itk

#endif

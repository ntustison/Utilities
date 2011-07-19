/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMinimalPathImageFunction.hxx,v $
  Language:  C++
  Date:      $Date: 2009/03/21 01:33:09 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkMinimalPathImageFunction_hxx
#define __itkMinimalPathImageFunction_hxx

#include "itkMinimalPathImageFunction.h"

#include "itkCastImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhood.h"

#include "vnl/vnl_math.h"

namespace itk
{

// Constructor
template<class TInputImage>
MinimalPathImageFunction<TInputImage>
::MinimalPathImageFunction()
{
  this->m_UseFaceConnectedness = true;
  this->m_UseImageSpacing = true;

  this->m_MaskImage = NULL;
  this->m_InsideMaskPixelValue = NumericTraits<MaskPixelType>::One;

  this->m_AnchorSeed.Fill( 0 );
}
 
// Destructor
template<class TInputImage>
MinimalPathImageFunction<TInputImage>
::~MinimalPathImageFunction()
{
}

/**
 * Initialize by setting the input image
 */
template <class TInputImage>
void
MinimalPathImageFunction<TInputImage>
::SetInputImage( const InputImageType * ptr )
{
  Superclass::SetInputImage( ptr ); 
  this->GeneratePathDirectionImage();
}

template <class TInputImage>
void
MinimalPathImageFunction<TInputImage>
::GeneratePathDirectionImage()
{
  if ( !this->IsInsideBuffer( this->m_AnchorSeed ) )
    {
    itkExceptionMacro( "m_AnchorSeed is not inside buffer." );
    }

  /**
   * Initialize data structures
   */
  typedef Image<bool, ImageDimension> BooleanImageType;
  typename BooleanImageType::Pointer expanded = BooleanImageType::New();
  expanded->SetOrigin( this->GetInputImage()->GetOrigin() );
  expanded->SetSpacing( this->GetInputImage()->GetSpacing() );
  expanded->SetRegions( this->GetInputImage()->GetRequestedRegion() );
  expanded->Allocate();
  expanded->FillBuffer( false );
  expanded->SetPixel( this->m_AnchorSeed, true );

  this->m_PathDirectionImage = OffsetImageType::New();
  this->m_PathDirectionImage->SetOrigin( this->GetInputImage()->GetOrigin() );
  this->m_PathDirectionImage->SetSpacing( this->GetInputImage()->GetSpacing() );
  this->m_PathDirectionImage->SetRegions( this->GetInputImage()->GetRequestedRegion() );
  this->m_PathDirectionImage->Allocate();

  typename ConstNeighborhoodIterator<InputImageType>::RadiusType radius;
  radius.Fill( 1 );
  ConstNeighborhoodIterator<BooleanImageType> It( radius, expanded, 
    expanded->GetRequestedRegion() );
  unsigned int numberOfNeighbors = It.GetNeighborhood().Size();

  Neighborhood<RealType, ImageDimension> scaleFactors;
  scaleFactors.SetRadius( radius );
  for ( unsigned int n = 0; n < numberOfNeighbors; n++ )
    {
    scaleFactors[n] = 0;
    if ( n == static_cast<unsigned int>( 0.5 * numberOfNeighbors ) )
      {
      continue; 
      } 
    typename Neighborhood<RealType, ImageDimension>::OffsetType offset 
      = scaleFactors.GetOffset( n );
   
				bool isFaceConnected = true;
				unsigned int sumOffset = 0; 
				for ( unsigned int d = 0; d < ImageDimension; d++ )
						{
						sumOffset += vnl_math_abs( offset[d] );
						if ( this->m_UseImageSpacing )
								{
								scaleFactors[n] += ( static_cast<RealType>( offset[d] * offset[d] ) 
										* this->GetInputImage()->GetSpacing()[d] 
          * this->GetInputImage()->GetSpacing()[d] );  
								}
						else
								{
								scaleFactors[n] += static_cast<RealType>( offset[d] * offset[d] );  
								}
						if ( sumOffset > 1 && this->m_UseFaceConnectedness )
								{
								isFaceConnected = false;
								break;
								}
      }      
				if ( !isFaceConnected )
						{
						scaleFactors[n] = 0;
						}
    if ( scaleFactors[n] > 0 )
      {
      scaleFactors[n] = vcl_sqrt( scaleFactors[n] );
      } 
    }   
  /**
   * Generate the PathDirectionImage with Dijkstra's algorithm
   */

  PriorityQueueElementType anchorElement( this->m_AnchorSeed, 0.0 );

  if ( this->m_MaskImage && 
       this->m_MaskImage->GetPixel( this->m_AnchorSeed ) 
         != this->m_InsideMaskPixelValue )
    {
    itkWarningMacro( "The anchor seed is outside the user-defined mask region." ); 
    return;
    }    

  typename PriorityQueueType::Pointer Q = PriorityQueueType::New();
  Q->Initialize();
  Q->Push( anchorElement );

  while ( !Q->Empty() )
    {
    PriorityQueueElementType centerElement = Q->Peek(); 
    Q->Pop();

    expanded->SetPixel( centerElement.m_Element, true );
    It.SetLocation( centerElement.m_Element ); 
    PointType centerPoint;
    this->GetInputImage()->TransformIndexToPhysicalPoint( 
      centerElement.m_Element, centerPoint );

    for ( unsigned int n = 0; n < numberOfNeighbors; n++ )
      {
      if ( scaleFactors[n] == 0 )
        {
        continue;
        } 
      bool inBounds;
      bool isExpanded = It.GetPixel( n, inBounds );
      if ( isExpanded || !inBounds )
        {
        continue;
        }  

						IndexType neighborIndex = It.GetIndex( n );

						if ( this->m_MaskImage && this->m_MaskImage->GetPixel( 
						  neighborIndex ) != this->m_InsideMaskPixelValue )
								{
								continue;
								}    

						RealType neighborCost 
						  = centerElement.m_Priority + scaleFactors[n] *
        this->GetInputImage()->GetPixel( neighborIndex );

					 PriorityQueueElementType neighborElement( neighborIndex, neighborCost );

						Q->Push( neighborElement ); 

						this->m_PathDirectionImage->SetPixel( neighborIndex, 
								centerElement.m_Element - neighborElement.m_Element ); 

						PriorityQueueElementType element = Q->Peek();
						while ( expanded->GetPixel( element.m_Element ) )
								{
								Q->Pop();
								element = Q->Peek();
								} 
      }  
    }  
}

template <class TInputImage>
typename MinimalPathImageFunction<TInputImage>::OutputType::Pointer 
MinimalPathImageFunction<TInputImage>
::EvaluateAtIndex( const IndexType &index ) const
{ 
  if ( !this->IsInsideBuffer( index ) )
    {
    itkWarningMacro( "Requested index is not inside buffer." );
    return NULL;
    }

  if ( this->m_MaskImage && 
       this->m_MaskImage->GetPixel( index ) 
         != this->m_InsideMaskPixelValue )
    {
    itkWarningMacro( "The index is outside the user-defined mask region." ); 
    return NULL;
    }    

  typename OutputType::Pointer output = OutputType::New();
  output->Initialize();

  IndexType currentIndex = index;
  while ( currentIndex != this->m_AnchorSeed )
    {
    output->AddVertex( VertexType( currentIndex ) );
    currentIndex += this->m_PathDirectionImage->GetPixel( currentIndex );
    } 
  output->AddVertex( VertexType( currentIndex ) );

  return output;
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputImage>
void
MinimalPathImageFunction<TInputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "AnchorSeed: " 
     << this->m_AnchorSeed << std::endl;
  os << indent << "UseImageSpacing: " 
     << this->m_UseImageSpacing << std::endl;
  os << indent << "UseFaceConnectedness: " 
     << this->m_UseFaceConnectedness << std::endl;

  os << indent << "MaskImage"
     << this->m_MaskImage << std::endl;
  os << indent << "InsideMaskPixelValue"
     << this->m_InsideMaskPixelValue << std::endl;
  
}
  
} // end namespace itk

#endif

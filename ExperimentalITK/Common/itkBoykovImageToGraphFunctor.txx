/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBoykovImageToGraphFunctor.txx,v $
  Language:  C++
  Date:      $Date: 2008/11/11 03:08:24 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkBoykovImageToGraphFunctor_txx
#define _itkBoykovImageToGraphFunctor_txx

#include "itkBoykovImageToGraphFunctor.h"
#include "vnl/vnl_math.h"

namespace itk
{

template<typename TInputImage, typename TOutputGraph>
BoykovImageToGraphFunctor<TInputImage, TOutputGraph>
::BoykovImageToGraphFunctor()
{
  this->SetNumberOfInputs( 3 );

  m_SourceIndexContainer.clear();
  m_SinkIndexContainer.clear();

  m_Lambda = 1.0;
  m_Sigma = 10.0;
}

template<typename TInputImage, typename TOutputGraph>
typename BoykovImageToGraphFunctor<TInputImage, TOutputGraph>::NodeWeightType
BoykovImageToGraphFunctor<TInputImage, TOutputGraph>
::GetSourceDataTerm( IndexType idx )
{
  typename LikelihoodImageType::Pointer source =
     const_cast<LikelihoodImageType*>( static_cast<const LikelihoodImageType*>
        ( this->ProcessObject::GetInput( 1 ) ) );
  if ( !source )
    {
    itkExceptionMacro( "Likelihood image not specified.");
    }

  RealType q = vnl_math_max( source->GetPixel( idx ), 1e-10 );
  return static_cast<NodeWeightType>( -m_Lambda * vcl_log( q ) );
}

template<typename TInputImage, typename TOutputGraph>
typename BoykovImageToGraphFunctor<TInputImage, TOutputGraph>::NodeWeightType
BoykovImageToGraphFunctor<TInputImage, TOutputGraph>
::GetSinkDataTerm( IndexType idx )
{
  typename LikelihoodImageType::Pointer sink =
    const_cast<LikelihoodImageType*>(static_cast<const LikelihoodImageType*>
    ( this->ProcessObject::GetInput( 2 ) ) );

  if ( !sink )
    {
    itkExceptionMacro( << "Likelihood image not specified.");
    }

  RealType r = vnl_math_max( sink->GetPixel( idx ), 1e-10 );
  return static_cast<NodeWeightType>( -m_Lambda * vcl_log( r ) );
}

template<typename TInputImage, typename TOutputGraph>
typename BoykovImageToGraphFunctor<TInputImage, TOutputGraph>::EdgeWeightType
BoykovImageToGraphFunctor<TInputImage, TOutputGraph>
::GetEdgeWeight( IndexType idx1, IndexType idx2 )
{
  typename InputImageType::Pointer input =
     const_cast<InputImageType*>( this->GetInput( 0 ) );

  RealType q = static_cast<RealType>(input->GetPixel( idx1 ) );
  RealType r = static_cast<RealType>(input->GetPixel( idx2 ) );

  RealType PixelDistance = 0.0;
  for ( unsigned int d = 0; d < InputImageType::ImageDimension; d++ )
    {
    PixelDistance += vnl_math_sqr( static_cast<RealType>( idx1[d]-idx2[d] ) );
    }
  return static_cast<EdgeWeightType>(
    static_cast<RealType>( this->GetNumberOfPixelsInNeighborhood() )
    * vcl_exp( -0.5* vnl_math_sqr( (q-r) / m_Sigma ) ) /sqrt( PixelDistance )
  );
}

template<typename TInputImage, typename TOutputGraph>
void
BoykovImageToGraphFunctor<TInputImage, TOutputGraph>
::NormalizeGraph( NodeImageType *image, OutputGraphType *graph )
{
  /** Find the reverse edges */
  graph->SetAllReverseEdges();

  if ( this->m_SourceIndexContainer.size() > 0 ||
    this->m_SinkIndexContainer.size() > 0 )
    {
    typename EdgeIdentifierContainerType::const_iterator it;

    /** Calculate K */
    NodeWeightType K = static_cast<NodeWeightType>( 0 );
    NodePointerType node;
    NodeIteratorType NIt(graph);
    for ( NIt.GoToBegin(); !NIt.IsAtEnd(); ++NIt )
      {
      node = NIt.GetPointer();
      NodeWeightType B = static_cast<NodeWeightType>( 0 );
      EdgeIdentifierContainerType edges = node->OutgoingEdges;
      for (it = edges.begin(); it != edges.end(); ++it )
        {
        B += graph->GetEdgeWeight(*it);
        }
      if ( K < B + static_cast<NodeWeightType>( 1 ) )
        {
        K = B + static_cast<NodeWeightType>(1);
        }
      }

    /** Calculate and reassign the node weights based on whether the pixel
      * is hard-constrained as a "sink" or a "source" pixel.
      */

    typename IndexContainerType::const_iterator It;
    if ( this->m_SourceIndexContainer.size() > 0 )
      {
      for ( It = this->m_SourceIndexContainer.begin();
        It != this->m_SourceIndexContainer.end(); ++It )
        {
        if ( this->IsPixelANode( *It ) )
          {
          NodePointerType node = graph->GetNodePointer(
            image->GetPixel( *It ) );
          graph->SetNodeWeight( node, K );
          }
        }
      }
    if ( this->m_SinkIndexContainer.size() > 0 )
      {
      for ( It = this->m_SinkIndexContainer.begin();
        It != this->m_SinkIndexContainer.end(); ++It )
        {
        if ( this->IsPixelANode( *It ) )
          {
          NodePointerType node
            = graph->GetNodePointer( image->GetPixel( *It ) );
          graph->SetNodeWeight( node, -K );
          }
        }
      }
    }
}

template<typename TInputImage, typename TOutputGraph>
void
BoykovImageToGraphFunctor<TInputImage, TOutputGraph>
::SetSourceLikelihoodImage( const LikelihoodImageType *input )
{
  this->ProcessObject::SetNthInput( 1,
    const_cast<LikelihoodImageType *>( input ) );
}

template<typename TInputImage, typename TOutputGraph>
void
BoykovImageToGraphFunctor<TInputImage, TOutputGraph>
::SetSinkLikelihoodImage( const LikelihoodImageType *input )
{
  this->ProcessObject::SetNthInput( 2,
    const_cast<LikelihoodImageType *>( input ) );
}

template<typename TInputImage, typename TOutputGraph>
void
BoykovImageToGraphFunctor<TInputImage, TOutputGraph>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Lambda = " << m_Lambda << std::endl;
  os << indent << "Sigma = " << m_Sigma << std::endl;
}


} // end namespace itk

#endif

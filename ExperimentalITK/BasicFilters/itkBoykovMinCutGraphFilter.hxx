/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBoykovMinCutGraphFilter.hxx,v $
  Language:  C++
  Date:
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBoykovMinCutGraphFilter_hxx_
#define __itkBoykovMinCutGraphFilter_hxx_

#include "itkBoykovMinCutGraphFilter.h"
#include "itkGraphToGraphFilter.h"
#include "vnl/vnl_numeric_traits.h"

namespace itk
{

template <class TGraph>
BoykovMinCutGraphFilter<TGraph>
::BoykovMinCutGraphFilter()
{
  this->m_WeightZero = static_cast<WeightType>( 0 );

  this->m_TerminalEdge = NULL;
  this->m_OrphanEdge = NULL;
}

/** Generate the data */
template <class TGraph>
void
BoykovMinCutGraphFilter<TGraph>
::GenerateData()
{
  this->GenerateMinCut();
}

template <class TGraph>
void
BoykovMinCutGraphFilter<TGraph>
::GenerateMinCut()
{
  NodePointerType i, j, orphan, node = NULL;
  EdgePointerType edge = NULL;
  typename EdgeIdentifierContainerType::const_iterator it;

  this->Initialize();

  while( true )
    {
    i = node;

    if( i != NULL )
      {
      /** remove active flag */
      i->IsActive = false;
      if( i->Parent == NULL )
        {
        i = NULL;
        }
      }
    if( i == NULL )
      {
      i = this->GetNextActiveNode();
      if( i == NULL )
        {
        break;
        }
      }

    /** Growing step */
    if( !i->IsSink )
      {
      /* Grow source tree **/
      for( it = i->OutgoingEdges.begin(); it != i->OutgoingEdges.end(); ++it )
        {
        edge = this->m_Output->GetEdgePointer( *it );
        if( !this->IsEqual( this->m_Output->GetEdgeWeight( edge ),
          this->m_WeightZero ) )
          {
          j = this->m_Output->GetTargetNodePointer( edge );
          if( j->Parent == NULL )
            {
            j->IsSink = false;
            j->Parent = this->m_Output->GetReverseEdgePointer( edge );
            j->TimeStamp = i->TimeStamp;
            j->DistanceToTerminal = i->DistanceToTerminal + 1;
            this->SetActiveNode( j );
            }
          else if( j->IsSink )
            {
            break;
            }
          else if( j->TimeStamp <= i->TimeStamp &&
            j->DistanceToTerminal > i->DistanceToTerminal )
            {
            /*
             * heuristic - trying to shorten the distance from j to the source
             **/
            j->Parent = this->m_Output->GetReverseEdgePointer( edge );
            j->TimeStamp = i->TimeStamp;
            j->DistanceToTerminal = i->DistanceToTerminal + 1;
            }
          }
        }
      }
    else
      {
    /* Grow sink tree **/
      for( it = i->OutgoingEdges.begin(); it != i->OutgoingEdges.end(); ++it )
        {
        edge = this->m_Output->GetEdgePointer(*it);
        if( !this->IsEqual( this->m_Output->GetEdgeWeight(
            edge->ReverseEdgeIdentifier ), this->m_WeightZero ) )
          {
          j = this->m_Output->GetNodePointer( edge->TargetIdentifier );
          if( j->Parent == NULL )
            {
            j->IsSink = true;
            j->Parent = this->m_Output->GetEdgePointer(
              edge->ReverseEdgeIdentifier );
            j->TimeStamp = i->TimeStamp;
            j->DistanceToTerminal = i->DistanceToTerminal + 1;
            this->SetActiveNode( j );
            }
          else if( !j->IsSink )
            {
            edge = this->m_Output->GetReverseEdgePointer( edge );
            break;
            }
          else if( j->TimeStamp <= i->TimeStamp &&
                   j->DistanceToTerminal > i->DistanceToTerminal )
            {
            /* heuristic - try to shorten the distance from j to the source **/
            j->Parent = this->m_Output->GetReverseEdgePointer( edge );
            j->TimeStamp = i->TimeStamp;
            j->DistanceToTerminal = i->DistanceToTerminal + 1;
            }
          }
        }
      }

    this->m_GlobalTime++;

    if( it != i->OutgoingEdges.end() )
      {
      /** set active flag */
      i->IsActive = true;
      node = i;

      /** Augmentation step */
      this->Augment( edge );

      /** Adoption step */
      while( !this->m_Orphans.empty() )
        {
        orphan = this->m_Orphans.back();
        this->m_Orphans.pop_back();
        if( orphan->IsSink )
          {
          this->ProcessSinkOrphan( orphan );
          }
        else
          {
          this->ProcessSourceOrphan( orphan );
          }
        }
      }
    else
      {
      node = NULL;
      }
    }
}

template <class TGraph>
void
BoykovMinCutGraphFilter<TGraph>
::Initialize()
{
  this->AllocateOutputs();
  this->m_Output = this->GetOutput();

  this->m_ActiveNodes.clear();
  this->m_Orphans.clear();

  if( !this->m_TerminalEdge )
    {
    this->m_TerminalEdge = this->m_Output->CreateNewEdge();
    }
  if( !this->m_OrphanEdge )
    {
    this->m_OrphanEdge = this->m_Output->CreateNewEdge();
    }
  this->m_GlobalTime = 0;
  this->m_MaxFlow = 0;

  /** Set other node parameters */
  NodeIteratorType It( this->m_Output );
  NodePointerType node;
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    node = It.GetPointer();
    node->IsActive = false;
    node->TimeStamp = 0;

    /* node is connected to the source **/
    if( this->m_Output->GetNodeWeight( node ) > this->m_WeightZero )
      {
      node->IsSink = false;
      node->Parent = this->m_TerminalEdge;
      this->SetActiveNode( node );
      node->DistanceToTerminal = 1;
      }
    /* node is connected to the sink **/
    else if( this->m_Output->GetNodeWeight( node ) < this->m_WeightZero )
      {
      node->IsSink = true;
      node->Parent = this->m_TerminalEdge;
      this->SetActiveNode( node );
      node->DistanceToTerminal = 1;
      }
    else
      {
      node->Parent = NULL;
      }
    }
}

template <class TGraph>
void
BoykovMinCutGraphFilter<TGraph>
::SetActiveNode(NodePointerType node)
{
  if( !node->IsActive )
    {
    node->IsActive = true;
    this->m_ActiveNodes.push_back( node );
    }
}

template <class TGraph>
typename BoykovMinCutGraphFilter<TGraph>::NodePointerType
BoykovMinCutGraphFilter<TGraph>
::GetNextActiveNode()
{
  NodePointerType node;

  while( true )
    {
    if( this->m_ActiveNodes.empty() )
      {
      return NULL;
      }

    /* remove the node from the active list **/
    node = this->m_ActiveNodes.front();
    node->IsActive = false;
    this->m_ActiveNodes.pop_front();

    /* a node is active iff it has a Parent **/
    if( node->Parent )
      {
      return node;
      }
    }
}

template <class TGraph>
void
BoykovMinCutGraphFilter<TGraph>
::Augment( EdgePointerType middle )
{
  NodePointerType node;
  EdgePointerType edge;
  WeightType bottleneck;
  WeightType weight;

  /* 1. find the bottleneck capacity **/
  /* the source tree **/
  bottleneck = this->m_Output->GetEdgeWeight( middle );
  for( node = this->m_Output->GetSourceNodePointer( middle ); ;
    node = this->m_Output->GetTargetNodePointer( edge ) )
    {
    edge = node->Parent;
    if( edge == this->m_TerminalEdge )
      {
      break;
      }
    weight = this->m_Output->GetEdgeWeight(
      this->m_Output->GetReverseEdgePointer( edge ) );
    if( bottleneck > weight )
      {
      bottleneck = weight;
      }
    }
  if( bottleneck > this->m_Output->GetNodeWeight( node ) )
    {
    bottleneck = this->m_Output->GetNodeWeight( node );
    }

  /* the sink tree **/
  for( node = this->m_Output->GetTargetNodePointer( middle ); ;
    node = this->m_Output->GetTargetNodePointer( edge ) )
    {
    edge = node->Parent;
    if( edge == this->m_TerminalEdge )
      {
      break;
      }
    if( bottleneck > this->m_Output->GetEdgeWeight( edge ) )
      {
      bottleneck = this->m_Output->GetEdgeWeight(edge);
      }
    }
  if( bottleneck > -this->m_Output->GetNodeWeight( node ) )
    {
    bottleneck = -this->m_Output->GetNodeWeight( node );
    }

  /* 2. Augmenting **/
  /* the source tree **/
  this->m_Output->AddEdgeWeight(
    this->m_Output->GetReverseEdgePointer( middle ), bottleneck );
  this->m_Output->AddEdgeWeight( middle, -bottleneck );
  for( node = this->m_Output->GetNodePointer( middle->SourceIdentifier ); ;
    node = this->m_Output->GetNodePointer( edge->TargetIdentifier ) )
    {
    edge = node->Parent;
    if( edge == this->m_TerminalEdge )
      {
      break;
      }
    this->m_Output->AddEdgeWeight( edge, bottleneck );
    this->m_Output->AddEdgeWeight(
      this->m_Output->GetReverseEdgePointer( edge ), -bottleneck );
    if( this->IsEqual( this->m_Output->GetEdgeWeight(
      this->m_Output->GetReverseEdgePointer( edge ) ), this->m_WeightZero ) )
      {
      /* add node to the adoption list */
      node->Parent = this->m_OrphanEdge;
      this->m_Orphans.push_back( node );
      }
    }
  this->m_Output->AddNodeWeight( node, -bottleneck );
  if( this->IsEqual( this->m_Output->GetNodeWeight( node ),
    this->m_WeightZero ) )
    {
    /* add node to the adoption list */
    node->Parent = this->m_OrphanEdge;
    this->m_Orphans.push_back( node );
    }

  /* the sink tree **/
  for( node = this->m_Output->GetTargetNodePointer( middle ); ;
    node = this->m_Output->GetTargetNodePointer( edge ) )
    {
    edge = node->Parent;
    if( edge == this->m_TerminalEdge )
      {
      break;
      }
    this->m_Output->AddEdgeWeight(
      this->m_Output->GetReverseEdgePointer( edge->Identifier ), bottleneck );
    this->m_Output->AddEdgeWeight( edge, -bottleneck );
    if( this->IsEqual( this->m_Output->GetEdgeWeight( edge ),
      this->m_WeightZero ) )
      {
      /* add node to the adoption list */
      node->Parent = this->m_OrphanEdge;
      this->m_Orphans.push_back( node );
      }
    }
  this->m_Output->AddNodeWeight( node, bottleneck );
  if( this->IsEqual(
    this->m_Output->GetNodeWeight( node ), this->m_WeightZero ) )
    {
    /* add node to the adoption list */
    node->Parent = this->m_OrphanEdge;
    this->m_Orphans.push_back( node );
    }
  this->m_MaxFlow += bottleneck;
}

template <class TGraph>
void
BoykovMinCutGraphFilter<TGraph>
::ProcessSourceOrphan( NodePointerType orphan )
{
  NodePointerType node;
  EdgePointerType edge_min = NULL;
  EdgePointerType edge;
  int distance;
  int distance_min = NumericTraits<int>::max();
  WeightType weight;

  /* trying to find a new Parent */
  typename EdgeIdentifierContainerType::iterator it;
  for( it = orphan->OutgoingEdges.begin();
    it != orphan->OutgoingEdges.end(); ++it)
    {
    weight = this->m_Output->GetEdgeWeight(
      this->m_Output->GetReverseEdgePointer( *it ) );
    if( !this->IsEqual( weight, this->m_WeightZero ) )
      {
      node = this->m_Output->GetTargetNodePointer(*it);
      edge = node->Parent;
      if( !node->IsSink && edge != NULL )
        {
        /* checking the origin of node **/
        distance = 0;
        while( true )
          {
          if( node->TimeStamp == this->m_GlobalTime )
            {
            distance += node->DistanceToTerminal;
            break;
            }
          edge = node->Parent;
          distance++;
          if( edge == this->m_TerminalEdge )
            {
            node->TimeStamp = this->m_GlobalTime;
            node->DistanceToTerminal = 1;
            break;
            }
          if( edge == this->m_OrphanEdge )
            {
            distance = NumericTraits<int>::max();
            break;
            }
          node = this->m_Output->GetTargetNodePointer( edge );
          }

        /* node originates from the source - done **/
        if( distance < NumericTraits<int>::max() )
          {
          if( distance < distance_min )
            {
            edge_min = this->m_Output->GetEdgePointer( *it );
            distance_min = distance;
            }
          /* set marks along the path */
          for( node = this->m_Output->GetTargetNodePointer( *it );
            node->TimeStamp != this->m_GlobalTime;
            node = this->m_Output->GetTargetNodePointer( node->Parent ) )
            {
            node->TimeStamp = this->m_GlobalTime;
            node->DistanceToTerminal = distance--;
            }
          }
        }
      }
    }

  orphan->Parent = edge_min;
  if( orphan->Parent )
    {
    orphan->TimeStamp = this->m_GlobalTime;
    orphan->DistanceToTerminal = distance_min + 1;
    }
  else
    {
    /* no parent is found */
    orphan->TimeStamp = 0;

    /* process neighbors */
    for( it = orphan->OutgoingEdges.begin();
      it != orphan->OutgoingEdges.end(); ++it )
      {
      node = this->m_Output->GetTargetNodePointer( *it );
      edge = node->Parent;
      if(!node->IsSink && edge != NULL)
        {
        weight = this->m_Output->GetEdgeWeight(
          this->m_Output->GetReverseEdgePointer( *it ) );
        if( !IsEqual( weight, this->m_WeightZero ) )
          {
          this->SetActiveNode(node);
          }
        if( edge != this->m_TerminalEdge && edge != this->m_OrphanEdge
           && this->m_Output->GetTargetNodePointer( edge ) == orphan )
          {
          /* add node to the adoption list */
          node->Parent = this->m_OrphanEdge;
          this->m_Orphans.push_back(node);
          }
        }
      }
    }
}

template <class TGraph>
void
BoykovMinCutGraphFilter<TGraph>
::ProcessSinkOrphan( NodePointerType orphan )
{
  NodePointerType node;
  EdgePointerType edge_min = NULL, edge;
  int distance;
  int distance_min = NumericTraits<int>::max();
  WeightType weight;

  /* trying to find a new Parent */
  typename EdgeIdentifierContainerType::iterator it;
  for( it = orphan->OutgoingEdges.begin();
    it != orphan->OutgoingEdges.end(); ++it )
    {
    weight = this->m_Output->GetEdgeWeight( *it );
    if( !this->IsEqual( weight, this->m_WeightZero ) )
      {
      node = this->m_Output->GetTargetNodePointer( *it );
      if( node->IsSink && ( edge = node->Parent ) )
        {
        /* checking the origin of node **/
        distance = 0;
        while( true )
          {
          if( node->TimeStamp == this->m_GlobalTime )
            {
            distance += node->DistanceToTerminal;
            break;
            }
          edge = node->Parent;
          distance++;
          if( edge == this->m_TerminalEdge )
            {
            node->TimeStamp = this->m_GlobalTime;
            node->DistanceToTerminal = 1;
            break;
            }
          if( edge == this->m_OrphanEdge )
            {
            distance = NumericTraits<int>::max();
            break;
            }
          node = this->m_Output->GetTargetNodePointer(edge);
          }

        /* node originates from the sink - done **/
        if( distance < NumericTraits<int>::max() )
          {
          if( distance < distance_min )
            {
            edge_min = this->m_Output->GetEdgePointer(*it);
            distance_min = distance;
            }
          /* set marks along the path */
          for( node = this->m_Output->GetTargetNodePointer( *it );
            node->TimeStamp != this->m_GlobalTime;
            node = this->m_Output->GetTargetNodePointer( node->Parent ) )
            {
            node->TimeStamp = this->m_GlobalTime;
            node->DistanceToTerminal = distance--;
            }
          }
        }
      }
    }

  orphan->Parent = edge_min;
  if( orphan->Parent )
    {
    orphan->TimeStamp = this->m_GlobalTime;
    orphan->DistanceToTerminal = distance_min + 1;
    }
  else
    {
    /* no parent is found */
    orphan->TimeStamp = 0;

    /* process neighbors */
    for(it = orphan->OutgoingEdges.begin();
      it != orphan->OutgoingEdges.end(); ++it)
      {
      node = this->m_Output->GetTargetNodePointer( *it );
      edge = node->Parent;
      if( node->IsSink && edge != NULL )
        {
        weight = this->m_Output->GetEdgeWeight( *it );
        if( !this->IsEqual( weight, this->m_WeightZero ) )
          {
          this->SetActiveNode( node );
          }
        if( edge != this->m_TerminalEdge && edge != this->m_OrphanEdge
           && this->m_Output->GetTargetNodePointer( edge ) == orphan )
          {
          /* add node to the adoption list */
          node->Parent = this->m_OrphanEdge;
          this->m_Orphans.push_back( node );
          }
        }
      }
    }
}

template <class TGraph>
void
BoykovMinCutGraphFilter<TGraph>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

} // end namespace itk

#endif



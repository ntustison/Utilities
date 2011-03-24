/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBoykovMinCutGraphFilter.h,v $
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
#ifndef __itkBoykovMinCutGraphFilter_h_
#define __itkBoykovMinCutGraphFilter_h_

#include "itkInPlaceGraphFilter.h"
#include "vnl/vnl_math.h"

#include <deque>

namespace itk
{

/** \class itkBoykovMinCutGraphFilter
 * \brief Partitions a graph into two subsets based on the min-cut criterion.
 *
 * \par
 * This class implements the min-cut/max-flow algorithm of Boykov and
 * Kolmogorov given in the reference below.  Given a weighted, directed
 * graph with a source node and a sink node, finding the max-flow
 * through the specified network is equivalent to finding the cut which
 * minimizes the sum of the weights between the two subsets associated
 * with the sink and the source.
 *
 * \par
 * Other max-flow/min-cut algorithms exist such as Goldberg-Tarjan style
 * "push-relabel" methods and Ford-Fulkerson "augmenting paths" methods.
 * Even though these other algorithms have lower theoretical complexity
 * (e.g. O(|E||V|^2) for the Dinic algorithm versus O(|E||V|^2|C|) for
 * the Boykov-Kolmogorov method where |C| is the cost of the minimum cut),
 * the Boykov-Kolmogorov method has been tailored to typical problems in
 * vision making its experimental performance better than other
 * algorithms.
 *
 * \par
 * This class is derived from the InPlaceGraphFilter since the output
 * is simply a labeling of the input (each node is designated to be
 * associated with either the source or the sink terminal.  Also,
 * instead of actually creating the terminal nodes and associated
 * links, memory usage is minimized by maintaing the difference between
 * the terminal link weights as the node weight.
 *
 * \par REFERENCE
 * Y. Boykov and V. Kolmogorov, "An Experimental Comparison of Min-
 * Cut/Max-Flow Algorithms for Energy Minimization in Vision,"
 * IEEE-PAMI, 26(9):1124-1137, 2004.
 *
 **/

template<class TGraph>
class ITK_EXPORT BoykovMinCutGraphFilter
: public InPlaceGraphFilter<TGraph>
{
public:
  /** Standard "Self" typedef. */
  typedef BoykovMinCutGraphFilter Self;
  typedef InPlaceGraphFilter<TGraph> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( BoykovMinCutGraphFilter, GraphToGraphFilter );

  /** Hold on to the type information specified by the template parameters. */
  typedef TGraph GraphType;
  typedef typename GraphType::GraphTraitsType       GraphTraitsType;
  typedef typename GraphType::NodeIterator          NodeIteratorType;
  typedef typename GraphType::EdgeIterator          EdgeIteratorType;
  typedef typename GraphTraitsType::NodeType        NodeType;
  typedef typename GraphTraitsType::EdgeType        EdgeType;
  typedef typename GraphTraitsType::NodePointerType NodePointerType;
  typedef typename GraphTraitsType::EdgePointerType EdgePointerType;
  typedef typename GraphTraitsType::NodeWeightType  WeightType;
  typedef typename GraphTraitsType
     ::EdgeIdentifierContainerType                  EdgeIdentifierContainerType;

  /** Define other types */
  typedef double                                    RealType;
  typedef std::deque<NodePointerType>               NodeListType;

  itkGetMacro( MaxFlow, WeightType );

protected:
  BoykovMinCutGraphFilter();
  ~BoykovMinCutGraphFilter() {}
  void PrintSelf( std::ostream& os, Indent indent ) const;

  void GenerateData();
  virtual void GenerateOutputInformation(){}; // do nothing

private:
  BoykovMinCutGraphFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  /** Private functions for processing the graph */
  void GenerateMinCut( void );
  void Initialize( void );
  void SetActiveNode( NodePointerType );
  NodePointerType GetNextActiveNode( void );
  void Augment( EdgePointerType );
  void ProcessSourceOrphan( NodePointerType );
  void ProcessSinkOrphan( NodePointerType );

  inline bool IsEqual( WeightType m, WeightType n )
    { return ( vnl_math_abs( static_cast<RealType>( m-n ) ) < 1e-10 ); }

  WeightType m_WeightZero;

  /** Dynamic lists of active nodes and orphans */
  NodeListType m_ActiveNodes;
  NodeListType m_Orphans;

  /** Two constant edges used for finding the min-cut/max-flow */
  EdgePointerType m_TerminalEdge;
  EdgePointerType m_OrphanEdge;

  /** Other variables */
  int m_GlobalTime;
  WeightType m_MaxFlow;

  typename GraphType::Pointer m_Output;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBoykovMinCutGraphFilter.txx"
#endif

#endif

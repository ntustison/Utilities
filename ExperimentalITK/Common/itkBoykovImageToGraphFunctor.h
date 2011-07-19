/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBoykovImageToGraphFunctor.h,v $
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
#ifndef __itkBoykovImageToGraphTraits_h
#define __itkBoykovImageToGraphTraits_h

#include "itkDefaultImageToGraphFunctor.h"

namespace itk
{

/** \class itkBoykovImageToGraphFunctor
 * \brief Class which defines node/edge weighting in constructing a
 *        graph from an image.
 *
 * \par
 * Using the ImageToGraphFilter class, one wishes to construct a
 * graph from a given input image.  This functor determines whether
 * or not a given pixel constitutes a node, what weight that node
 * should have, and what edge weight should be assigned an edge
 * between two nodes.  The weighting scheme is based on the
 * reference below.  Also note that the weights of the nodes are
 * used to n-link affinities whereas the edge weights store the
 * t-link affinities.  This eliminates unnecessary memory usage.
 *
 * \par INPUTS
 * This class should be used in conjunction with the
 * itkBoykovFilter class.  As such, the resulting graph is
 * going to have a binary labeling (sink vs. source).  The input
 * consists partly of 2 probability (likelihood) images of the
 * same dimension and size as the input image.  These can be
 * derived from non-parametric techniques such as Parzen windowing,
 * etc.  Optional: 2 IndexContainerTypes of pixel indices can be
 * specified to "hard-constrain" those pixels to be of a specific
 * labeling.  Other parameters include lambda and sigma which are
 * described below.
 *
 * \par REFERENCE
 * Y. Boykov and M.-P. Jolly, "Interactive Graph Cuts for Optimal Boundary
 * & Region Segmentation of Objects in N-D Images," ICCV, 2001, 105-112.
 *
 **/

template<typename TInputImage, typename TOutputGraph>
class BoykovImageToGraphFunctor
: public DefaultImageToGraphFunctor<TInputImage, TOutputGraph>
{
public:
  /** Standard class typedefs. */
  typedef BoykovImageToGraphFunctor Self;
  typedef DefaultImageToGraphFunctor<TInputImage, TOutputGraph> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  typedef TInputImage InputImageType;
  typedef TOutputGraph OutputGraphType;

  typedef typename Superclass::IndexType            IndexType;
  typedef typename Superclass::PixelType            PixelType;
  typedef typename Superclass::NodeType             NodeType;
  typedef typename Superclass::EdgeType             EdgeType;
  typedef typename Superclass::NodeIteratorType     NodeIteratorType;
  typedef typename Superclass::EdgeIteratorType     EdgeIteratorType;
  typedef typename Superclass::NodePointerType      NodePointerType;
  typedef typename Superclass::EdgePointerType      EdgePointerType;
  typedef typename Superclass::NodeWeightType       NodeWeightType;
  typedef typename Superclass::EdgeWeightType       EdgeWeightType;
  typedef typename Superclass::NodeImageType        NodeImageType;
  typedef typename Superclass::EdgeIdentifierContainerType
                                                    EdgeIdentifierContainerType;

  typedef double RealType;
  typedef std::vector<IndexType> IndexContainerType;

  /** define virtual functions */
  virtual EdgeWeightType GetEdgeWeight( IndexType, IndexType );
  virtual NodeWeightType GetNodeWeight( IndexType idx )
    { return this->GetSinkDataTerm( idx ) - this->GetSourceDataTerm( idx ); }
  virtual void NormalizeGraph( NodeImageType *, OutputGraphType * );

  NodeWeightType GetSourceDataTerm( IndexType );
  NodeWeightType GetSinkDataTerm( IndexType );
  EdgeWeightType GetSmoothnessTerm( IndexType idx1, IndexType idx2 )
    { return this->GetEdgeWeight( idx1, idx2 ); }

  /** lambda - factor which specifies the relative weighting
    * between the regional properties and the boundary properties.
    */
  itkGetMacro( Lambda, RealType );
  itkSetMacro( Lambda, RealType );

  /** sigma - standard deviation associated with the weighting
    * of the distances between neighboring pixels/
    */
  itkGetMacro( Sigma, RealType );
  itkSetMacro( Sigma, RealType );

  /** Declare the probability image type */
  typedef Image<RealType, InputImageType::ImageDimension> LikelihoodImageType;

  void SetSourceLikelihoodImage( const LikelihoodImageType * );
  void SetSinkLikelihoodImage( const LikelihoodImageType * );

  void SetSourceIndexContainer( const IndexContainerType source )
    { this->m_SourceIndexContainer = source; }
  void SetSinkIndexContainer( const IndexContainerType sink )
    { this->m_SinkIndexContainer = sink; }

protected:
  BoykovImageToGraphFunctor();
  ~BoykovImageToGraphFunctor() {}
  void PrintSelf( std::ostream& os, Indent indent ) const;

private:
  BoykovImageToGraphFunctor( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  RealType m_Lambda;
  RealType m_Sigma;

  IndexContainerType m_SourceIndexContainer;
  IndexContainerType m_SinkIndexContainer;
};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBoykovImageToGraphFunctor.hxx"
#endif

#endif

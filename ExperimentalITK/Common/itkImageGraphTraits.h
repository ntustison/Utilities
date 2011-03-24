/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImageGraphTraits.h,v $
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
#ifndef __itkImageGraphTraits_h
#define __itkImageGraphTraits_h

#include "itkDefaultGraphTraits.h"
#include "itkIndex.h"

namespace itk
{

/** \class ImageGraphTraits
 *  \brief Base class that associates a graph node with an image pixel.
 *
 *  Many graph algorithms require the representation of an image pixel
 *  by a node in the graph.  This is the base class for graph traits
 *  classes that are used in such algorithms.  Each node structure
 *  has an associate IndexType which contains the pixel index of the
 *  pixel the node represents.
 */

template <typename TWeight, unsigned int VImageDimension>
class ImageGraphTraits : public DefaultGraphTraits<TWeight, TWeight>
{
public:
  typedef ImageGraphTraits Self;
  typedef DefaultGraphTraits<TWeight, TWeight> Superclass;

  typedef Index<VImageDimension> IndexType;
  typedef TWeight NodeWeightType;
  typedef TWeight EdgeWeightType;
  typedef typename Superclass::NodeIdentifierType NodeIdentifierType;
  typedef typename Superclass::EdgeIdentifierType EdgeIdentifierType;
  typedef typename Superclass::EdgeType           EdgeType;
  typedef typename Superclass::EdgePointerType    EdgePointerType;

  typedef typename Superclass::EdgeIdentifierContainerType
                                               EdgeIdentifierContainerType;

  struct NodeType;
  typedef NodeType* NodePointerType;

  struct NodeType
    {
    NodeIdentifierType Identifier;
    EdgeIdentifierContainerType IncomingEdges;
    EdgeIdentifierContainerType OutgoingEdges;
    NodeWeightType Weight;
    IndexType ImageIndex;
    };
};

} // end namespace itk

#endif


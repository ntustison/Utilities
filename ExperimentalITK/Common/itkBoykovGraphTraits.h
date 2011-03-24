/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBoykovGraphTraits.h,v $
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
#ifndef __itkBoykovGraphTraits_h
#define __itkBoykovGraphTraits_h

#include "itkImageGraphTraits.h"

namespace itk
{

/**
 * Graph traits class for use with the BoykovMinCutFilter class.
 */

template <typename TWeight = short, unsigned int VImageDimension = 3>
class BoykovGraphTraits : public ImageGraphTraits<TWeight, VImageDimension>
{
public:
  typedef BoykovGraphTraits Self;
  typedef ImageGraphTraits<TWeight, VImageDimension> Superclass;

  typedef TWeight NodeWeightType;
  typedef TWeight EdgeWeightType;
  typedef typename Superclass::NodeIdentifierType NodeIdentifierType;
  typedef typename Superclass::EdgeIdentifierType EdgeIdentifierType;
  typedef typename Superclass::EdgeIdentifierContainerType
                                                  EdgeIdentifierContainerType;
  typedef typename Superclass::IndexType          IndexType;
  typedef typename Superclass::EdgeType           EdgeType;
  typedef typename Superclass::EdgePointerType    EdgePointerType;

  struct NodeType;
  typedef NodeType* NodePointerType;

  struct NodeType
    {
    NodeIdentifierType Identifier;
    EdgeIdentifierContainerType IncomingEdges;
    EdgeIdentifierContainerType OutgoingEdges;
    EdgePointerType Parent;
    int TimeStamp;
    int DistanceToTerminal;
    bool IsSink;
    bool IsActive;
    NodeWeightType Weight;
    IndexType ImageIndex;
    };
};

} // end namespace itk

#endif


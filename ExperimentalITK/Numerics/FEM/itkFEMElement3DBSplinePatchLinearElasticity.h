/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMElement3DBSplinePatchLinearElasticity.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:22:48 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkFEMElement3DBSplinePatchLinearElasticity_h
#define __itkFEMElement3DBSplinePatchLinearElasticity_h

#include "itkFEMElement3DBSplinePatch.h"
#include "itkFEMElement3DLinearElasticity.h"

namespace itk {
namespace fem {




/**
 * \class Element3DBSplinePatchLinearElasticity
 * \brief n-noded finite element class in 3D space for linear elasticity problem
 */
class Element3DBSplinePatchLinearElasticity : public Element3DLinearElasticity<Element3DBSplinePatch>
{
FEM_CLASS(Element3DBSplinePatchLinearElasticity,Element3DLinearElasticity<Element3DBSplinePatch>)
public:

  HANDLE_ELEMENT_LOADS();

  /**
   * Default constructor only clears the internal storage
   */
  Element3DBSplinePatchLinearElasticity();

  /**
   * Construct an element by specifying pointers to
   * an array of n points and a material.
   */
  Element3DBSplinePatchLinearElasticity(
      NodeIDType ns_[], 
      Material::ConstPointer p_ );

}; // class Element3DBSplinePatchLinearElasticity

FEM_CLASS_INIT(Element3DBSplinePatchLinearElasticity)




}} // end namespace itk::fem

#endif  // #ifndef __itkFEMElement3DBSplinePatchLinearElasticity_h

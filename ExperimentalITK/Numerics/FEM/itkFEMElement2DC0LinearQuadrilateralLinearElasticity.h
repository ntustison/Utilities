/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMElement2DC0LinearQuadrilateralLinearElasticity.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:22:47 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkFEMElement2DC0LinearQuadrilateralLinearElasticity_h
#define __itkFEMElement2DC0LinearQuadrilateralLinearElasticity_h

#include "itkFEMElement2DC0LinearQuadrilateral.h"
#include "itkFEMElement2DLinearElasticity.h"

namespace itk {
namespace fem {




/**
 * \class Element2DC0LinearQuadrilateralLinearElasticity
 * \brief 4-noded finite element class in 2D space for linear elasticity problem
 */
class Element2DC0LinearQuadrilateralLinearElasticity : public Element2DLinearElasticity<Element2DC0LinearQuadrilateral>
{
FEM_CLASS(Element2DC0LinearQuadrilateralLinearElasticity,Element2DLinearElasticity<Element2DC0LinearQuadrilateral>)
public:

  HANDLE_ELEMENT_LOADS();

  /**
   * Default constructor only clears the internal storage
   */
  Element2DC0LinearQuadrilateralLinearElasticity();

  /**
   * Construct an element by specifying pointers to
   * 4 points and a material.
   */
  Element2DC0LinearQuadrilateralLinearElasticity(
      NodeIDType n1_, 
      NodeIDType n2_,
      NodeIDType n3_,
      NodeIDType n4_,
      Material::ConstPointer p_ );

}; // class Element2DC0LinearQuadrilateralLinearElasticity

FEM_CLASS_INIT(Element2DC0LinearQuadrilateralLinearElasticity)




}} // end namespace itk::fem

#endif  // #ifndef __itkFEMElement2DC0LinearQuadrilateralLinearElasticity_h

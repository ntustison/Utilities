/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMElement2DBSplinePatchMembrane.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:22:47 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkFEMElement2DBSplinePatchMembrane_h
#define __itkFEMElement2DBSplinePatchMembrane_h

#include "itkFEMElement2DBSplinePatch.h"
#include "itkFEMElement2DMembrane.h"

namespace itk {
namespace fem {




/**
 * \class Element2DBSplinePatchMembrane
 * \brief n-noded finite element class in 2D space for linear elasticity problem
 */
class Element2DBSplinePatchMembrane : public Element2DMembrane<Element2DBSplinePatch>
{
FEM_CLASS(Element2DBSplinePatchMembrane, Element2DMembrane<Element2DBSplinePatch>)
public:

  HANDLE_ELEMENT_LOADS();

  /**
   * Default constructor only clears the internal storage
   */
  Element2DBSplinePatchMembrane();

  /**
   * Construct an element by specifying pointers to
   * n points and a material.
   */
  Element2DBSplinePatchMembrane(
      NodeIDType ns_[],
      Material::ConstPointer p_ );

}; // class Element2DBSplinePatchMembrane

FEM_CLASS_INIT(Element2DBSplinePatchMembrane)




}} // end namespace itk::fem

#endif  // #ifndef __itkFEMElement2DBSplinePatchMembrane_h

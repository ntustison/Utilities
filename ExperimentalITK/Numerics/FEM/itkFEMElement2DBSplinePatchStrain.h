/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMElement2DBSplinePatchStrain.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:22:47 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkFEMElement2DBSplinePatchStrain_h
#define __itkFEMElement2DBSplinePatchStrain_h

#include "itkFEMElement2DBSplinePatch.h"
#include "itkFEMElement2DStrain.h"

namespace itk {
namespace fem {




/**
 * \class Element2DBSplinePatchStrain
 * \brief n-noded finite element class in 2D space for linear elasticity problem
 */
class Element2DBSplinePatchStrain : public Element2DStrain<Element2DBSplinePatch>
{
FEM_CLASS(Element2DBSplinePatchStrain,Element2DStrain<Element2DBSplinePatch>)
public:

  HANDLE_ELEMENT_LOADS();

  /**
   * Default constructor only clears the internal storage
   */
  Element2DBSplinePatchStrain();

  /**
   * Construct an element by specifying pointers to
   * n points and a material.
   */
  Element2DBSplinePatchStrain(
      NodeIDType ns_[],
      Material::ConstPointer m_);

}; // class Element2DBSplinePatchStrain

FEM_CLASS_INIT(Element2DBSplinePatchStrain)




}} // end namespace itk::fem

#endif  // #ifndef __itkFEMElement2DBSplinePatchStrain_h

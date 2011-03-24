/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMElement3DBSplinePatchMembrane.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:22:48 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkFEMElement3DBSplinePatchMembrane_h
#define __itkFEMElement3DBSplinePatchMembrane_h

#include "itkFEMElement3DBSplinePatch.h"
#include "itkFEMElement3DMembrane.h"

namespace itk {
namespace fem {




/**
 * \class Element3DBSplinePatchMembrane
 * \brief n-noded finite element class in 3D space for linear elasticity problem
 */
class Element3DBSplinePatchMembrane : public Element3DMembrane<Element3DBSplinePatch>
{
FEM_CLASS(Element3DBSplinePatchMembrane,Element3DMembrane<Element3DBSplinePatch>)
public:

  HANDLE_ELEMENT_LOADS();

  /**
   * Default constructor only clears the internal storage
   */
  Element3DBSplinePatchMembrane();

  /**
   * Construct an element by specifying pointers to
   * an array of n points and a material.
   */
  Element3DBSplinePatchMembrane(
      NodeIDType ns_[], 
      Material::ConstPointer p_ );

}; // class Element3DBSplinePatchMembrane

FEM_CLASS_INIT(Element3DBSplinePatchMembrane)




}} // end namespace itk::fem

#endif  // #ifndef __itkFEMElement3DBSplinePatchMembrane_h

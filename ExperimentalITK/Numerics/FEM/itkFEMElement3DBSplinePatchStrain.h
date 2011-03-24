/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMElement3DBSplinePatchStrain.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:22:48 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkFEMElement3DBSplinePatchStrain_h
#define __itkFEMElement3DBSplinePatchStrain_h

#include "itkFEMElement3DBSplinePatch.h"
#include "itkFEMElement3DStrain.h"

namespace itk {
namespace fem {




/**
 * \class Element3DBSplinePatchStrain
 * \brief n-noded finite element class in 3D space for linear elasticity problem
 */
class Element3DBSplinePatchStrain : public Element3DStrain<Element3DBSplinePatch>
{
FEM_CLASS(Element3DBSplinePatchStrain,Element3DStrain<Element3DBSplinePatch>)
public:

  HANDLE_ELEMENT_LOADS();

  /**
   * Default constructor only clears the internal storage
   */
  Element3DBSplinePatchStrain();

  /**
   * Construct an element by specifying pointers to
   * an array of n points and a material.
   */
  Element3DBSplinePatchStrain(
      NodeIDType ns_[], 
      Material::ConstPointer p_ );

}; // class Element3DBSplinePatchStrain

FEM_CLASS_INIT(Element3DBSplinePatchStrain)




}} // end namespace itk::fem

#endif  // #ifndef __itkFEMElement3DBSplinePatchStrain_h

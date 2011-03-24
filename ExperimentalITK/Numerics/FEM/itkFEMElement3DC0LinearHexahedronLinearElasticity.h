/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMElement3DC0LinearHexahedronLinearElasticity.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:22:48 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkFEMElement3DC0LinearHexahedronLinearElasticity_h
#define __itkFEMElement3DC0LinearHexahedronLinearElasticity_h

#include "itkFEMElement3DC0LinearHexahedron.h"
#include "itkFEMElement3DLinearElasticity.h"

namespace itk {
namespace fem {
 



/**
 * \class Element3DC0LinearHexahedronLinearElasticity
 * \brief 8-noded finite element class in 3D space for linear elasticity problem
 */
class Element3DC0LinearHexahedronLinearElasticity : public Element3DLinearElasticity<Element3DC0LinearHexahedron>
{
FEM_CLASS(Element3DC0LinearHexahedronLinearElasticity,Element3DLinearElasticity<Element3DC0LinearHexahedron>)
public:

  HANDLE_ELEMENT_LOADS();

  /**
   * Default constructor only clears the internal storage
   */
  Element3DC0LinearHexahedronLinearElasticity();

  /**
   * Construct an element by specifying pointers to
   * an array of 8 points and a material.
   */
  Element3DC0LinearHexahedronLinearElasticity(
      NodeIDType ns_[], 
      Material::ConstPointer p_ );

}; // class Element3DC0LinearHexahedronLinearElasticity

FEM_CLASS_INIT(Element3DC0LinearHexahedronLinearElasticity)




}} // end namespace itk::fem

#endif  // #ifndef __itkFEMElement3DC0LinearHexahedronLinearElasticity_h

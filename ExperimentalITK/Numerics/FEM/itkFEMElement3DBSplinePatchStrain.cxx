/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMElement3DBSplinePatchStrain.cxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:22:48 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

// disable debug warnings in MS compiler
#ifdef _MSC_VER
#pragma warning(disable: 4786)
#endif

#include "itkFEMElement3DBSplinePatchStrain.h"

namespace itk {
namespace fem {


Element3DBSplinePatchStrain
::Element3DBSplinePatchStrain() : Superclass()
{
}

Element3DBSplinePatchStrain
::Element3DBSplinePatchStrain(
      NodeIDType ns_[],
      Material::ConstPointer m_) : Superclass()
{
  // Set the geometrical points
  for (unsigned int k = 0; k < this->GetNumberOfNodes(); k++)
    this->SetNode( k, ns_[k] );

  /*
   * Initialize the pointer to material object and check that
   * we were given the pointer to the right class.
   * If the material class was incorrect an exception is thrown.
   */
  if( (m_mat=dynamic_cast<const MaterialLinearElasticity*>(&*m_)) == 0 )
  {
    throw FEMExceptionWrongClass(__FILE__,__LINE__,"Element3DBSplinePatchStrain::Element3DBSplinePatchStrain()");
  }
}




FEM_CLASS_REGISTER(Element3DBSplinePatchStrain)




}} // end namespace itk::fem

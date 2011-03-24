/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMElements.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:22:49 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/


/**
 * \file itkFEMElements.h
 * \brief Include all finite element classes defined in FEM toolkit.
 *
 * To make sure you have everything, just include this header file.
 */
#include "itkFEMElement2DC0LinearLineStress.h"
#include "itkFEMElement2DC1Beam.h"
#include "itkFEMElement2DC0LinearTriangularMembrane.h"
#include "itkFEMElement2DC0LinearTriangularStress.h"
#include "itkFEMElement2DC0LinearTriangularStrain.h"
#include "itkFEMElement2DC0LinearQuadrilateralMembrane.h"
#include "itkFEMElement2DC0LinearQuadrilateralStress.h"
#include "itkFEMElement2DC0LinearQuadrilateralStrain.h"
#include "itkFEMElement3DC0LinearTetrahedronMembrane.h"
#include "itkFEMElement3DC0LinearTetrahedronStrain.h"
#include "itkFEMElement3DC0LinearHexahedronMembrane.h"
#include "itkFEMElement3DC0LinearHexahedronStrain.h"
#include "itkFEMElement2DC0QuadraticTriangularStress.h"
#include "itkFEMElement2DC0QuadraticTriangularStrain.h"

#include "itkFEMElement2DBSplinePatchLinearElasticity.h"
#include "itkFEMElement2DBSplinePatchMembrane.h"
#include "itkFEMElement2DBSplinePatchStress.h"
#include "itkFEMElement2DBSplinePatchStrain.h"
#include "itkFEMElement3DBSplinePatchLinearElasticity.h"
#include "itkFEMElement3DBSplinePatchMembrane.h"
#include "itkFEMElement3DBSplinePatchStrain.h"

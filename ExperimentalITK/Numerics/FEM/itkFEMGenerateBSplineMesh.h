/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMGenerateBSplineMesh.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:22:49 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkFEMGenerateBSplineMesh_h
#define __itkFEMGenerateBSplineMesh_h

#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "itkFEMSolver.h"
#include "itkFEMBSplineShapeFunctions.h"
#include "itkImage.h"
#include "itkPoint.h"

#define CLAMPED_BSPLINES

namespace itk {
namespace fem {

/**
 * \function Generate2DRectilinearBSplineMesh
 * \brief Use this function to generate 2D BSpline meshes in Solver.
 *
 * This function uses the generic quadrilateral elements
 * to build meshes that can be used with specific elements for solving 
 * membrane or linear elasticity problems.
 *
 * See other functions if you need to constuct the mesh from other types
 * of elements.
 *
 * \note All elements will be created by copying the existing element which
 *       is passed to the function. Only number and node pointers will
 *       be changed in copied element. Make sure that this element has material
 *       class and any other properties defined before generating a mesh.
 *
 * \sa Generate3DRectilinearBSplineMesh
 */

/**
 * Generate a rectangular mesh of quadrilateral elements
 */
 
void Generate2DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, vnl_vector<double>& orig, vnl_vector<double>& size, vnl_vector<double>& Nel, 
                                      unsigned int BSplineOrder);

void Generate2DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, vnl_vector<double>& orig, vnl_vector<double>& size, vnl_vector<double>& Nel, 
                                      unsigned int BSplineOrder, vnl_vector<double>& knots0, vnl_vector<double>& knots1);

void Generate2DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, vnl_vector<double>& orig, vnl_vector<double>& size, vnl_vector<double>& Nel, 
                                      unsigned int BSplineOrder, Image<char, 2>::Pointer maskImage);

void Generate2DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, vnl_vector<double>& orig, vnl_vector<double>& size, vnl_vector<double>& Nel, 
                                      unsigned int BSplineOrder, vnl_vector<double>& knots0, vnl_vector<double>& knots1, Image<char, 2>::Pointer maskImage);
         
/**
 * \fucntion Generate3DRectilinearBSplineMesh
 * \brief Use this function to generate 3D B-Spline meshes in Solver.
 *
 * Generate a rectangular mesh of hexahedral B-spline elements.
 *
 * \sa Generate3DRectilinearBSplineMesh
 */
// Assume the knot vector will be uniform 
void Generate3DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, vnl_vector<double>& orig, vnl_vector<double>& size, vnl_vector<double>& Nel, 
                                      unsigned int BSplineOrder);

void Generate3DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, vnl_vector<double>& orig, vnl_vector<double>& size, vnl_vector<double>& Nel, 
                                      unsigned int BSplineOrder, vnl_vector<double>& knots0, vnl_vector<double>& knots1, vnl_vector<double>& knots2);

void Generate3DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, vnl_vector<double>& orig, vnl_vector<double>& size, vnl_vector<double>& Nel, 
                                      unsigned int BSplineOrder, Image<char, 3>::Pointer maskImage);

void Generate3DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, vnl_vector<double>& orig, vnl_vector<double>& size, vnl_vector<double>& Nel, 
                                      unsigned int BSplineOrder, vnl_vector<double>& knots0, vnl_vector<double>& knots1, vnl_vector<double>& knots2, Image<char, 3>::Pointer maskImage);
          
/** Helper function to calculate extended B-spline shape functions */
double calculate_e(vnl_vector<unsigned int> i, vnl_vector<unsigned int> j, vnl_vector<unsigned int> l, unsigned int n);


/** Dummy functions that shouldn't be executed*/

void Generate2DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, vnl_vector<double>& orig, vnl_vector<double>& size, vnl_vector<double>& Nel, 
                                      unsigned int BSplineOrder, Image<char, 3>::Pointer maskImage);

void Generate2DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, vnl_vector<double>& orig, vnl_vector<double>& size, vnl_vector<double>& Nel, 
                                      unsigned int BSplineOrder, vnl_vector<double>& knots0, vnl_vector<double>& knots1, Image<char, 3>::Pointer maskImage);

void Generate3DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, vnl_vector<double>& orig, vnl_vector<double>& size, vnl_vector<double>& Nel, 
                                      unsigned int BSplineOrder, Image<char, 2>::Pointer maskImage);

void Generate3DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, vnl_vector<double>& orig, vnl_vector<double>& size, vnl_vector<double>& Nel, 
                                      unsigned int BSplineOrder, vnl_vector<double>& knots0, vnl_vector<double>& knots1, vnl_vector<double>& knots2, Image<char, 2>::Pointer maskImage);

          
}} // end namespace itk::fem

#endif // #ifndef __itkFEMGenerateBSplineMesh_h

/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMBSplineShapeFunctions.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:22:47 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkFEMBSplineShapeFunctions_h
#define __itkFEMBSplineShapeFunctions_h

#include "vnl/vnl_real_polynomial.h"

namespace itk {
namespace fem {

/**
 * \class BSplineShapeFunctions
 * \brief Generates the B-spline shape functions given an order and knot vector.
 * This class contains the routines necessary for constructing the B-spline
 * shape functions.  Antecedent to creating the shape functions, the order and
 * the knot vector must be specified.  The class assumes that the specified order 
 * and knot vector will be correct, e.g. order >= 1 and knot vector is monotonically
 * increasing.  Extensible functionality will require additional routines.
 */
class BSplineShapeFunctions
{
public:
 
   typedef double                 RealType;
   typedef vnl_vector<RealType>   VectorType;
   typedef vnl_real_polynomial    PolynomialType;

   /**
    * Empty constructor and destructor
    */   
   BSplineShapeFunctions() {};
  ~BSplineShapeFunctions() {}; 

   /**
    * Input:  Order -
    *         whichBasisFunction - 
    *         whichPiece - Each basis function is a piecewise polynomial.  This parameter specifies which piece.
    *         knots -
    * Output: polynomial piece of the specified B-spline basis function
    * Note:  This function calls the recursive version of itself.
    */
   PolynomialType GenerateBSplineShapeFunction(unsigned short order, unsigned int whichBasisFunction, unsigned int whichPiece, VectorType knots);
   /**
    * Input:  Order -
    *         whichBasisFunction - 
    *         whichPiece - Each basis function is a piecewise polynomial.  This parameter specifies which piece.
    * Output: polynomial piece of the specified B-spline basis function
    * Note:  Recursive function based on the Cox-DeBoor recurrence relation.
    */
   PolynomialType GenerateBSplineShapeFunction(unsigned short order, unsigned int whichBasisFunction, unsigned int whichPiece);
   
   /**
    * Input:  value - Parametric value at which to evaluate the specified basis function
    *         Order -
    *         whichBasisFunction - 
    *         knots -
    * Output: polynomial piece of the specified B-spline basis function
    */
   RealType CalculateBSplineShapeFunctionAtValue(RealType value, unsigned short order, unsigned int whichBasisFunction, VectorType knots);
   RealType CalculateBSplineShapeFunctionAtValue(RealType value, unsigned short order, unsigned int whichBasisFunction);

   void SetOrder(unsigned short order) { m_order = order; };
   void SetKnots(VectorType knots) { m_knots = knots; };   
   
private:  

   unsigned short m_order;   
   VectorType m_knots;
};   

}} // end namespace itk::fem

#endif  // #ifndef __itkFEMBSplineShapeFunctions_h

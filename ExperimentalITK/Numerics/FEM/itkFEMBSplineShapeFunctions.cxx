/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMBSplineShapeFunctions.cxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:22:47 $
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

#include "itkFEMBSplineShapeFunctions.h"

namespace itk {
namespace fem {

BSplineShapeFunctions::PolynomialType 
BSplineShapeFunctions::
GenerateBSplineShapeFunction(unsigned short order, unsigned int whichBasisFunction, unsigned int whichPiece, VectorType knots)
{
   m_knots = knots;
   return this->GenerateBSplineShapeFunction(order, whichBasisFunction, whichPiece);
}

BSplineShapeFunctions::PolynomialType 
BSplineShapeFunctions::
GenerateBSplineShapeFunction(unsigned short order, unsigned int whichBasisFunction, unsigned int whichPiece)
{
   VectorType tmp(2);
   PolynomialType poly1(0.0), poly2(0.0);
   RealType den;
   unsigned short p = order-1;
   unsigned short i = whichBasisFunction;   
   
   if (p == 0 && whichBasisFunction == whichPiece)
   {
      PolynomialType poly(1.0);
      return poly;
   }          
     
   // Term 1
   den = m_knots(i+p)-m_knots(i);
   if (den == 0.0) 
   {
      PolynomialType poly(0.0);
      poly1 = poly;
   }
   else
   {
      tmp(0) = 1.0;
      tmp(1) = -m_knots(i);
      tmp /= den;
      poly1 = PolynomialType(tmp) * this->GenerateBSplineShapeFunction(order-1, i, whichPiece);   
   }
   
   // Term 2
   den = m_knots(i+p+1)-m_knots(i+1);
   if (den == 0.0) 
   {
      PolynomialType poly(0.0);
      poly2 = poly;
   }
   else
   {
      tmp(0) = -1.0;
      tmp(1) = m_knots(i+p+1);
      tmp /= den;
      poly2 = PolynomialType(tmp) * this->GenerateBSplineShapeFunction(order-1, i+1, whichPiece);
   }   
   
   return (poly1 + poly2);

}
   
BSplineShapeFunctions::RealType 
BSplineShapeFunctions::
CalculateBSplineShapeFunctionAtValue(RealType value, unsigned short order, unsigned int whichBasisFunction, VectorType knots)
{
   m_knots = knots;
   return this->CalculateBSplineShapeFunctionAtValue(value, order, whichBasisFunction);
}

BSplineShapeFunctions::RealType 
BSplineShapeFunctions::
CalculateBSplineShapeFunctionAtValue(RealType value, unsigned short order, unsigned int whichBasisFunction)
{
   unsigned short p = order-1;
   unsigned short i = whichBasisFunction;  
   unsigned short m = m_knots.size();
   RealType u = value, Uright, Uleft, temp, saved; 
   VectorType N(order);   
   
   if ((i == 0 && u == m_knots[0]) || (i == m-p-1 && u == m_knots[m]))
   {
      return 1.0;
   }
   
   if (u < m_knots[i] || u >= m_knots[i+p+1])
   {
      return 0.0;
   }
   
   for (unsigned short j = 0; j <= p; j++)
      if (u >= m_knots[i+j] && u < m_knots[i+j+1])  N[j] = 1.0;
        else          N[j] = 0.0;
   for (unsigned k = 1; k <= p; k++)
   {
      if (N[0] == 0.0)  saved = 0.0;
        else   saved = ((u-m_knots[i])*N[0])/(m_knots[i+k]-m_knots[i]);
      for (unsigned short j = 0; j < p-k+1; j++)
      {
         Uleft = m_knots[i+j+1];
         Uright = m_knots[i+j+k+1];
  if (N[j+1] == 0.0)
  {
     N[j] = saved;   saved = 0.0;
  }
  else
  {
     temp = N[j+1]/(Uright-Uleft);
     N[j] = saved+(Uright-u)*temp;
     saved = (u-Uleft)*temp;
  }
      }
   }   
   return N[0];
}


}} // end namespace itk::fem

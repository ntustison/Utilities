/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMElement2DBSplinePatch.cxx,v $
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

#include "itkFEMElement2DBSplinePatch.h"
#include "vnl/vnl_math.h"
#include "vnl/algo/vnl_rpoly_roots.h"
#include "vnl/algo/vnl_rnpoly_solve.h"

namespace itk {
namespace fem {




void
Element2DBSplinePatch
::GetIntegrationPointAndWeight(unsigned int i, VectorType& pt, Float& w, unsigned int order) const
{
  // FIXME: range checking

  // default integration order
  //if (order==0) { order=DefaultIntegrationOrder; }
  //order = static_cast<unsigned int>(0.5*static_cast<float>(m_BSplineOrder)+0.5);
  order = m_BSplineOrder;

  pt.set_size(NumberOfSpatialDimensions);

  pt[0]=gaussPoint[order][i%order];
  pt[1]=gaussPoint[order][i/order];
  
  w=gaussWeight[order][i%order]*gaussWeight[order][i/order];
  
}

unsigned int
Element2DBSplinePatch
::GetNumberOfIntegrationPoints(unsigned int order) const
{
  // FIXME: range checking

  // default integration order
  //if (order==0) { order=DefaultIntegrationOrder; }
  //order = static_cast<unsigned int>(0.5*static_cast<float>(m_BSplineOrder)+0.5);
  order = m_BSplineOrder;

  return order*order;
}

void
Element2DBSplinePatch
::GenerateShapeFunctions()
{
  unsigned int m = NumberOfSpatialDimensions*m_BSplineOrder;
  vnl_vector<double> knots(m);
  
  for (unsigned int i = 0; i < m; i++)    
    knots[i] = static_cast<double>(i)/static_cast<double>(m-1);

  BSplineShapeFunctions shape;  shape.SetKnots(knots);
  
  vnl_matrix<double> shapeFunctions(NumberOfSpatialDimensions*m_BSplineOrder, m_BSplineOrder);

  for(unsigned int f = 0; f < m_BSplineOrder; f++) 
  {
    shapeFunctions.set_row(f                 , 
      (shape.GenerateBSplineShapeFunction(m_BSplineOrder, f, m_BSplineOrder-1)).coefficients());
    shapeFunctions.set_row(f+m_BSplineOrder  , 
      (shape.GenerateBSplineShapeFunction(m_BSplineOrder, f, m_BSplineOrder-1)).coefficients());   
  }
  this->SetBSplineShapeFunctions(shapeFunctions);      
}

Element2DBSplinePatch::VectorType
Element2DBSplinePatch
::ShapeFunctions( const VectorType& pt ) const
{
  Float a, b, u, v;
  a = m_localCoordinateBoundaries(0);
  b = m_localCoordinateBoundaries(1);
  u = 0.5*(b-a)*(pt[0]+1.0) + a;

  a = m_localCoordinateBoundaries(2);
  b = m_localCoordinateBoundaries(3);
  v = 0.5*(b-a)*(pt[1]+1.0) + a;

  VectorType shapeF(m_NumberOfNodes);
  
  if (!m_ExtendedBSplineStatus)   
  {
    for (unsigned int j = 0; j < m_BSplineOrder; j++)
      for (unsigned int i = 0; i < m_BSplineOrder; i++)
        shapeF[j*m_BSplineOrder+i] = PolynomialType(m_BSplineShapeFunctions.get_row(i               )).evaluate(u)*
                                 PolynomialType(m_BSplineShapeFunctions.get_row(j+m_BSplineOrder)).evaluate(v);
  }
  else 
  { 
    unsigned int index_i, index_j;
    
    for (unsigned int i = 0; i < m_BSplineShapeFunctionsIndexI.size(); i++)
    {
      if (m_BSplineShapeFunctionsIndexI(i) == -1)
      {
        shapeF[i] = 0.0;
      }
      else
      {
        index_j = static_cast<unsigned int>(static_cast<Float>(m_BSplineShapeFunctionsIndexI(i))/static_cast<Float>(m_BSplineOrder));
 index_i = m_BSplineShapeFunctionsIndexI(i) - index_j*m_BSplineOrder;  
        shapeF[i] = PolynomialType(m_BSplineShapeFunctions.get_row(index_i)).evaluate(u)*
                  PolynomialType(m_BSplineShapeFunctions.get_row(index_j+m_BSplineOrder)).evaluate(v);
      }  
      for (unsigned int j = 0; j < m_BSplineShapeFunctionsIndexJi[i]->size(); j++)
      {
        index_j = static_cast<unsigned int>(static_cast<Float>(m_BSplineShapeFunctionsIndexJi[i]->get(j))/static_cast<Float>(m_BSplineOrder));
 index_i = m_BSplineShapeFunctionsIndexJi[i]->get(j) - index_j*m_BSplineOrder; 
 
 shapeF[i] += m_BSplineShapeFunctionsJi_e[i]->get(j) *
                    PolynomialType(m_BSplineShapeFunctions.get_row(index_i)).evaluate(u)*
                   PolynomialType(m_BSplineShapeFunctions.get_row(index_j+m_BSplineOrder)).evaluate(v);
      }
    }
  }           
  return shapeF;
}

void
Element2DBSplinePatch
::ShapeFunctionDerivatives( const VectorType& pt, MatrixType& shapeD ) const
{
  Float a, b, u, v;
  a = m_localCoordinateBoundaries(0);
  b = m_localCoordinateBoundaries(1);
  u = 0.5*(b-a)*(pt[0]+1.0) + a;

  a = m_localCoordinateBoundaries(2);
  b = m_localCoordinateBoundaries(3);
  v = 0.5*(b-a)*(pt[1]+1.0) + a;

  /** functions at directions u and v.  */
  shapeD.set_size(NumberOfSpatialDimensions, m_NumberOfNodes);
  
  if (!m_ExtendedBSplineStatus)   
  {
    for (unsigned int j = 0; j < m_BSplineOrder; j++)
       for (unsigned int i = 0; i < m_BSplineOrder; i++)
          shapeD[0][j*m_BSplineOrder+i] = PolynomialType(m_BSplineShapeFunctions.get_row(i               )).devaluate(u)*
                               PolynomialType(m_BSplineShapeFunctions.get_row(j+m_BSplineOrder)).evaluate(v);

    for (unsigned int j = 0; j < m_BSplineOrder; j++)
       for (unsigned int i = 0; i < m_BSplineOrder; i++)
          shapeD[1][j*m_BSplineOrder+i] = PolynomialType(m_BSplineShapeFunctions.get_row(i               )).evaluate(u)*
                            PolynomialType(m_BSplineShapeFunctions.get_row(j+m_BSplineOrder)).devaluate(v);
  } 
  else 
  { 
    unsigned int index_i, index_j;
    for (unsigned int i = 0; i < m_BSplineShapeFunctionsIndexI.size(); i++)
    {
      if (m_BSplineShapeFunctionsIndexI(i) != -1)
      {
        shapeD[0][i] = shapeD[1][i] = 0.0;
      }
      else
      {
        index_j = static_cast<unsigned int>(static_cast<Float>(m_BSplineShapeFunctionsIndexI(i))/static_cast<Float>(m_BSplineOrder));
 index_i = m_BSplineShapeFunctionsIndexI(i) - index_j*m_BSplineOrder;  
 shapeD[0][i] = PolynomialType(m_BSplineShapeFunctions.get_row(index_i)).devaluate(u)*
                       PolynomialType(m_BSplineShapeFunctions.get_row(index_j+m_BSplineOrder)).evaluate(v);
 shapeD[1][i] = PolynomialType(m_BSplineShapeFunctions.get_row(index_i)).evaluate(u)*
                       PolynomialType(m_BSplineShapeFunctions.get_row(index_j+m_BSplineOrder)).devaluate(v);
      }      
      for (unsigned int j = 0; j < m_BSplineShapeFunctionsIndexJi[i]->size(); j++)
      {
        index_j = static_cast<unsigned int>(static_cast<Float>(m_BSplineShapeFunctionsIndexJi[i]->get(j))/static_cast<Float>(m_BSplineOrder));
 index_i = m_BSplineShapeFunctionsIndexJi[i]->get(j) - index_j*m_BSplineOrder; 
 shapeD[0][i] += m_BSplineShapeFunctionsJi_e[i]->get(j) *
                       PolynomialType(m_BSplineShapeFunctions.get_row(index_i)).devaluate(u)*
                      PolynomialType(m_BSplineShapeFunctions.get_row(index_j+m_BSplineOrder)).evaluate(v);
 shapeD[1][i] += m_BSplineShapeFunctionsJi_e[i]->get(j) *
                       PolynomialType(m_BSplineShapeFunctions.get_row(index_i)).evaluate(u)*
                      PolynomialType(m_BSplineShapeFunctions.get_row(index_j+m_BSplineOrder)).devaluate(v);
      }
    }
  }
  shapeD *= 0.5*(b-a);           
}

const Element2DBSplinePatch::VectorType
Element2DBSplinePatch
::GetCornerCoordinates(unsigned int n) const
{
   VectorType localPt;
   localPt.set_size(NumberOfSpatialDimensions);
   
   switch (n) 
   {
      case 0:
  localPt[0] = -1.0;
  localPt[1] = -1.0;
  break;
      case 1:
  localPt[0] =  1.0;
  localPt[1] = -1.0;
  break;
      case 2:
  localPt[0] = 1.0;
  localPt[1] = 1.0;
  break;
      case 3:
  localPt[0] = -1.0;
  localPt[1] =  1.0;
  break;
    }  
    
    return this->GetGlobalFromLocalCoordinates(localPt);
}

bool
Element2DBSplinePatch
::GetLocalFromGlobalCoordinates(const VectorType& globalPt, VectorType& localPt) const
{

  // If we can assume that the grid is structured such that the global coordinates
  //  are functions of the separated local coordinates, i.e. x = x(u) and y = y(v)
  //  instead of x = x(u, v) and y = y(u, v) then we can solve for the local coordinates
  //  by finding the appropriate roots of two separate univariate polynomials.

  Float x, y;
  localPt.set_size(NumberOfSpatialDimensions);
  localPt.fill(0.0);
  
  VectorType coeffs(m_BSplineOrder, 0.0), roots0, roots1;
  Float val0 = 0.5*(m_localCoordinateBoundaries(0)+m_localCoordinateBoundaries(1));
  Float val1 = 0.5*(m_localCoordinateBoundaries(2)+m_localCoordinateBoundaries(3));
  Float val;  

  for (unsigned int j = 0; j < m_BSplineOrder; j++)
  {
     val = PolynomialType(m_BSplineShapeFunctions.get_row(j+m_BSplineOrder)).evaluate(val1);     
     for (unsigned int i = 0; i < m_BSplineOrder; i++) 
     {
        x = m_node[j*m_BSplineOrder+i]->GetCoordinates()[0];
        coeffs += m_BSplineShapeFunctions.get_row(i)*val*x;
     } 
  }
  coeffs(m_BSplineOrder-1) -= globalPt(0);
  while (coeffs.size() > 0) 
  {
     if (fabs(coeffs[0]) < 1e-10)  
        coeffs = coeffs.extract(coeffs.size()-1, 1);
     else 
        break;
  }  
  vnl_rpoly_roots rps0(coeffs);
  roots0 = rps0.realroots(1e-10);  
  if (roots0.size() == 0)  return false;


  coeffs.set_size(m_BSplineOrder);  
  coeffs.fill(0.0);
  for (unsigned int i = 0; i < m_BSplineOrder; i++)
  {
     val = PolynomialType(m_BSplineShapeFunctions.get_row(i)).evaluate(val0);     
     for (unsigned int j = 0; j < m_BSplineOrder; j++) 
     {
        y = m_node[j*m_BSplineOrder+i]->GetCoordinates()[1];
        coeffs += m_BSplineShapeFunctions.get_row(j+m_BSplineOrder)*val*y; 
     } 
  }
  coeffs(m_BSplineOrder-1) -= globalPt(1);
  while (coeffs.size() > 0) 
  {
     if (fabs(coeffs[0]) < 1e-10)  
        coeffs = coeffs.extract(coeffs.size()-1, 1);
     else 
        break;
  }  
  vnl_rpoly_roots rps1(coeffs);
  roots1 = rps1.realroots(1e-10);  
  if (roots1.size() == 0)  return false;
     
          
  localPt[0] = roots0[0];
  for (unsigned int i = 1; i < roots0.size(); i++)
     if (fabs(val0-roots0[i]) < fabs(val0-localPt[0]))
        localPt[0] = roots0[i];
 
  localPt[1] = roots1[0];
  for (unsigned int i = 1; i < roots1.size(); i++)
     if (fabs(val1-roots1[i]) < fabs(val1-localPt[1]))
        localPt[1] = roots1[i];
  
  bool isInside = false;  
  // Rescale and translate the local points to be in the range [-1.0, 1.0]
  localPt[0] = 2.0*(localPt[0]-m_localCoordinateBoundaries(0))/(m_localCoordinateBoundaries(1)-m_localCoordinateBoundaries(0)) - 1.0;
  localPt[1] = 2.0*(localPt[1]-m_localCoordinateBoundaries(2))/(m_localCoordinateBoundaries(3)-m_localCoordinateBoundaries(2)) - 1.0;

  if (localPt[0] >= -(1.0+1.e-10) && localPt[0] <= (1.0+1.e-10) && localPt[1] >= -(1.0+1.e-10) && localPt[1] <= (1.0+1.e-10))
     isInside = true;

  return isInside;

  // The code below works for the general case

/*
  Float x, y;
  
  localPt.set_size(NumberOfSpatialDimensions);
  localPt.fill(0.0);
  
  VectorType coeffs(m_BSplineOrder*m_BSplineOrder, 0.0);
  vnl_matrix<int> expnts(m_BSplineOrder*m_BSplineOrder, 2);
  vcl_vector<vnl_real_npolynomial*> npoly_vector(2);  
  vcl_vector<vnl_vector<double>*> roots;

  for (unsigned int i = 0; i < m_BSplineOrder; i++)
  {
     for (unsigned int j = 0; j < m_BSplineOrder; j++)
     {
        expnts(i*m_BSplineOrder+j, 0) = m_BSplineOrder-i-1;
        expnts(i*m_BSplineOrder+j, 1) = m_BSplineOrder-j-1;             
     } 
  }   
  npoly_vector[0] = new vnl_real_npolynomial;  npoly_vector[0]->set(coeffs, expnts);
  npoly_vector[1] = new vnl_real_npolynomial;  npoly_vector[1]->set(coeffs, expnts);
  
  for (unsigned int j = 0; j < m_BSplineOrder; j++)
  {
     for (unsigned int i = 0; i < m_BSplineOrder; i++)
     {
        for (unsigned int n = 0; n < m_BSplineOrder; n++)
           for (unsigned int m = 0; m < m_BSplineOrder; m++)
       coeffs(n*m_BSplineOrder+m) = m_BSplineShapeFunctions(j+m_BSplineOrder, m)*m_BSplineShapeFunctions(i, n); 

        x = m_node[j*m_BSplineOrder+i]->GetCoordinates()[0];
        y = m_node[j*m_BSplineOrder+i]->GetCoordinates()[1];

 npoly_vector[0]->set(npoly_vector[0]->coefficients()+x*coeffs, expnts);
 npoly_vector[1]->set(npoly_vector[1]->coefficients()+y*coeffs, expnts);
     }   
  } 
  coeffs = npoly_vector[0]->coefficients();
  coeffs.put(coeffs.size()-1, coeffs.get(coeffs.size()-1) - globalPt[0]);
  npoly_vector[0]->set(coeffs, expnts);

  coeffs = npoly_vector[1]->coefficients();
  coeffs.put(coeffs.size()-1, coeffs.get(coeffs.size()-1) - globalPt[1]);
  npoly_vector[1]->set(coeffs, expnts);
  
  vnl_rnpoly_solve rps(npoly_vector);
  roots = rps.realroots();

  delete npoly_vector[0];
  delete npoly_vector[1];

  bool isInside = false;  
  for (unsigned int i = 0; i < roots.size(); i++)
  {
     localPt = *(roots[i]);
     
     // Rescale and translate the local points to be in the range [-1.0, 1.0]
     localPt[0] = 2.0*(localPt[0]-m_localCoordinateBoundaries(0))/(m_localCoordinateBoundaries(1)-m_localCoordinateBoundaries(0)) - 1.0;
     localPt[1] = 2.0*(localPt[1]-m_localCoordinateBoundaries(2))/(m_localCoordinateBoundaries(3)-m_localCoordinateBoundaries(2)) - 1.0;

     if (localPt[0] >= -1.0 && localPt[0] <= 1.0 && localPt[1] >= -1.0 && localPt[1] <= 1.0)
     {     
        isInside = true;
 break;
     } 
  }    

  return isInside;
*/  
}

/*
 * Draw the element on device context pDC.
 */
#ifdef FEM_BUILD_VISUALIZATION
void
Element2DBSplinePatch
::Draw(CDC* pDC, Solution::ConstPointer sol) const
{  
  VectorType localPt(NumberOfSpatialDimensions);

  localPt[0] = -1.0;  localPt[1] = -1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->MoveTo(globalPt[0], globalPt[1]);  

  localPt[0] =  1.0;  localPt[1] = -1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->LineTo(globalPt[0], globalPt[1]);

  localPt[0] =  1.0;  localPt[1] =  1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->LineTo(globalPt[0], globalPt[1]);

  localPt[0] = -1.0;  localPt[1] =  1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->LineTo(globalPt[0], globalPt[1]);

  localPt[0] = -1.0;  localPt[1] = -1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->LineTo(globalPt[0], globalPt[1]);  
}
#endif

}} // end namespace itk::fem

/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMElement3DBSplinePatch.cxx,v $
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

#include "itkFEMElement3DBSplinePatch.h"
#include "itkFEMBSplineShapeFunctions.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_real_polynomial.h"
#include "vnl/algo/vnl_rpoly_roots.h"
#include "vnl/algo/vnl_rnpoly_solve.h"

namespace itk {
namespace fem {


void
Element3DBSplinePatch
::GetIntegrationPointAndWeight(unsigned int i, VectorType& pt, Float& w, unsigned int order) const
{
  // FIXME: range checking

  // default integration order
  //if (order==0) { order=DefaultIntegrationOrder; }
  //order = static_cast<unsigned int>(0.5*static_cast<float>(m_BSplineOrder)+0.5);
  order = m_BSplineOrder;

  pt.set_size(NumberOfSpatialDimensions);
  pt[0] = gaussPoint[order][i%order];
  pt[1] = gaussPoint[order][(i/order)%order];
  pt[2] = gaussPoint[order][(i/(order*order))];

  w=gaussWeight[order][i%order]*gaussWeight[order][(i/order)%order]*gaussWeight[order][(i/(order*order))];
}

unsigned int
Element3DBSplinePatch
::GetNumberOfIntegrationPoints(unsigned int order) const
{
  // FIXME: range checking

  //if (order==0) { order=DefaultIntegrationOrder; }
  //order = static_cast<unsigned int>(0.5*static_cast<float>(m_BSplineOrder)+0.5);
  order = m_BSplineOrder;

  return order*order*order;
}

void
Element3DBSplinePatch
::GenerateShapeFunctions()
{
  unsigned int m = 2*m_BSplineOrder;
  vnl_vector<double> knots(m);
  
  for (unsigned int i = 0; i < m; i++)    
    knots[i] = static_cast<double>(i)/static_cast<double>(m-1);

  BSplineShapeFunctions shape;  shape.SetKnots(knots);
  
  vnl_matrix<double> shapeFunctions(3*m_BSplineOrder, m_BSplineOrder);

  for(unsigned int f = 0; f < m_BSplineOrder; f++) 
  {
    shapeFunctions.set_row(f                 , 
      (shape.GenerateBSplineShapeFunction(m_BSplineOrder, f, m_BSplineOrder-1)).coefficients());
    shapeFunctions.set_row(f+m_BSplineOrder  , 
      (shape.GenerateBSplineShapeFunction(m_BSplineOrder, f, m_BSplineOrder-1)).coefficients());   
    shapeFunctions.set_row(f+m_BSplineOrder*2, 
      (shape.GenerateBSplineShapeFunction(m_BSplineOrder, f, m_BSplineOrder-1)).coefficients());   
  }
  this->SetBSplineShapeFunctions(shapeFunctions);      
}

Element3DBSplinePatch::VectorType
Element3DBSplinePatch
::ShapeFunctions( const VectorType& pt ) const
{
  Float a, b, u, v, w;
  a = m_localCoordinateBoundaries(0);
  b = m_localCoordinateBoundaries(1);
  u = 0.5*(b-a)*(pt[0]+1.0) + a;

  a = m_localCoordinateBoundaries(2);
  b = m_localCoordinateBoundaries(3);
  v = 0.5*(b-a)*(pt[1]+1.0) + a;

  a = m_localCoordinateBoundaries(4);
  b = m_localCoordinateBoundaries(5);
  w = 0.5*(b-a)*(pt[2]+1.0) + a;
   
  VectorType shapeF(m_NumberOfNodes);

  if (!m_ExtendedBSplineStatus)   
  {
    for (unsigned int k = 0; k < m_BSplineOrder; k++)
      for (unsigned int j = 0; j < m_BSplineOrder; j++)
        for (unsigned int i = 0; i < m_BSplineOrder; i++)
          shapeF[k*m_BSplineOrder*m_BSplineOrder+j*m_BSplineOrder+i] = 
                                vnl_real_polynomial(m_BSplineShapeFunctions.get_row(i)).evaluate(u)*
                                   vnl_real_polynomial(m_BSplineShapeFunctions.get_row(j+m_BSplineOrder)).evaluate(v)*
           vnl_real_polynomial(m_BSplineShapeFunctions.get_row(k+2*m_BSplineOrder)).evaluate(w);
  } 
  else 
  { 
    unsigned int index_i, index_j, index_k;
    for (unsigned int i = 0; i < m_BSplineShapeFunctionsIndexI.size(); i++)
    {
      if (m_BSplineShapeFunctionsIndexI(i) == -1)
      {
        shapeF[i] = 0.0;
      }
      else
      {
        index_k = static_cast<unsigned int>(static_cast<Float>(m_BSplineShapeFunctionsIndexI(i))/static_cast<Float>(m_BSplineOrder*m_BSplineOrder));
        index_j = static_cast<unsigned int>(static_cast<Float>(m_BSplineShapeFunctionsIndexI(i)-index_k*m_BSplineOrder)/static_cast<Float>(m_BSplineOrder));
 index_i = m_BSplineShapeFunctionsIndexI(i) - index_k*m_BSplineOrder*m_BSplineOrder - index_j*m_BSplineOrder;  
        shapeF[i] = vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_i)).evaluate(u)*
                vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_j+m_BSplineOrder)).evaluate(v)*
      vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_k+2*m_BSplineOrder)).evaluate(w);
      }      
      for (unsigned int j = 0; j < m_BSplineShapeFunctionsIndexJi[i]->size(); j++)
      {
        index_k = static_cast<unsigned int>(static_cast<Float>(m_BSplineShapeFunctionsIndexJi[i]->get(j))/static_cast<Float>(m_BSplineOrder*m_BSplineOrder));
        index_j = static_cast<unsigned int>(static_cast<Float>(m_BSplineShapeFunctionsIndexJi[i]->get(j)-index_k*m_BSplineOrder)/static_cast<Float>(m_BSplineOrder));
 index_i = (m_BSplineShapeFunctionsIndexJi[i]->get(j)) - index_k*m_BSplineOrder*m_BSplineOrder - index_j*m_BSplineOrder;  
        shapeF[i] += vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_i)).evaluate(u)*
                 vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_j+m_BSplineOrder)).evaluate(v)*
       vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_k+2*m_BSplineOrder)).evaluate(w);
      }
    }
  }           
  return shapeF;
}

void
Element3DBSplinePatch
::ShapeFunctionDerivatives( const VectorType& pt, MatrixType& shapeD ) const
{

  Float a, b, u, v, w;
  a = m_localCoordinateBoundaries(0);
  b = m_localCoordinateBoundaries(1);
  u = 0.5*(b-a)*(pt[0]+1.0) + a;

  a = m_localCoordinateBoundaries(2);
  b = m_localCoordinateBoundaries(3);
  v = 0.5*(b-a)*(pt[1]+1.0) + a;

  a = m_localCoordinateBoundaries(4);
  b = m_localCoordinateBoundaries(5);
  w = 0.5*(b-a)*(pt[2]+1.0) + a;

  /** functions at directions r and s.  */
  shapeD.set_size(NumberOfSpatialDimensions, m_NumberOfNodes);
  
  if (!m_ExtendedBSplineStatus)   
  {
    for (unsigned int k = 0; k < m_BSplineOrder; k++)
       for (unsigned int j = 0; j < m_BSplineOrder; j++)
          for (unsigned int i = 0; i < m_BSplineOrder; i++)
             shapeD[0][k*m_BSplineOrder*m_BSplineOrder+j*m_BSplineOrder+i] = 
                             vnl_real_polynomial(m_BSplineShapeFunctions.get_row(i)).devaluate(u)*
                             vnl_real_polynomial(m_BSplineShapeFunctions.get_row(j+m_BSplineOrder)).evaluate(v)*
                             vnl_real_polynomial(m_BSplineShapeFunctions.get_row(k+m_BSplineOrder*2)).evaluate(w);

    for (unsigned int k = 0; k < m_BSplineOrder; k++)
       for (unsigned int j = 0; j < m_BSplineOrder; j++)
          for (unsigned int i = 0; i < m_BSplineOrder; i++)
             shapeD[1][k*m_BSplineOrder*m_BSplineOrder+j*m_BSplineOrder+i] = 
                             vnl_real_polynomial(m_BSplineShapeFunctions.get_row(i)).evaluate(u)*
                             vnl_real_polynomial(m_BSplineShapeFunctions.get_row(j+m_BSplineOrder)).devaluate(v)*
                             vnl_real_polynomial(m_BSplineShapeFunctions.get_row(k+m_BSplineOrder*2)).evaluate(w);

    for (unsigned int k = 0; k < m_BSplineOrder; k++)
       for (unsigned int j = 0; j < m_BSplineOrder; j++)
          for (unsigned int i = 0; i < m_BSplineOrder; i++)
             shapeD[2][k*m_BSplineOrder*m_BSplineOrder+j*m_BSplineOrder+i] = 
                             vnl_real_polynomial(m_BSplineShapeFunctions.get_row(i)).evaluate(u)*
                             vnl_real_polynomial(m_BSplineShapeFunctions.get_row(j+m_BSplineOrder)).evaluate(v)*
                             vnl_real_polynomial(m_BSplineShapeFunctions.get_row(k+m_BSplineOrder*2)).devaluate(w);
  }
  else 
  { 
    unsigned int index_i, index_j, index_k;
    for (unsigned int i = 0; i < m_BSplineShapeFunctionsIndexI.size(); i++)
    {
      if (m_BSplineShapeFunctionsIndexI(i) != -1)
      {
        shapeD[0][i] = shapeD[1][i] = shapeD[2][i] = 0.0;
      }
      else
      {
        index_k = static_cast<unsigned int>(static_cast<Float>(m_BSplineShapeFunctionsIndexI(i))/static_cast<Float>(m_BSplineOrder*m_BSplineOrder));
        index_j = static_cast<unsigned int>(static_cast<Float>(m_BSplineShapeFunctionsIndexI(i)-index_k*m_BSplineOrder)/static_cast<Float>(m_BSplineOrder));
 index_i = m_BSplineShapeFunctionsIndexI(i) - index_k*m_BSplineOrder*m_BSplineOrder - index_j*m_BSplineOrder;  
        shapeD[0][i] = vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_i)).devaluate(u)*
                   vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_j+m_BSplineOrder)).evaluate(v)*
         vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_k+2*m_BSplineOrder)).evaluate(w);
        shapeD[1][i] = vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_i)).evaluate(u)*
                   vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_j+m_BSplineOrder)).devaluate(v)*
         vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_k+2*m_BSplineOrder)).evaluate(w);
        shapeD[2][i] = vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_i)).evaluate(u)*
                   vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_j+m_BSplineOrder)).evaluate(v)*
         vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_k+2*m_BSplineOrder)).devaluate(w);
      }      
      for (unsigned int j = 0; j < m_BSplineShapeFunctionsIndexJi[i]->size(); j++)
      {
        index_k = static_cast<unsigned int>(static_cast<Float>((m_BSplineShapeFunctionsIndexJi[i]->get(j)))/static_cast<Float>(m_BSplineOrder*m_BSplineOrder));
        index_j = static_cast<unsigned int>(static_cast<Float>((m_BSplineShapeFunctionsIndexJi[i]->get(j))-index_k*m_BSplineOrder)/static_cast<Float>(m_BSplineOrder));
 index_i = (m_BSplineShapeFunctionsIndexJi[i]->get(j)) - index_k*m_BSplineOrder*m_BSplineOrder - index_j*m_BSplineOrder;  
        shapeD[0][i] += vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_i)).devaluate(u)*
                    vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_j+m_BSplineOrder)).evaluate(v)*
   vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_k+2*m_BSplineOrder)).evaluate(w);
        shapeD[1][i] += vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_i)).evaluate(u)*
                    vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_j+m_BSplineOrder)).devaluate(v)*
   vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_k+2*m_BSplineOrder)).evaluate(w);
        shapeD[2][i] += vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_i)).evaluate(u)*
             vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_j+m_BSplineOrder)).evaluate(v)*
   vnl_real_polynomial(m_BSplineShapeFunctions.get_row(index_k+2*m_BSplineOrder)).devaluate(w);
      }
    }
  }           
  shapeD *= 0.5*(b-a);           
}

const Element3DBSplinePatch::VectorType
Element3DBSplinePatch
::GetCornerCoordinates(unsigned int n) const
{
   VectorType localPt;
   localPt.set_size(NumberOfSpatialDimensions);
   
   switch (n) 
   {
      case 0:
  localPt[0] = -1.0;
  localPt[1] = -1.0;
  localPt[2] = -1.0;
  break;
      case 1:
  localPt[0] =  1.0;
  localPt[1] = -1.0;
  localPt[2] = -1.0;
  break;
      case 2:
  localPt[0] =  1.0;
  localPt[1] =  1.0;
  localPt[2] = -1.0;
  break;
      case 3:
  localPt[0] = -1.0;
  localPt[1] =  1.0;
  localPt[2] = -1.0;
  break;
      case 4:
  localPt[0] = -1.0;
  localPt[1] = -1.0;
  localPt[2] =  1.0;
  break;
      case 5:
  localPt[0] =  1.0;
  localPt[1] = -1.0;
  localPt[2] =  1.0;
  break;
      case 6:
  localPt[0] =  1.0;
  localPt[1] =  1.0;
  localPt[2] =  1.0;
  break;
      case 7:
  localPt[0] = -1.0;
  localPt[1] =  1.0;
  localPt[2] =  1.0;
  break;
    }  
    
    return this->GetGlobalFromLocalCoordinates(localPt);
}

bool
Element3DBSplinePatch
::GetLocalFromGlobalCoordinates( const VectorType& globalPt, VectorType& localPt) const
{
   
  // If we can assume that the grid is structured such that the global coordinates
  //  are functions of the separated local coordinates, i.e. x = x(u) and y = y(v)
  //  instead of x = x(u, v) and y = y(u, v) then we can solve for the local coordinates
  //  by finding the appropriate roots of two separate univariate polynomials.

  Float x, y, z;
  localPt.set_size(NumberOfSpatialDimensions);
  localPt.fill(0.0);
  
  VectorType coeffs(m_BSplineOrder, 0.0), roots0, roots1, roots2;
  Float val0 = 0.5*(m_localCoordinateBoundaries(0)+m_localCoordinateBoundaries(1));
  Float val1 = 0.5*(m_localCoordinateBoundaries(2)+m_localCoordinateBoundaries(3));
  Float val2 = 0.5*(m_localCoordinateBoundaries(4)+m_localCoordinateBoundaries(5));
  Float v1, v2;  

  for (unsigned int k = 0; k < m_BSplineOrder; k++)
  {
     v2 = vnl_real_polynomial(m_BSplineShapeFunctions.get_row(k+2*m_BSplineOrder)).evaluate(val2);     
     for (unsigned int j = 0; j < m_BSplineOrder; j++)
     {
        v1 = vnl_real_polynomial(m_BSplineShapeFunctions.get_row(j+m_BSplineOrder)).evaluate(val1);     
 for (unsigned int i = 0; i < m_BSplineOrder; i++) 
 {
           x = m_node[k*m_BSplineOrder*m_BSplineOrder+j*m_BSplineOrder+i]->GetCoordinates()[0];
           coeffs += m_BSplineShapeFunctions.get_row(i)*v1*v2*x;
 } 
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
     v2 = vnl_real_polynomial(m_BSplineShapeFunctions.get_row(i)).evaluate(val0);     
     for (unsigned int k = 0; k < m_BSplineOrder; k++)
     {
        v1 = vnl_real_polynomial(m_BSplineShapeFunctions.get_row(k+2*m_BSplineOrder)).evaluate(val2);     
 for (unsigned int j = 0; j < m_BSplineOrder; j++) 
 {
           y = m_node[k*m_BSplineOrder*m_BSplineOrder+j*m_BSplineOrder+i]->GetCoordinates()[1];
           coeffs += m_BSplineShapeFunctions.get_row(j+m_BSplineOrder)*v1*v2*y;
 } 
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

  coeffs.set_size(m_BSplineOrder);  
  coeffs.fill(0.0);
  
  for (unsigned int j = 0; j < m_BSplineOrder; j++)
  {
     v2 = vnl_real_polynomial(m_BSplineShapeFunctions.get_row(j+2*m_BSplineOrder)).evaluate(val1);     
     for (unsigned int i = 0; i < m_BSplineOrder; i++)
     {
        v1 = vnl_real_polynomial(m_BSplineShapeFunctions.get_row(i)).evaluate(val0);     
 for (unsigned int k = 0; k < m_BSplineOrder; k++) 
 {
           z = m_node[k*m_BSplineOrder*m_BSplineOrder+j*m_BSplineOrder+i]->GetCoordinates()[2];
           coeffs += m_BSplineShapeFunctions.get_row(k+2*m_BSplineOrder)*v1*v2*z;
 } 
     }
  }   
  coeffs(m_BSplineOrder-1) -= globalPt(2);
  while (coeffs.size() > 0) 
  {
     if (fabs(coeffs[0]) < 1e-10)  
        coeffs = coeffs.extract(coeffs.size()-1, 1);
     else 
        break;
  }  
  vnl_rpoly_roots rps2(coeffs);
  roots2 = rps2.realroots(1e-10);  
  if (roots2.size() == 0)  return false;

  coeffs.set_size(m_BSplineOrder);  
  coeffs.fill(0.0);
     
          
  localPt[0] = roots0[0];
  for (unsigned int i = 1; i < roots0.size(); i++)
     if (fabs(val0-roots0[i]) < fabs(val0-localPt[0]))
        localPt[0] = roots0[i];
 
  localPt[1] = roots1[0];
  for (unsigned int i = 1; i < roots1.size(); i++)
     if (fabs(val1-roots1[i]) < fabs(val1-localPt[1]))
        localPt[1] = roots1[i];
  
  localPt[2] = roots2[0];
  for (unsigned int i = 1; i < roots2.size(); i++)
     if (fabs(val2-roots2[i]) < fabs(val2-localPt[2]))
        localPt[2] = roots2[i];

  bool isInside = false;  
  // Rescale and translate the local points to be in the range [-1.0, 1.0]
  localPt[0] = 2.0*(localPt[0]-m_localCoordinateBoundaries(0))/(m_localCoordinateBoundaries(1)-m_localCoordinateBoundaries(0)) - 1.0;
  localPt[1] = 2.0*(localPt[1]-m_localCoordinateBoundaries(2))/(m_localCoordinateBoundaries(3)-m_localCoordinateBoundaries(2)) - 1.0;
  localPt[2] = 2.0*(localPt[2]-m_localCoordinateBoundaries(4))/(m_localCoordinateBoundaries(5)-m_localCoordinateBoundaries(4)) - 1.0;

  if (localPt[0] >= -(1.0+1.e-10) && localPt[0] <= (1.0+1.e-10) 
      && localPt[1] >= -(1.0+1.e-10) && localPt[1] <= (1.0+1.e-10) 
      && localPt[2] >= -(1.0+1.e-10) && localPt[2] <= (1.0+1.e-10))
     isInside = true;

  return isInside;
  
  
  
   // The code below works for the general case

/*
  Float x, y, z;
  
  localPt.set_size(NumberOfSpatialDimensions);
  localPt.fill(0.0);
  
  VectorType coeffs(m_BSplineOrder*m_BSplineOrder*m_BSplineOrder, 0.0);
  vnl_matrix<int> expnts(m_BSplineOrder*m_BSplineOrder*m_BSplineOrder, 3);
  vcl_vector<vnl_real_npolynomial*> npoly_vector(NumberOfSpatialDimensions);  
  vcl_vector<vnl_vector<double>*> roots;
  
  unsigned int index;

  for (unsigned int i = 0; i < m_BSplineOrder; i++)
  {
     for (unsigned int j = 0; j < m_BSplineOrder; j++)
     {
        for (unsigned int k = 0; k < m_BSplineOrder; k++)
        {
    index = i*m_BSplineOrder*m_BSplineOrder+j*m_BSplineOrder+k;
           expnts(index, 0) = m_BSplineOrder-i-1;
           expnts(index, 1) = m_BSplineOrder-j-1;             
           expnts(index, 2) = m_BSplineOrder-k-1;             
 }   
     } 
  }   
  npoly_vector[0] = new vnl_real_npolynomial;  npoly_vector[0]->set(coeffs, expnts);
  npoly_vector[1] = new vnl_real_npolynomial;  npoly_vector[1]->set(coeffs, expnts);
  npoly_vector[2] = new vnl_real_npolynomial;  npoly_vector[2]->set(coeffs, expnts);

  for (unsigned int k = 0; k < m_BSplineOrder; k++)
  {
     for (unsigned int j = 0; j < m_BSplineOrder; j++)
     {
 for (unsigned int i = 0; i < m_BSplineOrder; i++)
 {
           for (unsigned int n = 0; n < m_BSplineOrder; n++)
              for (unsigned int m = 0; m < m_BSplineOrder; m++)
          for (unsigned int l = 0; l < m_BSplineOrder; l++)
      coeffs(n*m_BSplineOrder*m_BSplineOrder+m*m_BSplineOrder+l) 
         =  m_BSplineShapeFunctions(k+2*m_BSplineOrder, l)
           *m_BSplineShapeFunctions(j+m_BSplineOrder, m)
    *m_BSplineShapeFunctions(i, n); 

           x = m_node[k*m_BSplineOrder*m_BSplineOrder+j*m_BSplineOrder+i]->GetCoordinates()[0];
           y = m_node[k*m_BSplineOrder*m_BSplineOrder+j*m_BSplineOrder+i]->GetCoordinates()[1];
           z = m_node[k*m_BSplineOrder*m_BSplineOrder+j*m_BSplineOrder+i]->GetCoordinates()[2];

    npoly_vector[0]->set(npoly_vector[0]->coefficients()+x*coeffs, expnts);
    npoly_vector[1]->set(npoly_vector[1]->coefficients()+y*coeffs, expnts);
    npoly_vector[2]->set(npoly_vector[2]->coefficients()+z*coeffs, expnts);
 }    
     }   
  } 
  coeffs = npoly_vector[0]->coefficients();
  coeffs(coeffs.size()-1) = globalPt[0];
  npoly_vector[0]->set(coeffs, expnts);

  coeffs = npoly_vector[1]->coefficients();
  coeffs(coeffs.size()-1) = globalPt[1];
  npoly_vector[1]->set(coeffs, expnts);
  
  coeffs = npoly_vector[2]->coefficients();
  coeffs(coeffs.size()-1) = globalPt[2];
  npoly_vector[2]->set(coeffs, expnts);
  
  vnl_rnpoly_solve rps(npoly_vector);
  roots = rps.realroots();

  delete npoly_vector[0];
  delete npoly_vector[1];
  delete npoly_vector[2];

  bool isInside = false;  
  for (unsigned int i = 0; i < roots.size(); i++)
  {
     localPt = *(roots[i]);
     // Rescale and translate the local points to be in the range [-1.0, 1.0]
     localPt[0] = 2.0*(localPt[0]-m_localCoordinateBoundaries(0))/(m_localCoordinateBoundaries(1)-m_localCoordinateBoundaries(0)) - 1.0;
     localPt[1] = 2.0*(localPt[1]-m_localCoordinateBoundaries(2))/(m_localCoordinateBoundaries(3)-m_localCoordinateBoundaries(2)) - 1.0;
     localPt[2] = 2.0*(localPt[2]-m_localCoordinateBoundaries(4))/(m_localCoordinateBoundaries(5)-m_localCoordinateBoundaries(4)) - 1.0;
     
     if (localPt[0] > -1.0 && localPt[0] < 1.0 && localPt[1] > -1.0 && localPt[1] < 1.0 && localPt[2] > -1.0 && localPt[2] < 1.0)
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
Element3DBSplinePatch
::Draw(CDC* pDC, Solution::ConstPointer sol) const
{
  VectorType localPt(NumberOfSpatialDimensions);

  // Draw Bottom Quadrilateral

  localPt[0] = -1.0;  localPt[1] = -1.0;  localPt[2] = -1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->MoveTo(globalPt[0], globalPt[1]);  

  localPt[0] =  1.0;  localPt[1] = -1.0;  localPt[2] = -1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->LineTo(globalPt[0], globalPt[1]);

  localPt[0] =  1.0;  localPt[1] =  1.0;  localPt[2] = -1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->LineTo(globalPt[0], globalPt[1]);

  localPt[0] = -1.0;  localPt[1] =  1.0;  localPt[2] = -1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->LineTo(globalPt[0], globalPt[1]);

  localPt[0] = -1.0;  localPt[1] = -1.0;  localPt[2] = -1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->LineTo(globalPt[0], globalPt[1]);  

  // Draw Top Quadrilateral

  localPt[0] = -1.0;  localPt[1] = -1.0;  localPt[2] =  1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->MoveTo(globalPt[0], globalPt[1]);  

  localPt[0] =  1.0;  localPt[1] = -1.0;  localPt[2] =  1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->LineTo(globalPt[0], globalPt[1]);

  localPt[0] =  1.0;  localPt[1] =  1.0;  localPt[2] =  1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->LineTo(globalPt[0], globalPt[1]);

  localPt[0] = -1.0;  localPt[1] =  1.0;  localPt[2] =  1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->LineTo(globalPt[0], globalPt[1]);

  localPt[0] = -1.0;  localPt[1] = -1.0;  localPt[2] =  1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->LineTo(globalPt[0], globalPt[1]);  

  // Draw Sides

  localPt[0] = -1.0;  localPt[1] = -1.0;  localPt[2] = -1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->MoveTo(globalPt[0], globalPt[1]);  

  localPt[0] = -1.0;  localPt[1] = -1.0;  localPt[2] =  1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->LineTo(globalPt[0], globalPt[1]);  

  localPt[0] =  1.0;  localPt[1] = -1.0;  localPt[2] = -1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->MoveTo(globalPt[0], globalPt[1]);  

  localPt[0] =  1.0;  localPt[1] = -1.0;  localPt[2] =  1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->LineTo(globalPt[0], globalPt[1]);  

  localPt[0] =  1.0;  localPt[1] =  1.0;  localPt[2] = -1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->MoveTo(globalPt[0], globalPt[1]);  

  localPt[0] =  1.0;  localPt[1] =  1.0;  localPt[2] =  1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->LineTo(globalPt[0], globalPt[1]);  

  localPt[0] = -1.0;  localPt[1] =  1.0;  localPt[2] = -1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->MoveTo(globalPt[0], globalPt[1]);  

  localPt[0] = -1.0;  localPt[1] =  1.0;  localPt[2] =  1.0;
  globalPt = this->GetSolutionGlobalFromLocalCoordinates(localPt, sol)*DC_Scale;
  pDC->LineTo(globalPt[0], globalPt[1]);  
}
#endif




}} // end namespace itk::fem

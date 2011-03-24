/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMElement3DC0LinearHexahedron.cxx,v $
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

#include "itkFEMElement3DC0LinearHexahedron.h"
#include "vnl/vnl_math.h"
#include "vnl/algo/vnl_rnpoly_solve.h"
#include "vcl_vector.h"

namespace itk {
namespace fem {




void
Element3DC0LinearHexahedron
::GetIntegrationPointAndWeight(unsigned int i, VectorType& pt, Float& w, unsigned int order) const
{
  // FIXME: range checking

  // default integration order=2
  if (order==0) { order=2; }

  pt.set_size(3);
  pt[0] = gaussPoint[order][i%order];
  pt[1] = gaussPoint[order][(i/order)%order];
  pt[2] = gaussPoint[order][(i/(order*order))];

  w=gaussWeight[order][i%order]*gaussWeight[order][(i/order)%order]*gaussWeight[order][(i/(order*order))];

}




unsigned int
Element3DC0LinearHexahedron
::GetNumberOfIntegrationPoints(unsigned int order) const
{
  // FIXME: range checking

  // default integration order=2
  if (order==0) { order=2; }

  return order*order*order;
}




Element3DC0LinearHexahedron::VectorType
Element3DC0LinearHexahedron
::ShapeFunctions( const VectorType& pt ) const
{
  /* Linear hexahedral element has eight shape functions  */
  VectorType shapeF(8);

  /*
   * Linear hexahedral element has local coordinates
   *  (-1,-1,-1), (1,-1,-1), (1,1,-1), (-1,1,-1), (-1,-1,1), (1,-1,1), (1,1,1), (-1,1,1)
   */
   
   Float sign[8][3] = { {-1.0,-1.0,-1.0}, {1.0,-1.0,-1.0}, {1.0,1.0,-1.0}, {-1.0,1.0,-1.0}, 
                        {-1.0,-1.0, 1.0}, {1.0,-1.0, 1.0}, {1.0,1.0, 1.0}, {-1.0,1.0, 1.0} };
    
  /* given local point x=(r,s,t), where -1 <= r,s,t <= 1 and */
  
  for (unsigned int i = 0; i < 8; i++)
  {
     shapeF[i] = 0.125;
     for (unsigned int j = 0; j < 3; j++)
        shapeF[i] *= (1.0 + sign[i][j]*pt[j]);
  }

  return shapeF;
}



void
Element3DC0LinearHexahedron
::ShapeFunctionDerivatives( const VectorType& pt, MatrixType& shapeD ) const
{
  /** functions at directions r and s.  */
  shapeD.set_size(3,8);

  /*
   * Linear hexahedral element has local coordinates
   *  (-1,-1,-1), (1,-1,-1), (1,1,-1), (-1,1,-1), (-1,-1,1), (1,-1,1), (1,1,1), (-1,1,1)
   */
   
   Float sign[8][3] = { {-1.0,-1.0,-1.0}, {1.0,-1.0,-1.0}, {1.0,1.0,-1.0}, {-1.0,1.0,-1.0}, 
                        {-1.0,-1.0, 1.0}, {1.0,-1.0, 1.0}, {1.0,1.0, 1.0}, {-1.0,1.0, 1.0} };
    
  for (unsigned int i = 0; i < 8; i++)
     shapeD.put(0, i, sign[i][0] * (1.0 + sign[i][1]*pt[1]) * (1.0 + sign[i][2]*pt[2]) * 0.125);

  for (unsigned int i = 0; i < 8; i++)
     shapeD.put(1, i, sign[i][1] * (1.0 + sign[i][2]*pt[2]) * (1.0 + sign[i][0]*pt[0]) * 0.125);

  for (unsigned int i = 0; i < 8; i++)
     shapeD.put(2, i, sign[i][2] * (1.0 + sign[i][0]*pt[0]) * (1.0 + sign[i][1]*pt[1]) * 0.125);

}


bool
Element3DC0LinearHexahedron
::GetLocalFromGlobalCoordinates( const VectorType& globalPt , VectorType& localPt ) const
{

  Float sign[8][3] = { {-1.0,-1.0,-1.0}, {1.0,-1.0,-1.0}, {1.0,1.0,-1.0}, {-1.0,1.0,-1.0}, 
                       {-1.0,-1.0, 1.0}, {1.0,-1.0, 1.0}, {1.0,1.0, 1.0}, {-1.0,1.0, 1.0} };

  Float x, y, z;
  
  localPt.set_size(3);
  localPt.fill(0.0);
  
  VectorType coeffs(8, 0.0);
  vnl_matrix<unsigned int> expnts(8, 3);
  vcl_vector<vnl_real_npolynomial*> npoly_vector(3);  
  vcl_vector<vnl_vector<double>*> roots;
  
  unsigned int index;

  for (unsigned int i = 0; i < 2; i++)
  {
     for (unsigned int j = 0; j < 2; j++)
     {
        for (unsigned int k = 0; k < 2; k++)
        {
    index = i*4+j*2+k;
           expnts(index, 0) = 2-i-1;
           expnts(index, 1) = 2-j-1;             
           expnts(index, 2) = 2-k-1;             
 }   
     } 
  }   
  npoly_vector[0] = new vnl_real_npolynomial;  npoly_vector[0]->set(coeffs, expnts);
  npoly_vector[1] = new vnl_real_npolynomial;  npoly_vector[1]->set(coeffs, expnts);
  npoly_vector[2] = new vnl_real_npolynomial;  npoly_vector[2]->set(coeffs, expnts);

  for (unsigned int k = 0; k < 2; k++)
  {
     for (unsigned int j = 0; j < 2; j++)
     {
 for (unsigned int i = 0; i < 2; i++)
 {
    index = k*4+j*2+i;      
           coeffs[0] = sign[index][0]*sign[index][1]*sign[index][2];
           coeffs[1] = sign[index][0]*sign[index][1];
           coeffs[2] = sign[index][0]*sign[index][2];
           coeffs[3] = sign[index][0];
           coeffs[4] = sign[index][1]*sign[index][2];
           coeffs[5] = sign[index][1];
           coeffs[6] = sign[index][2];
    coeffs[7] = 1.0;

           x = this->m_node[index]->GetCoordinates()[0];
           y = this->m_node[index]->GetCoordinates()[1];
           z = this->m_node[index]->GetCoordinates()[2];

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
     if (localPt[0] > -1.0 && localPt[0] < 1.0 && localPt[1] > -1.0 && localPt[1] < 1.0 && localPt[2] > -1.0 && localPt[2] < 1.0)
     {
        isInside = true;
 break;
     } 
  }    
  return isInside;

}




/*
 * Draw the element on device context pDC.
 */
#ifdef FEM_BUILD_VISUALIZATION
void
Element3DC0LinearHexahedron
::Draw(CDC* pDC, Solution::ConstPointer sol) const 
{

  int x1=m_node[0]->GetCoordinates()[0]*DC_Scale;
  int y1=m_node[0]->GetCoordinates()[1]*DC_Scale;
  int z1=m_node[0]->GetCoordinates()[2]*DC_Scale;

  int x2=m_node[1]->GetCoordinates()[0]*DC_Scale;
  int y2=m_node[1]->GetCoordinates()[1]*DC_Scale;
  int z2=m_node[1]->GetCoordinates()[2]*DC_Scale;
  
  int x3=m_node[2]->GetCoordinates()[0]*DC_Scale;
  int y3=m_node[2]->GetCoordinates()[1]*DC_Scale;
  int z3=m_node[2]->GetCoordinates()[2]*DC_Scale;
  
  int x4=m_node[3]->GetCoordinates()[0]*DC_Scale;
  int y4=m_node[3]->GetCoordinates()[1]*DC_Scale;
  int z4=m_node[3]->GetCoordinates()[2]*DC_Scale;

  int x5=m_node[4]->GetCoordinates()[0]*DC_Scale;
  int y5=m_node[4]->GetCoordinates()[1]*DC_Scale;
  int z5=m_node[4]->GetCoordinates()[2]*DC_Scale;

  int x6=m_node[5]->GetCoordinates()[0]*DC_Scale;
  int y6=m_node[5]->GetCoordinates()[1]*DC_Scale;
  int z6=m_node[5]->GetCoordinates()[2]*DC_Scale;
  
  int x7=m_node[6]->GetCoordinates()[0]*DC_Scale;
  int y7=m_node[6]->GetCoordinates()[1]*DC_Scale;
  int z7=m_node[6]->GetCoordinates()[2]*DC_Scale;
  
  int x8=m_node[7]->GetCoordinates()[0]*DC_Scale;
  int y8=m_node[7]->GetCoordinates()[1]*DC_Scale;
  int z8=m_node[7]->GetCoordinates()[2]*DC_Scale;

  x1+=sol->GetSolutionValue(this->m_node[0]->GetDegreeOfFreedom(0))*DC_Scale;
  y1+=sol->GetSolutionValue(this->m_node[0]->GetDegreeOfFreedom(1))*DC_Scale;
  z1+=sol->GetSolutionValue(this->m_node[0]->GetDegreeOfFreedom(2))*DC_Scale;

  x2+=sol->GetSolutionValue(this->m_node[1]->GetDegreeOfFreedom(0))*DC_Scale;
  y2+=sol->GetSolutionValue(this->m_node[1]->GetDegreeOfFreedom(1))*DC_Scale;
  z2+=sol->GetSolutionValue(this->m_node[1]->GetDegreeOfFreedom(2))*DC_Scale;

  x3+=sol->GetSolutionValue(this->m_node[2]->GetDegreeOfFreedom(0))*DC_Scale;
  y3+=sol->GetSolutionValue(this->m_node[2]->GetDegreeOfFreedom(1))*DC_Scale;
  z3+=sol->GetSolutionValue(this->m_node[2]->GetDegreeOfFreedom(2))*DC_Scale;

  x4+=sol->GetSolutionValue(this->m_node[3]->GetDegreeOfFreedom(0))*DC_Scale;
  y4+=sol->GetSolutionValue(this->m_node[3]->GetDegreeOfFreedom(1))*DC_Scale;
  z4+=sol->GetSolutionValue(this->m_node[3]->GetDegreeOfFreedom(2))*DC_Scale;

  x5+=sol->GetSolutionValue(this->m_node[4]->GetDegreeOfFreedom(0))*DC_Scale;
  y5+=sol->GetSolutionValue(this->m_node[4]->GetDegreeOfFreedom(1))*DC_Scale;
  z5+=sol->GetSolutionValue(this->m_node[4]->GetDegreeOfFreedom(2))*DC_Scale;

  x6+=sol->GetSolutionValue(this->m_node[5]->GetDegreeOfFreedom(0))*DC_Scale;
  y6+=sol->GetSolutionValue(this->m_node[5]->GetDegreeOfFreedom(1))*DC_Scale;
  z6+=sol->GetSolutionValue(this->m_node[5]->GetDegreeOfFreedom(2))*DC_Scale;

  x7+=sol->GetSolutionValue(this->m_node[6]->GetDegreeOfFreedom(0))*DC_Scale;
  y7+=sol->GetSolutionValue(this->m_node[6]->GetDegreeOfFreedom(1))*DC_Scale;
  z7+=sol->GetSolutionValue(this->m_node[6]->GetDegreeOfFreedom(2))*DC_Scale;

  x8+=sol->GetSolutionValue(this->m_node[7]->GetDegreeOfFreedom(0))*DC_Scale;
  y8+=sol->GetSolutionValue(this->m_node[7]->GetDegreeOfFreedom(1))*DC_Scale;
  z8+=sol->GetSolutionValue(this->m_node[7]->GetDegreeOfFreedom(2))*DC_Scale;

  // FIXME: this isn't the correct drawing scheme
/*  pDC->MoveTo(x1,y1,z1);
  pDC->LineTo(x2,y2,z2);
  pDC->LineTo(x3,y3,z3);
  pDC->LineTo(x4,y4,z4);
  pDC->LineTo(x5,y5,z5);
  pDC->LineTo(x6,y6,z6);
  pDC->LineTo(x7,y7,z7);
  pDC->LineTo(x8,y8,z8);
  */

}
#endif




}} // end namespace itk::fem

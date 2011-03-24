/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMGenerateBSplineMesh.cxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:22:49 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itkFEMGenerateBSplineMesh.h"
#include "itkFEMElement2DBSplinePatch.h"
#include "itkFEMElement3DBSplinePatch.h"

#include "vcl_vector.h"

#include <math.h>

#define CLAMPED_BSPLINES


namespace itk {
namespace fem {
/*
 * Generate a rectangular mesh of quadrilateral elements
 */
 
void Generate2DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, 
                                      vnl_vector<double>& orig, vnl_vector<double>& size, 
                                      vnl_vector<double>& Nel, unsigned int BSplineOrder)
{
  // Check for correct number of dimensions
  if(orig.size() != Element2DBSplinePatch::NumberOfSpatialDimensions ||
     size.size() != Element2DBSplinePatch::NumberOfSpatialDimensions ||
     Nel.size()  != Element2DBSplinePatch::NumberOfSpatialDimensions)
  {
    throw FEMException(__FILE__, __LINE__, "GenerateBSplineMesh::Rectangular");
  }

  const unsigned int order = BSplineOrder;
  unsigned int m0 = static_cast<unsigned int>(Nel[0])+2*order-1;
  unsigned int m1 = static_cast<unsigned int>(Nel[1])+2*order-1;
  vnl_vector<double> knots0(m0), knots1(m1);

  #ifdef CLAMPED_BSPLINES
    for (unsigned int i = 0; i < order; i++)
    {
      knots0[i] = 0.0;  knots0[m0-i-1] = 1.0;  
      knots1[i] = 0.0;  knots1[m1-i-1] = 1.0;  
    }
    for (unsigned int i = order; i < m0-order; i++)  knots0[i] = static_cast<double>(i-order+1)/static_cast<double>(m0-2*order+1);
    for (unsigned int i = order; i < m1-order; i++)  knots1[i] = static_cast<double>(i-order+1)/static_cast<double>(m1-2*order+1);
  #else  
    for (unsigned int i = 0; i < m0; i++)    knots0[i] = static_cast<double>(i)/static_cast<double>(m0-1);
    for (unsigned int i = 0; i < m1; i++)    knots1[i] = static_cast<double>(i)/static_cast<double>(m1-1);
  #endif

  Generate2DRectilinearBSplineMesh(e0, S, orig, size, Nel, order, knots0, knots1);
}

void Generate2DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, 
                                      vnl_vector<double>& orig, vnl_vector<double>& size, 
                                      vnl_vector<double>& Nel, unsigned int BSplineOrder, 
                                      vnl_vector<double>& knots0, vnl_vector<double>& knots1)
{

  // Check for correct number of dimensions
  if(orig.size() != Element2DBSplinePatch::NumberOfSpatialDimensions ||
     size.size() != Element2DBSplinePatch::NumberOfSpatialDimensions ||
     Nel.size()  != Element2DBSplinePatch::NumberOfSpatialDimensions)
  {
    throw FEMException(__FILE__, __LINE__, "GenerateBSplineMesh::Rectangular");
  }

  // Clear existing elements, loads and nodes in Solver
  S.load.clear();
  S.el.clear();
  S.node.clear();

  // Number of elements in each dimension
  Nel[0]=static_cast<double>(floor(Nel[0]));
  Nel[1]=static_cast<double>(floor(Nel[1]));
  double m0 = static_cast<double>(knots0.size());
  double m1 = static_cast<double>(knots1.size());
  
  // The order and number of control points in each dimension are found based on the following relationships
  //   between the number of control points, the order, and the number of spans
  //
  //   number_of_spans = number_of_control_points - order + 1
  //   number_of_knots = number_of_control_points + order
  //   ->  number_of_control_points = 0.5*(Nel+number_of_knots-1)
  
  //Number of control points (nodes) in each dimension
  double Ncps0 = 0.5*(Nel[0]+m0-1.0);
  double Ncps1 = 0.5*(Nel[1]+m1-1.0);
  
  unsigned int order = static_cast<int>(m0-Ncps0); // = static_cast<int>(m1-Ncps1)
  if (order != BSplineOrder)
  {
     throw FEMException(__FILE__, __LINE__, "GenerateBSplineMesh::Rectangular");
  }

  // Create nodes
  Node::Pointer n;
  int gn=0; // number of node
  for(double j=0; j< Ncps1; j++)
  {
    for(double i=0; i< Ncps0; i++)  
    {
      #ifdef CLAMPED_BSPLINES
        n = new Node(orig[0]+i*size[0]/(Ncps0-1.0), orig[1]+j*size[1]/(Ncps1-1.0));
      #else
        n = new Node(orig[0]-(0.5*static_cast<double>(order-2)-i)*size[0]/(Ncps0-static_cast<double>(order-1)), 
              orig[1]-(0.5*static_cast<double>(order-2)-j)*size[1]/(Ncps1-static_cast<double>(order-1)));
      #endif 
      
      n->GN = gn;
      gn++;
      S.node.push_back(FEMP<Node>(n));
    }
  }

  unsigned int t0 = order*order;                                           // Nodes per element
  unsigned int t1 = order*(order-1);
  unsigned int t2 = (order-1)*(order-1);
  unsigned int f0 = (t0*(t0-1) >> 1);                                        
  unsigned int f1 = (t1*(t1-1) >> 1);                                        
  unsigned int f2 = (t2*(t2-1) >> 1);                                        
  
  unsigned int ne0 = f0;
  unsigned int ne1 = f0 - f1;
  unsigned int ne2 = f0 - f1 - (f1 - f2);

  long nze = 2*(gn + 2*(ne0 + ne1*(static_cast<int>(Nel[0]+Nel[1])-2) + ne2*(static_cast<int>(Nel[0])-1)*(static_cast<int>(Nel[1])-1)));

  S.SetMaximumNumberOfNonZeroElements(nze);

  vnl_vector<double> boundaries(4);
  vnl_matrix<double> shapeFunctions(2*order, order);

  BSplineShapeFunctions shape0;  shape0.SetKnots(knots0);
  BSplineShapeFunctions shape1;  shape1.SetKnots(knots1);
  
  // Create elements  
  gn=0; // global number of the element
  
  Element2DBSplinePatch::Pointer e1;
  for(unsigned int j=0; j<Nel[1]; j++)
  {
    boundaries[2] =  knots1[order-1+j];
    boundaries[3] =  knots1[order+j];
    for(unsigned int i=0; i<Nel[0]; i++)
    {
      boundaries[0] =  knots0[order-1+i];
      boundaries[1] =  knots0[order+i];
      
      e1 = dynamic_cast<Element2DBSplinePatch*>(e0->Clone());
      e1->SetLocalCoordinateBoundaries(boundaries);
      e1->SetBSplineOrder(order);
      e1->GN=gn;        
     
      for(unsigned int h=0; h<order; h++)
        for(unsigned int g=0; g<order; g++)
           e1->SetNode(g+h*order,S.node.Find((unsigned int) ((i+g)+Ncps0*(j+h))));

      for(unsigned int g=0; g<order; g++)
         shapeFunctions.set_row(g, (shape0.GenerateBSplineShapeFunction(order, i+g, i+order-1)).coefficients());
      for(unsigned int h=0; h<order; h++)
         shapeFunctions.set_row(h+order, (shape0.GenerateBSplineShapeFunction(order, j+h, j+order-1)).coefficients());  
  
      e1->SetBSplineShapeFunctions(shapeFunctions);      
      gn++;
      S.el.push_back(FEMP<Element>(e1));
    }
  }
}

void Generate2DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, 
                                      vnl_vector<double>& orig, vnl_vector<double>& size, 
                                      vnl_vector<double>& Nel, unsigned int BSplineOrder, 
                                      Image<char, 2>::Pointer maskImage)
{
  // Check for correct number of dimensions
  if(orig.size() != Element2DBSplinePatch::NumberOfSpatialDimensions ||
     size.size() != Element2DBSplinePatch::NumberOfSpatialDimensions ||
     Nel.size()  != Element2DBSplinePatch::NumberOfSpatialDimensions)
  {
    throw FEMException(__FILE__, __LINE__, "GenerateBSplineMesh::Rectangular");
  }

  const unsigned int order = BSplineOrder;
  unsigned int m0 = static_cast<unsigned int>(Nel[0])+2*order-1;
  unsigned int m1 = static_cast<unsigned int>(Nel[1])+2*order-1;
  vnl_vector<double> knots0(m0), knots1(m1);
  
  #ifdef CLAMPED_BSPLINES
    for (unsigned int i = 0; i < order; i++)
    {
      knots0[i] = 0.0;  knots0[m0-i-1] = 1.0;  
      knots1[i] = 0.0;  knots1[m1-i-1] = 1.0;  
    }

    for (unsigned int i = order; i < m0-order; i++)  knots0[i] = static_cast<double>(i-order+1)/static_cast<double>(m0-2*order+1);
    for (unsigned int i = order; i < m1-order; i++)  knots1[i] = static_cast<double>(i-order+1)/static_cast<double>(m1-2*order+1);
  #else  
    for (unsigned int i = 0; i < m0; i++)    knots0[i] = static_cast<double>(i)/static_cast<double>(m0-1);
    for (unsigned int i = 0; i < m1; i++)    knots1[i] = static_cast<double>(i)/static_cast<double>(m1-1);
  #endif

  Generate2DRectilinearBSplineMesh(e0, S, orig, size, Nel, order, knots0, knots1, maskImage);
}

void Generate2DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, 
                                      vnl_vector<double>& orig, vnl_vector<double>& size, 
                                      vnl_vector<double>& Nel, unsigned int BSplineOrder, 
                                      vnl_vector<double>& knots0, vnl_vector<double>& knots1, 
                                     Image<char, 2>::Pointer maskImage)
{
  // Check for correct number of dimensions
  if(orig.size() != Element2DBSplinePatch::NumberOfSpatialDimensions ||
     size.size() != Element2DBSplinePatch::NumberOfSpatialDimensions ||
     Nel.size()  != Element2DBSplinePatch::NumberOfSpatialDimensions)
  {
    throw FEMException(__FILE__, __LINE__, "GenerateBSplineMesh::Rectangular");
  }

  // Clear existing elements, loads and nodes in Solver
  S.load.clear();
  S.el.clear();
  S.node.clear();

  // Number of elements in each dimension
  Nel[0]=static_cast<double>(floor(Nel[0]));
  Nel[1]=static_cast<double>(floor(Nel[1]));
  double m0 = static_cast<double>(knots0.size());
  double m1 = static_cast<double>(knots1.size());
  
  // The order and number of control points in each dimension are found based on the following relationships
  //   between the number of control points, the order, and the number of spans
  //
  //   number_of_spans = number_of_control_points - order + 1
  //   number_of_knots = number_of_control_points + order
  //   ->  number_of_control_points = 0.5*(Nel+number_of_knots-1)
  
  //Number of control points (nodes) in each dimension
  double Ncps0 = 0.5*(Nel[0]+m0-1.0);
  double Ncps1 = 0.5*(Nel[1]+m1-1.0);
  
  unsigned int order = static_cast<int>(m0-Ncps0); // = static_cast<int>(m1-Ncps1)
  if (order != BSplineOrder)
  {
     throw FEMException(__FILE__, __LINE__, "GenerateBSplineMesh::Rectangular");
  }

  // Create nodes
  Node::Pointer n;
  int gn=0; // number of node
  for(double j=0; j< Ncps1; j++)
  {
    for(double i=0; i< Ncps0; i++)  
    {
      #ifdef CLAMPED_BSPLINES
        n = new Node(orig[0]+i*size[0]/(Ncps0-1.0), orig[1]+j*size[1]/(Ncps1-1.0));
      #else
        n = new Node(orig[0]-(0.5*static_cast<double>(order-2)-i)*size[0]/(Ncps0-static_cast<double>(order-1)), 
              orig[1]-(0.5*static_cast<double>(order-2)-j)*size[1]/(Ncps1-static_cast<double>(order-1)));
      #endif 

      n->GN=gn;
      gn++;
      S.node.push_back(FEMP<Node>(n));
    }
  }

  unsigned int t0 = order*order;                                           // Nodes per element
  unsigned int t1 = order*(order-1);
  unsigned int t2 = (order-1)*(order-1);
  unsigned int f0 = (t0*(t0-1) >> 1);                                        
  unsigned int f1 = (t1*(t1-1) >> 1);                                        
  unsigned int f2 = (t2*(t2-1) >> 1);                                        
  
  unsigned int ne0 = f0;
  unsigned int ne1 = f0 - f1;
  unsigned int ne2 = f0 - f1 - (f1 - f2);

  long nze = 2*(gn + 2*(ne0 + ne1*(static_cast<int>(Nel[0]+Nel[1])-2) + ne2*(static_cast<int>(Nel[0])-1)*(static_cast<int>(Nel[1])-1)));

  S.SetMaximumNumberOfNonZeroElements(nze);

  vnl_vector<double> boundaries(4);
  vnl_matrix<double> shapeFunctions(2*order, order);

  BSplineShapeFunctions shape0;  shape0.SetKnots(knots0);
  BSplineShapeFunctions shape1;  shape1.SetKnots(knots1);

  vnl_vector<double> Corner(2);
  bool OneCornerInside, OneCornerOutside;
  
  // Lay out the elements to form the mesh
  // Label each element as `inner,' `outer,' or `border'  
  gn=0; // global number of the element
  Element2DBSplinePatch::Pointer e1;
  for(unsigned int j=0; j<Nel[1]; j++)
  {
    boundaries[2] =  knots1[order-1+j];
    boundaries[3] =  knots1[order+j];
    for(unsigned int i=0; i<Nel[0]; i++)
    {
      boundaries[0] =  knots0[order-1+i];
      boundaries[1] =  knots0[order+i];
      
      e1 = dynamic_cast<Element2DBSplinePatch*>(e0->Clone());
      e1->SetLocalCoordinateBoundaries(boundaries);
      e1->SetExtendedBSplineStatus(false);
      e1->SetNumberOfNodes(order*order);
      e1->GN=gn;
     
      for(unsigned int h=0; h<order; h++)
        for(unsigned int g=0; g<order; g++)
           e1->SetNode(g+h*order,S.node.Find((unsigned int) ((i+g)+Ncps0*(j+h))));

      for(unsigned int g=0; g<order; g++)
         shapeFunctions.set_row(g, (shape0.GenerateBSplineShapeFunction(order, i+g, i+order-1)).coefficients());
      for(unsigned int h=0; h<order; h++)
         shapeFunctions.set_row(h+order, (shape0.GenerateBSplineShapeFunction(order, j+h, j+order-1)).coefficients());  
  
      e1->SetBSplineShapeFunctions(shapeFunctions); 
            
      OneCornerInside = OneCornerOutside = false;
      for(unsigned int c = 0; c < e1->GetNumberOfCorners(); c++)
      {
        Corner = e1->GetCornerCoordinates(c);
 Image<char, 2>::IndexType idx;
 idx[0] = static_cast<unsigned int>(Corner[0]+0.5);      
 idx[1] = static_cast<unsigned int>(Corner[1]+0.5);  

 if (maskImage->GetPixel(idx) == 0)   // Assume 'background' pixels are 0
   OneCornerOutside = true;
 else
   OneCornerInside = true;
      }   
      if (OneCornerInside && OneCornerOutside)
 e1->SetPartitionStatus(Element2DBSplinePatch::BORDER);
      else if (OneCornerInside && !OneCornerOutside)
 e1->SetPartitionStatus(Element2DBSplinePatch::INNER);
      else
 e1->SetPartitionStatus(Element2DBSplinePatch::OUTER);
 
      e1->SetExtendedBSplineStatus(true);  

      gn++;      
      S.el.push_back(FEMP<Element>(e1));
    }
  }  
  // Classify the B-splines as `inner,' `outer,' or `border'
  
  vnl_matrix<unsigned int> I(static_cast<unsigned int>(Ncps0)*static_cast<unsigned int>(Ncps1), 2);
  vnl_matrix<unsigned int> J(static_cast<unsigned int>(Ncps0)*static_cast<unsigned int>(Ncps1), 2);
  unsigned int countI = 0, countJ = 0;
  int tmpStatus;
   
  for(unsigned int j = 0; j < static_cast<unsigned int>(Nel[1])-order+1; j++)
  {
    for(unsigned int i = 0; i < static_cast<unsigned int>(Nel[0])-order+1; i++)
    {
      bool isJ = false;
      bool isI = false;

      for (unsigned int n = 0; n < order; n++) 
      {
 for (unsigned int m = 0; m < order; m++) 
 {
          e1 = dynamic_cast<Element2DBSplinePatch*>(S.el.Find((unsigned int) ((i+m)+static_cast<unsigned int>(Nel[0])*(j+n))));
          tmpStatus = static_cast<int>(e1->GetPartitionStatus());

   if (tmpStatus == -1)       // Inner element
          {
     isI = true;
     break;
   } 
   else if (tmpStatus == 0)  // Border element
   {
            isJ = true;
   }
        }
 if (isI) 
 {
   break;
 }
      }
      if (isI)
      {
 I(countI, 0) = i;  
 I(countI++, 1) = j;
      }
      else if (isJ)
      {
 J(countJ, 0) = i;  
 J(countJ++, 1) = j;
      }
    }
  }
  I = I.get_n_rows(0, countI) + order - 1;
  J = J.get_n_rows(0, countJ) + order - 1;  
  
  // Reset the nodes so that the Solver only includes the  
  // inner nodes in the K matrix.
  
  S.node.clear();
  for (unsigned int i = 0; i < I.rows(); i++) 
  {
    n = new Node(orig[0]+static_cast<double>(I(i, 0))*size[0]/(Ncps0-1.0), 
                 orig[1]+static_cast<double>(I(i, 1))*size[1]/(Ncps1-1.0));
    n->GN = i;
    S.node.push_back(FEMP<Node>(n));
  }

  // Determine the index matrices L, Ji, and Ij
  
  vnl_matrix<unsigned int> L(I.rows(), 2);
  vnl_vector<unsigned int> tmpI(2), tmpJ(2);
  unsigned countL = 0;
  bool flag, found;

  for (unsigned int i = 0; i < I.rows(); i++)
  {
    flag = true;
    for (unsigned int n = 0; n < order; n++)
    {
      for (unsigned int m = 0; m < order; m++)
      { 
 tmpI(0) = I(i, 0)+m;
 tmpI(1) = I(i, 1)+n;

        found = false;
        for (unsigned int ii = 0; ii < I.rows(); ii++)
 {         
   if (tmpI(0) == I(ii, 0) && tmpI(1) == I(ii, 1))
   {
     found = true;
     break;
          }
 }
 if (!found)
 {
   flag = false;
   break;
 }
      }
    } 
    if (flag == true)
    {
      L.set_row(countL++, I.get_row(i));
    }
  }   
  L = L.get_n_rows(0, countL);
  
  vnl_vector<int> Ij(J.rows());
  vnl_vector<unsigned int> Jv, Lv;
  vcl_vector<vnl_vector<unsigned int>*> Ji(I.rows());
  vnl_vector<unsigned int> countJi(I.rows());
  unsigned int mindist, dist;

  countJi = 0;
  for (unsigned int i = 0; i < I.rows(); i++)
  {
     Ji[i] = new vnl_vector<unsigned int>;  Ji[i]->set_size(I.rows());
  } 
  
  for (unsigned int j = 0; j < J.rows(); j++)       
  {
    mindist = INT_MAX;
    Ij(j) = -1;
    Jv = J.get_row(j);
    for (unsigned int l = 0; l < L.rows(); l++)
    {
      Lv = L.get_row(l);
      dist = (Jv-Lv).squared_magnitude();
      if (dist < mindist)
      {
        mindist = dist;
 Ij(j) = l;
      }
    }
    
    for (unsigned int i = 0; i < I.rows(); i++)
    {
      if ((I(i, 0) - L(Ij(j), 0)) >= 0 && (I(i, 0) - L(Ij(j), 0)) < order &&
          (I(i, 1) - L(Ij(j), 1)) >= 0 && (I(i, 1) - L(Ij(j), 1)) < order )
      {    
        Ji[i]->put(countJi[i]++, j);
      }
    }  
  }
  
  for (unsigned int i = 0; i < I.rows(); i++)
  {
    (*Ji[i]) = Ji[i]->extract(countJi[i]);
  }
  
  // Calculate the modified (extended B-spline) basis functions and 
  //  assign the nodes to the proper elements.
   
  vnl_vector<int> elementIndexI(I.rows());  
  vcl_vector<vnl_vector<unsigned int>*> elementIndexJi;
  vcl_vector<vnl_vector<double>*> elementJi_e;
  double e;   
 
  for(unsigned int j=0; j < Nel[1]; j++)
  {
    for(unsigned int i=0; i < Nel[0]; i++)
    {
      e1 = dynamic_cast<Element2DBSplinePatch*>(S.el.Find((unsigned int) (i+static_cast<unsigned int>(Nel[0])*j)));
      
      countJi = 0;
      countI = 0;
      elementIndexJi.resize(I.rows());
      elementJi_e.resize(I.rows());
      for (unsigned int ii = 0; ii < I.rows(); ii++)
      {

 tmpI(0) = I(ii, 0)-i;
 tmpI(1) = I(ii, 1)-j;
 if (tmpI(0) >= 0 && tmpI(0) < order && tmpI(1) >= 0 && tmpI(1) < order)
 {
   elementIndexI.put(countI, tmpI(1)*order+tmpI(0));
   elementIndexJi[countI] = new vnl_vector<unsigned int>;  elementIndexJi[countI]->set_size(Ji[ii]->size());
   elementJi_e[countI] = new vnl_vector<double>;         elementJi_e[countI]->set_size(Ji[ii]->size());
 
   for (unsigned int ji = 0; ji < Ji[ii]->size(); ji++)
   {
     tmpJ(0) = J(Ji[ii]->get(ji), 0)-i;
     tmpJ(1) = J(Ji[ii]->get(ji), 1)-j;
     if (tmpJ(0) >= 0 && tmpJ(0) < order && tmpJ(1) >= 0 && tmpJ(1) < order)
     {
        e = calculate_e(I.get_row(ii), J.get_row(Ji[ii]->get(ji)), L.get_row(Ij(Ji[ii]->get(ji))), order-1);
       if (fabs(e) > 0.0)
       {
         elementIndexJi[countI]->put(countJi[countI], tmpJ(1)*order+tmpJ(0));
         elementJi_e[countI]->put(countJi[countI]++, e);
              }
     }
   }
   countI++;
 } 
 else
 {
   bool isI = false;
   for (unsigned int ji = 0; ji < Ji[ii]->size(); ji++)
   {
     tmpJ(0) = J(Ji[ii]->get(ji), 0)-i;
     tmpJ(1) = J(Ji[ii]->get(ji), 1)-j;
     if (tmpJ(0) >= 0 && tmpJ(0) < order && tmpJ(1) >= 0 && tmpJ(1) < order)
     {
       isI = true;
       break;
     }
   }
   
   if (isI)
   {
           elementIndexI.put(countI, -1);
     elementIndexJi[countI] = new vnl_vector<unsigned int>;  elementIndexJi[countI]->set_size(Ji[ii]->size());
     elementJi_e[countI] = new vnl_vector<double>;         elementJi_e[countI]->set_size(Ji[ii]->size());
     
     for (unsigned int ji = 0; ji < Ji[ii]->size(); ji++)
     {
       tmpJ(0) = J(Ji[ii]->get(ji), 0)-i;
       tmpJ(1) = J(Ji[ii]->get(ji), 1)-j;
       if (tmpJ(0) >= 0 && tmpJ(0) < order && tmpJ(1) >= 0 && tmpJ(1) < order)
       {
   e = calculate_e(I.get_row(ii), J.get_row(Ji[ii]->get(ji)), L.get_row(Ij(Ji[ii]->get(ji))), order-1);
  if (fabs(e) > 0.0)  
  {
           elementIndexJi[countI]->put(countJi[countI], tmpJ(1)*order+tmpJ(0));
           elementJi_e[countI]->put(countJi[countI]++, e);
         }
       }
     }
     countI++;
   }  
 }
      }
      elementIndexJi.resize(countI);   
      elementJi_e.resize(countI);   
      for (unsigned int ii = 0; ii < countI; ii++)
      {
 (*elementIndexJi[ii]) = elementIndexJi[ii]->extract(countJi[ii]);
 (*elementJi_e[ii]) = elementJi_e[ii]->extract(countJi[ii]);
      }
      
      e1->SetBSplineShapeFunctionsIndexI(elementIndexI.extract(countI));
      e1->SetBSplineShapeFunctionsIndexJi(elementIndexJi);
      e1->SetBSplineShapeFunctionsJi_e(elementJi_e);
      e1->SetNumberOfNodes(countI);

      countI = 0;
      for (unsigned int ii = 0; ii < I.rows(); ii++)
      {
 tmpI(0) = I(ii, 0)-i;
 tmpI(1) = I(ii, 1)-j;
 if (tmpI(0) >= 0 && tmpI(0) < order && tmpI(1) >= 0 && tmpI(1) < order)
 {
   e1->SetNode(countI++, S.node.Find(ii));
 } 
 else
 {
   bool isI = false;
   for (unsigned int ji = 0; ji < Ji[ii]->size(); ji++)
   {
     tmpJ(0) = J(Ji[ii]->get(ji), 0)-i;
     tmpJ(1) = J(Ji[ii]->get(ji), 1)-j;
     if (tmpJ(0) >= 0 && tmpJ(0) < order && tmpJ(1) >= 0 && tmpJ(1) < order)
     {
       isI = true;
       break;
     }
   }
   
   if (isI)
   {
      e1->SetNode(countI++, S.node.Find(ii));
   }  
 }
      }   
    }  
  } 
}

/*
 * Generate a rectangular mesh of hexahedron elements
 */

void Generate3DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, 
                                      vnl_vector<double>& orig, vnl_vector<double>& size, 
                                      vnl_vector<double>& Nel, unsigned int BSplineOrder)
{
  // Check for correct number of dimensions
  if(orig.size() != Element3DBSplinePatch::NumberOfSpatialDimensions ||
     size.size() != Element3DBSplinePatch::NumberOfSpatialDimensions ||
     Nel.size()  != Element3DBSplinePatch::NumberOfSpatialDimensions)
  {
    throw FEMException(__FILE__, __LINE__, "GenerateBSplineMesh::Rectangular");
  }

  const unsigned int order = BSplineOrder;
  unsigned int m0 = static_cast<unsigned int>(Nel[0])+2*order-1;
  unsigned int m1 = static_cast<unsigned int>(Nel[1])+2*order-1;
  unsigned int m2 = static_cast<unsigned int>(Nel[2])+2*order-1;
  vnl_vector<double> knots0(m0), knots1(m1), knots2(m2);
  
  #ifdef CLAMPED_BSPLINES
    for (unsigned int i = 0; i < order; i++)
    {
      knots0[i] = 0.0;  knots0[m0-i-1] = 1.0;  
      knots1[i] = 0.0;  knots1[m1-i-1] = 1.0;  
      knots2[i] = 0.0;  knots2[m2-i-1] = 1.0;  
    }

    for (unsigned int i = order; i < m0-order; i++)  knots0[i] = static_cast<double>(i-order+1)/static_cast<double>(m0-2*order+1);
    for (unsigned int i = order; i < m1-order; i++)  knots1[i] = static_cast<double>(i-order+1)/static_cast<double>(m1-2*order+1);
    for (unsigned int i = order; i < m2-order; i++)  knots2[i] = static_cast<double>(i-order+1)/static_cast<double>(m2-2*order+1);
  #else  
    for (unsigned int i = 0; i < m0; i++)    knots0[i] = static_cast<double>(i)/static_cast<double>(m0-1);
    for (unsigned int i = 0; i < m1; i++)    knots1[i] = static_cast<double>(i)/static_cast<double>(m1-1);
    for (unsigned int i = 0; i < m2; i++)    knots2[i] = static_cast<double>(i)/static_cast<double>(m2-1);
  #endif
  
  Generate3DRectilinearBSplineMesh(e0, S, orig, size, Nel, order, knots0, knots1, knots2);
}
 
void Generate3DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S,  
                                      vnl_vector<double>& orig, vnl_vector<double>& size, 
                                      vnl_vector<double>& Nel, unsigned int BSplineOrder, 
                                      vnl_vector<double>& knots0, vnl_vector<double>& knots1, 
                                      vnl_vector<double>& knots2)
{

  // Check for correct number of dimensions
  if(orig.size() != Element3DBSplinePatch::NumberOfSpatialDimensions ||
     size.size() != Element3DBSplinePatch::NumberOfSpatialDimensions ||
     Nel.size()  != Element3DBSplinePatch::NumberOfSpatialDimensions)
  {
    throw FEMException(__FILE__, __LINE__, "GenerateBSplineMesh<Element3DBSplinePatch>::Rectangular");
  }

  // Number of nodes in each dimension
  Nel[0]=static_cast<double>(floor(Nel[0]));
  Nel[1]=static_cast<double>(floor(Nel[1]));
  Nel[2]=static_cast<double>(floor(Nel[2]));
  
  double m0 = static_cast<double>(knots0.size());
  double m1 = static_cast<double>(knots1.size());
  double m2 = static_cast<double>(knots2.size());
  
  // The order and number of control points in each dimension are found based on the following relationships
  //   between the number of control points, the order, and the number of spans
  //
  //   number_of_spans = number_of_control_points - order + 1
  //   number_of_knots = number_of_control_points + order
  //   ->  number_of_control_points = 0.5*(Nel+number_of_knots-1)
  
  //Number of control points (nodes) in each dimension
  double Ncps0 = 0.5*(Nel[0]+m0-1.0);
  double Ncps1 = 0.5*(Nel[1]+m1-1.0);
  double Ncps2 = 0.5*(Nel[2]+m2-1.0);
  
  unsigned int order = static_cast<int>(m0-Ncps0);
  if (order != BSplineOrder)
  {
    throw FEMException(__FILE__, __LINE__, "GenerateBSplineMesh::Rectangular");
  }  

  // Create nodes
  Node::Pointer n;
  int gn=0; // number of node
  for(double k=0; k<Ncps0; k++)
  {
    for(double j=0; j<Ncps1; j++)
    {
      for(double i=0; i<Ncps2; i++)
      {
 #ifdef CLAMPED_BSPLINES
          n = new Node(orig[0]+i*size[0]/(Ncps0-1.0), 
                orig[1]+j*size[1]/(Ncps1-1.0),
         orig[2]+k*size[2]/(Ncps2-1.0));
 #else
          n = new Node(orig[0]-(0.5*static_cast<double>(order-2)-i)*size[0]/(Ncps0-static_cast<double>(order-1)), 
                orig[1]-(0.5*static_cast<double>(order-2)-j)*size[1]/(Ncps1-static_cast<double>(order-1)),
         orig[2]-(0.5*static_cast<double>(order-2)-k)*size[2]/(Ncps2-static_cast<double>(order-1)));
 #endif 
 
        n->GN=gn;
        gn++;
        S.node.push_back(FEMP<Node>(n));
      }
    }
  }
  
  unsigned int t0 = order*order*order;                                           // Nodes per element
  unsigned int t1 = order*order*(order-1);
  unsigned int t2 = order*(order-1)*(order-1);
  unsigned int t3 = (order-1)*(order-1)*(order-1);
  unsigned int f0 = (t0*(t0-1) >> 1);                                        
  unsigned int f1 = (t1*(t1-1) >> 1);                                        
  unsigned int f2 = (t2*(t2-1) >> 1);                                        
  unsigned int f3 = (t3*(t3-1) >> 1);                                        
  
  unsigned int ne0 = f0;
  unsigned int ne1 = f0 - f1;
  unsigned int ne2 = f0 - f1 - (f1 - f2);
  unsigned int ne3 = f0 - f1 - (f1 - f2) - (f1 - f2 - f3);
  
  
  long nze = 3*(gn + 3*(ne0 + ne1*(static_cast<int>(Nel[0]+Nel[1]+Nel[2])-3) 
                            + ne2*((static_cast<int>(Nel[0])-1)*(static_cast<int>(Nel[1])-1)
             +(static_cast<int>(Nel[1])-1)*(static_cast<int>(Nel[2])-1)
      +(static_cast<int>(Nel[2])-1)*(static_cast<int>(Nel[0])-1))
                     + ne3*(static_cast<int>(Nel[0])-1)*(static_cast<int>(Nel[1])-1)*(static_cast<int>(Nel[2])-1)));
       
  S.SetMaximumNumberOfNonZeroElements(nze);
  
  vnl_vector<double> boundaries(6);
  vnl_matrix<double> shapeFunctions(3*order, order);

  // Create elements  
  BSplineShapeFunctions shape0;  shape0.SetKnots(knots0);
  BSplineShapeFunctions shape1;  shape1.SetKnots(knots1);
  BSplineShapeFunctions shape2;  shape2.SetKnots(knots2);  
  
  gn=0; // global number of the element
  Element3DBSplinePatch::Pointer e1;
  for(unsigned int k=0; k<Nel[2]; k++)
  {
    boundaries[4] =  knots2[order-1+k];
    boundaries[5] =  knots2[order+k];
    for(unsigned int j=0; j<Nel[1]; j++)
    {
      boundaries[2] =  knots1[order-1+j];
      boundaries[3] =  knots1[order+j];
      for(unsigned int i=0; i<Nel[0]; i++)
      {
        boundaries[0] =  knots0[order-1+i];
        boundaries[1] =  knots0[order+i];
 
        e1 = dynamic_cast<Element3DBSplinePatch*>(e0->Clone());
        e1->SetLocalCoordinateBoundaries(boundaries);
        e1->SetBSplineOrder(order);
        e1->GN=gn;
 
 for(unsigned int h=0; h<order; h++)
          for(unsigned int g=0; g<order; g++)
            for(unsigned int f=0; f<order; f++)   
              e1->SetNode(f+g*order+h*order*order,S.node.Find((unsigned int) ((i+f)+Ncps0*(j+g)+Ncps0*Ncps1*(k+h))));

        for(unsigned int f=0; f<order; f++)
           shapeFunctions.set_row(f, (shape0.GenerateBSplineShapeFunction(order, i+f, i+order-1)).coefficients());      
 for(unsigned int g=0; g<order; g++)
           shapeFunctions.set_row(g+order, (shape1.GenerateBSplineShapeFunction(order, j+g, j+order-1)).coefficients());   
 for(unsigned int h=0; h<order; h++)
           shapeFunctions.set_row(h+order*2, (shape2.GenerateBSplineShapeFunction(order, k+h, k+order-1)).coefficients());  

 e1->SetBSplineShapeFunctions(shapeFunctions);      
       
        gn++;
 S.el.push_back(FEMP<Element>(e1));
      }
    }
  }
}

void Generate3DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, 
                                      vnl_vector<double>& orig, vnl_vector<double>& size, 
                                      vnl_vector<double>& Nel, unsigned int BSplineOrder, 
                                      Image<char, 3>::Pointer maskImage)
{
  // Check for correct number of dimensions
  if(orig.size() != Element3DBSplinePatch::NumberOfSpatialDimensions ||
     size.size() != Element3DBSplinePatch::NumberOfSpatialDimensions ||
     Nel.size()  != Element3DBSplinePatch::NumberOfSpatialDimensions)
  {
    throw FEMException(__FILE__, __LINE__, "GenerateBSplineMesh::Rectangular");
  }

  const unsigned int order = BSplineOrder;
  unsigned int m0 = static_cast<unsigned int>(Nel[0])+2*order-1;
  unsigned int m1 = static_cast<unsigned int>(Nel[1])+2*order-1;
  unsigned int m2 = static_cast<unsigned int>(Nel[2])+2*order-1;
  vnl_vector<double> knots0(m0), knots1(m1), knots2(m2);
  
  #ifdef CLAMPED_BSPLINES
    for (unsigned int i = 0; i < order; i++)
    {
      knots0[i] = 0.0;  knots0[m0-i-1] = 1.0;  
      knots1[i] = 0.0;  knots1[m1-i-1] = 1.0;  
      knots2[i] = 0.0;  knots2[m2-i-1] = 1.0;  
    }

    for (unsigned int i = order; i < m0-order; i++)  knots0[i] = static_cast<double>(i-order+1)/static_cast<double>(m0-2*order+1);
    for (unsigned int i = order; i < m1-order; i++)  knots1[i] = static_cast<double>(i-order+1)/static_cast<double>(m1-2*order+1);
    for (unsigned int i = order; i < m2-order; i++)  knots2[i] = static_cast<double>(i-order+1)/static_cast<double>(m2-2*order+1);
  #else  
    for (unsigned int i = 0; i < m0; i++)    knots0[i] = static_cast<double>(i)/static_cast<double>(m0-1);
    for (unsigned int i = 0; i < m1; i++)    knots1[i] = static_cast<double>(i)/static_cast<double>(m1-1);
    for (unsigned int i = 0; i < m2; i++)    knots2[i] = static_cast<double>(i)/static_cast<double>(m2-1);
  #endif
  
  Generate3DRectilinearBSplineMesh(e0, S, orig, size, Nel, order, knots0, knots1, knots2, maskImage);
}

void Generate3DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, 
                                      vnl_vector<double>& orig, vnl_vector<double>& size, 
                                      vnl_vector<double>& Nel, unsigned int BSplineOrder, 
                                      vnl_vector<double>& knots0, vnl_vector<double>& knots1, 
                                      vnl_vector<double>& knots2, Image<char, 3>::Pointer maskImage)
{

  // Check for correct number of dimensions
  if(orig.size() != Element3DBSplinePatch::NumberOfSpatialDimensions ||
     size.size() != Element3DBSplinePatch::NumberOfSpatialDimensions ||
     Nel.size()  != Element3DBSplinePatch::NumberOfSpatialDimensions)
  {
    throw FEMException(__FILE__, __LINE__, "GenerateBSplineMesh<Element3DBSplinePatch>::Rectangular");
  }

  // Number of nodes in each dimension
  Nel[0]=static_cast<double>(floor(Nel[0]));
  Nel[1]=static_cast<double>(floor(Nel[1]));
  Nel[2]=static_cast<double>(floor(Nel[2]));
  
  double m0 = static_cast<double>(knots0.size());
  double m1 = static_cast<double>(knots1.size());
  double m2 = static_cast<double>(knots2.size());
  
  // The order and number of control points in each dimension are found based on the following relationships
  //   between the number of control points, the order, and the number of spans
  //
  //   number_of_spans = number_of_control_points - order + 1
  //   number_of_knots = number_of_control_points + order
  //   ->  number_of_control_points = 0.5*(Nel+number_of_knots-1)
  
  //Number of control points (nodes) in each dimension
  double Ncps0 = 0.5*(Nel[0]+m0-1.0);
  double Ncps1 = 0.5*(Nel[1]+m1-1.0);
  double Ncps2 = 0.5*(Nel[2]+m2-1.0);
  
  unsigned int order = static_cast<int>(m0-Ncps0);
  if (order != BSplineOrder)
  {
    throw FEMException(__FILE__, __LINE__, "GenerateBSplineMesh::Rectangular");
  }  

  // Create nodes
  Node::Pointer n;
  int gn=0; // number of node
  for(double k=0; k<Ncps0; k++)
  {
    for(double j=0; j<Ncps1; j++)
    {
      for(double i=0; i<Ncps2; i++)
      {
 #ifdef CLAMPED_BSPLINES
          n = new Node(orig[0]+i*size[0]/(Ncps0-1.0), 
                orig[1]+j*size[1]/(Ncps1-1.0),
         orig[2]+k*size[2]/(Ncps2-1.0));
 #else
          n = new Node(orig[0]-(0.5*static_cast<double>(order-2)-i)*size[0]/(Ncps0-static_cast<double>(order-1)), 
                orig[1]-(0.5*static_cast<double>(order-2)-j)*size[1]/(Ncps1-static_cast<double>(order-1)),
         orig[2]-(0.5*static_cast<double>(order-2)-k)*size[2]/(Ncps2-static_cast<double>(order-1)));
 #endif 

        n->GN=gn;
        gn++;
        S.node.push_back(FEMP<Node>(n));
      }
    }
  }
  
  unsigned int t0 = order*order*order;                                           // Nodes per element
  unsigned int t1 = order*order*(order-1);
  unsigned int t2 = order*(order-1)*(order-1);
  unsigned int t3 = (order-1)*(order-1)*(order-1);
  unsigned int f0 = (t0*(t0-1) >> 1);                                        
  unsigned int f1 = (t1*(t1-1) >> 1);                                        
  unsigned int f2 = (t2*(t2-1) >> 1);                                        
  unsigned int f3 = (t3*(t3-1) >> 1);                                        
  
  unsigned int ne0 = f0;
  unsigned int ne1 = f0 - f1;
  unsigned int ne2 = f0 - f1 - (f1 - f2);
  unsigned int ne3 = f0 - f1 - (f1 - f2) - (f1 - f2 - f3);
  
  
  long nze = 3*(gn + 2*(ne0 + ne1*(static_cast<int>(Nel[0]+Nel[1]+Nel[2])-3) 
                            + ne2*((static_cast<int>(Nel[0])-1)*(static_cast<int>(Nel[1])-1)
             +(static_cast<int>(Nel[1])-1)*(static_cast<int>(Nel[2])-1)
      +(static_cast<int>(Nel[2])-1)*(static_cast<int>(Nel[0])-1))
                     + ne3*(static_cast<int>(Nel[0])-1)*(static_cast<int>(Nel[1])-1)*(static_cast<int>(Nel[2])-1)));
       
  S.SetMaximumNumberOfNonZeroElements(nze);

  vnl_vector<double> boundaries(6);
  vnl_matrix<double> shapeFunctions(3*order, order);

  // Create elements  
  BSplineShapeFunctions shape0;  shape0.SetKnots(knots0);
  BSplineShapeFunctions shape1;  shape1.SetKnots(knots1);
  BSplineShapeFunctions shape2;  shape2.SetKnots(knots2);  
  
  vnl_vector<double> Corner(3);
  bool OneCornerInside, OneCornerOutside;
  
  gn=0; // global number of the element
  Element3DBSplinePatch::Pointer e1;
  for(unsigned int k=0; k<Nel[2]; k++)
  {
    boundaries[4] =  knots2[order-1+k];
    boundaries[5] =  knots2[order+k];
    for(unsigned int j=0; j<Nel[1]; j++)
    {
      boundaries[2] =  knots1[order-1+j];
      boundaries[3] =  knots1[order+j];
      for(unsigned int i=0; i<Nel[0]; i++)
      {
        boundaries[0] =  knots0[order-1+i];
        boundaries[1] =  knots0[order+i];
 
        e1 = dynamic_cast<Element3DBSplinePatch*>(e0->Clone());
        e1->SetLocalCoordinateBoundaries(boundaries);
    e1->SetExtendedBSplineStatus(false);
        e1->SetNumberOfNodes(order*order*order);
        e1->GN=gn;
 
 for(unsigned int h=0; h<order; h++)
          for(unsigned int g=0; g<order; g++)
            for(unsigned int f=0; f<order; f++)   
              e1->SetNode(f+g*order+h*order*order,S.node.Find((unsigned int) ((i+f)+Ncps0*(j+g)+Ncps0*Ncps1*(k+h))));

        for(unsigned int f=0; f<order; f++)
           shapeFunctions.set_row(f, (shape0.GenerateBSplineShapeFunction(order, i+f, i+order-1)).coefficients());      
 for(unsigned int g=0; g<order; g++)
           shapeFunctions.set_row(g+order, (shape1.GenerateBSplineShapeFunction(order, j+g, j+order-1)).coefficients());   
 for(unsigned int h=0; h<order; h++)
           shapeFunctions.set_row(h+order*2, (shape2.GenerateBSplineShapeFunction(order, k+h, k+order-1)).coefficients());  

 e1->SetBSplineShapeFunctions(shapeFunctions);      
 
 OneCornerInside = OneCornerOutside = false; 
 for(unsigned int c = 0; c < e1->GetNumberOfCorners(); c++)
 {
          Corner = e1->GetCornerCoordinates(c);
   Image<char, 3>::IndexType idx;
   idx[0] = static_cast<unsigned int>(Corner[0]+0.5);      
   idx[1] = static_cast<unsigned int>(Corner[1]+0.5);      
   idx[2] = static_cast<unsigned int>(Corner[2]+0.5);      
   if (maskImage->GetPixel(idx) == 0)    // Assume `background' pixels are 0
     OneCornerOutside = true;
   else
     OneCornerInside = true;
 }   
 if (OneCornerInside && OneCornerOutside)
    e1->SetPartitionStatus(Element3DBSplinePatch::BORDER);
 else if (OneCornerInside && !OneCornerOutside)
    e1->SetPartitionStatus(Element3DBSplinePatch::INNER);
 else
    e1->SetPartitionStatus(Element3DBSplinePatch::OUTER);
    
    e1->SetExtendedBSplineStatus(true);  
       
        gn++;
 S.el.push_back(FEMP<Element>(e1));
      }
    }
  }

  // Determine which B-splines are `inner,' `outer,' and `border'
  
  vnl_matrix<unsigned int> I(static_cast<unsigned int>(Ncps0)*static_cast<unsigned int>(Ncps1)*static_cast<unsigned int>(Ncps2), 3);
  vnl_matrix<unsigned int> J(static_cast<unsigned int>(Ncps0)*static_cast<unsigned int>(Ncps1)*static_cast<unsigned int>(Ncps2), 3);
  unsigned int countI = 0, countJ = 0;
  int tmpStatus;
   
  for(unsigned int k = 0; k < static_cast<unsigned int>(Nel[2])-order+1; k++)
  {
    for(unsigned int j = 0; j < static_cast<unsigned int>(Nel[1])-order+1; j++)
    {
      for(unsigned int i = 0; i < static_cast<unsigned int>(Nel[0])-order+1; i++)
      {
 bool isJ = false;
 bool isI = false;

 for (unsigned int n = 0; n < order; n++) 
 {
   for (unsigned int m = 0; m < order; m++) 
   {
     for (unsigned int l = 0; l < order; l++) 
     {
       e1 = dynamic_cast<Element3DBSplinePatch*>(S.el.Find((unsigned int) ((i+l)+static_cast<unsigned int>(Nel[0])
                                                 *(j+m)+static_cast<unsigned int>(Nel[0]*Nel[1])*(k+n))));
              tmpStatus = static_cast<int>(e1->GetPartitionStatus());

       if (tmpStatus == -1)  // Inner element
              {
  isI = true;
  break;
       } 
       else if (tmpStatus == 0)  // Border element
       {
         isJ = true;
       }
            }
     if (isI) 
     {
        break;
     }
   }
 }  
 if (isI)
 {
   I(countI, 0) = i;  
   I(countI, 1) = j;
   I(countI++, 2) = k;
 }
 else if (isJ)
 {
   J(countJ, 0) = i;  
   J(countJ, 1) = j;
   J(countJ++, 2) = k;
 }
      }
    }
  }  
  I = I.get_n_rows(0, countI) + order;
  J = J.get_n_rows(0, countJ) + order;
  

  // Re-create the nodes so that the Solver only includes the  
  // inner nodes in the K matrix.
  
  S.node.clear();
  for (unsigned int i = 0; i < I.rows(); i++) 
  {
    n = new Node(orig[0]+static_cast<double>(I(i, 0))*size[0]/(Ncps0-1.0), 
                 orig[1]+static_cast<double>(I(i, 1))*size[1]/(Ncps1-1.0),
   orig[2]+static_cast<double>(I(i, 2))*size[2]/(Ncps2-1.0));
    n->GN = i;
    S.node.push_back(FEMP<Node>(n));
  }

  
  // Determine the index matrices L, Ji, and Ij
  
  vnl_matrix<unsigned int> L(I.rows(), 3);
  vnl_vector<unsigned int> tmpI(3), tmpJ(3);
  unsigned countL = 0;
  bool flag, found;

  for (unsigned int i = 0; i < I.rows(); i++)
  {
    flag = true;
    for (unsigned int n = 0; n < order; n++)
    {
      for (unsigned int m = 0; m < order; m++)
      { 
 for (unsigned int l = 0; m < order; m++)
 { 
   tmpI(0) = I(i, 0)+l;
   tmpI(1) = I(i, 1)+m;
   tmpI(2) = I(i, 2)+n;

          found = false;
          for (unsigned int ii = 0; ii < I.rows(); ii++)
   {         
     if (tmpI(0) == I(ii, 0) && tmpI(1) == I(ii, 1) && tmpI(2) == I(ii, 2))
     {
       found = true;
       break;
            }
   }
   if (found)
   {
     flag = false;
     break;
   }
 }  
      }
    } 
    if (flag == true)
    {
      L.set_row(countL++, tmpI);
    }
  }   
  L = L.get_n_rows(0, countL);
 

  vnl_vector<int> Ij(J.rows());
  vnl_vector<unsigned int> Jv, Lv;
  vcl_vector<vnl_vector<unsigned int>*> Ji(I.rows());
  vnl_vector<unsigned int> countJi(I.rows());
  unsigned int mindist, dist;

  for (unsigned int i = 0; i < I.rows(); i++)
  {
     Ji[i] = new vnl_vector<unsigned int>;  Ji[i]->set_size(I.rows());
     Ji[i]->set_size(I.rows());
  } 
  
  for (unsigned int j = 0; j < J.rows(); j++)       
  {
    mindist = INT_MAX;
    Ij(j) = -1;
    Jv = J.get_row(j);
    for (unsigned int l = 0; l < L.rows(); l++)
    {
      Lv = L.get_row(l);
      dist = (Jv-Lv).squared_magnitude();
      if (dist < mindist)
      {
        mindist = dist;
 Ij(j) = l;
      }
    }
    
    countJi = 0;
    for (unsigned int i = 0; i < I.rows(); i++)
    {
      if ((I(i, 0) - L(Ij(j), 0)) >= 0 && (I(i, 0) - L(Ij(j), 0)) < order &&
          (I(i, 1) - L(Ij(j), 1)) >= 0 && (I(i, 1) - L(Ij(j), 1)) < order &&
          (I(i, 2) - L(Ij(j), 2)) >= 0 && (I(i, 2) - L(Ij(j), 2)) < order)
      {    
        Ji[i]->put(countJi[i]++, j);
      }
    }  
  }

  for (unsigned int i = 0; i < I.rows(); i++)
  {
    (*Ji[i]) = Ji[i]->extract(countJi[i]);
  }
  
  // Calculate the modified (extended B-spline) basis functions and 
  //  assign the nodes to the proper elements.
   
  vnl_vector<int> elementIndexI(I.rows());  
  vcl_vector<vnl_vector<unsigned int>*> elementIndexJi;
  vcl_vector<vnl_vector<double>*> elementJi_e;
  double e;   
 
  for(unsigned int k=0; k < Nel[2]; k++)
  {
    for(unsigned int j=0; j < Nel[1]; j++)
    {
      for(unsigned int i=0; i < Nel[0]; i++)
      {
 e1 = dynamic_cast<Element3DBSplinePatch*>(S.el.Find((unsigned int) (i+static_cast<unsigned int>(Nel[0])*j+static_cast<unsigned int>(Nel[0]*Nel[1])*k)));

 countJi = 0;
 countI = 0;
 elementIndexJi.resize(I.rows());
 elementJi_e.resize(I.rows());
 for (unsigned int ii = 0; ii < I.rows(); ii++)
 {

   tmpI(0) = I(ii, 0)-i;
   tmpI(1) = I(ii, 1)-j;
   tmpI(2) = I(ii, 2)-k;
   if (tmpI(0) >= 0 && tmpI(0) < order && tmpI(1) >= 0 && tmpI(1) < order && tmpI(2) >= 0 && tmpI(2) < order)
   {
      elementIndexI.put(countI, tmpI(2)*order*order+tmpI(1)*order+tmpI(0));
      elementIndexJi[countI] = new vnl_vector<unsigned int>;  elementIndexJi[countI]->set_size(Ji[ii]->size());
      elementJi_e[countI] = new vnl_vector<double>;         elementJi_e[countI]->set_size(Ji[ii]->size());

     for (unsigned int ji = 0; ji < Ji[ii]->size(); ji++)
     {
       tmpJ(0) = J(Ji[ii]->get(ji), 0)-i;
       tmpJ(1) = J(Ji[ii]->get(ji), 1)-j;
       tmpJ(2) = J(Ji[ii]->get(ji), 2)-k;
       if (tmpJ(0) >= 0 && tmpJ(0) < order && tmpJ(1) >= 0 && tmpJ(1) < order && tmpJ(2) >= 0 && tmpJ(2) < order)
       {
   e = calculate_e(I.get_row(ii), J.get_row(Ji[ii]->get(ji)), L.get_row(Ij(Ji[ii]->get(ji))), order-1);
  if (fabs(e) > 0.0)
  {
           elementIndexJi[countI]->put(countJi[countI], tmpJ(2)*order*order+tmpJ(1)*order+tmpJ(0));
           elementJi_e[countI]->put(countJi[countI]++, e);
          }
       }
     }
   } 
   else
   {
     bool isI = false;
     for (unsigned int ji = 0; ji < Ji[ii]->size(); ji++)
     {
       tmpJ(0) = J(Ji[ii]->get(ji), 0)-i;
       tmpJ(1) = J(Ji[ii]->get(ji), 1)-j;
       tmpJ(2) = J(Ji[ii]->get(ji), 2)-k;
       if (tmpJ(0) >= 0 && tmpJ(0) < order && tmpJ(1) >= 0 && tmpJ(1) < order && tmpJ(2) >= 0 && tmpJ(2) < order)
       {
  isI = true;
  break;
       }
     }

     if (isI)
     {
             elementIndexI.put(countI, -1);
       elementIndexJi[countI] = new vnl_vector<unsigned int>;  elementIndexJi[countI]->set_size(Ji[ii]->size());
       elementJi_e[countI] = new vnl_vector<double>;         elementJi_e[countI]->set_size(Ji[ii]->size());

       for (unsigned int ji = 0; ji < Ji[ii]->size(); ji++)
       {
  tmpJ(0) = J(Ji[ii]->get(ji), 0)-i;
  tmpJ(1) = J(Ji[ii]->get(ji), 1)-j;
  tmpJ(2) = J(Ji[ii]->get(ji), 2)-j;
  if (tmpJ(0) >= 0 && tmpJ(0) < order && tmpJ(1) >= 0 && tmpJ(1) < order && tmpJ(2) >= 0 && tmpJ(2) < order)
  {
     e = calculate_e(I.get_row(ii), J.get_row(Ji[ii]->get(ji)), L.get_row(Ij(Ji[ii]->get(ji))), order-1);
    if (fabs(e) > 0.0)
    {
             elementIndexJi[countI]->put(countJi[countI], tmpJ(2)*order*order+tmpJ(1)*order+tmpJ(0));
             elementJi_e[countI]->put(countJi[countI]++, e);
           }
  }
       }
     }  
   }
 }   

 elementIndexJi.resize(countI);   
 elementJi_e.resize(countI);   
 for (unsigned int ii = 0; ii < countI; ii++)
 {
   (*elementIndexJi[ii]) = elementIndexJi[ii]->extract(countJi[ii]);
   (*elementJi_e[ii]) = elementJi_e[ii]->extract(countJi[ii]);
 }
 
 e1->SetBSplineShapeFunctionsIndexI(elementIndexI.extract(countI));
 e1->SetBSplineShapeFunctionsIndexJi(elementIndexJi);
 e1->SetBSplineShapeFunctionsJi_e(elementJi_e);
 e1->SetNumberOfNodes(countI);

 countI = 0;
 for (unsigned int ii = 0; ii < I.rows(); ii++)
 {
   tmpI(0) = I(ii, 0)-i;
   tmpI(1) = I(ii, 1)-j;
   tmpI(2) = I(ii, 2)-k;
   if (tmpI(0) >= 0 && tmpI(0) < order && tmpI(1) >= 0 && tmpI(1) < order && tmpI(2) >= 0 && tmpI(2) < order)
   {
     e1->SetNode(countI++, S.node.Find(ii));
   } 
   else
   {
     bool isI = false;
     for (unsigned int ji = 0; ji < Ji[ii]->size(); ji++)
     {
       tmpJ(0) = J(Ji[ii]->get(ji), 0)-i;
       tmpJ(1) = J(Ji[ii]->get(ji), 1)-j;
       tmpJ(2) = J(Ji[ii]->get(ji), 2)-k;
       if (tmpJ(0) >= 0 && tmpJ(0) < order && tmpJ(1) >= 0 && tmpJ(1) < order && tmpJ(2) >= 0 && tmpJ(2) < order)
       {
  isI = true;
  break;
       }
     }

     if (isI)
     {
        e1->SetNode(countI++, S.node.Find(ii));
     }  
   }
 }   
      }  
    }    
  }  
}

double calculate_e(vnl_vector<unsigned int> i, vnl_vector<unsigned int> j, vnl_vector<unsigned int> l, unsigned int n)
{
  unsigned int m = i.size();
  double e = 1.0;

  for (unsigned int nu = 0; nu < m; nu++)
  {
    for (unsigned int mu = 0; mu <= n; mu++)
    {
      if ((l(nu)+mu) != i(nu))
      {
        e *= (static_cast<double>(j(nu))-static_cast<double>(l(nu))-static_cast<double>(mu));
 e /= (static_cast<double>(i(nu))-static_cast<double>(l(nu))-static_cast<double>(mu));
      }
    }
  }    

  return e;
}




/** Dummy functions that shouldn't be executed*/

void Generate2DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, vnl_vector<double>& orig, vnl_vector<double>& size, vnl_vector<double>& Nel, 
                                      unsigned int BSplineOrder, Image<char, 3>::Pointer maskImage)
{
  throw FEMException(__FILE__, __LINE__, "GenerateBSplineMesh::BSpline (mismatched dimensions)");
}

void Generate2DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, vnl_vector<double>& orig, vnl_vector<double>& size, vnl_vector<double>& Nel, 
                                      unsigned int BSplineOrder, vnl_vector<double>& knots0, vnl_vector<double>& knots1, Image<char, 3>::Pointer maskImage)
{
  throw FEMException(__FILE__, __LINE__, "GenerateBSplineMesh::BSpline (mismatched dimensions)");
}

void Generate3DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, vnl_vector<double>& orig, vnl_vector<double>& size, vnl_vector<double>& Nel, 
                                      unsigned int BSplineOrder, Image<char, 2>::Pointer maskImage)
{
  throw FEMException(__FILE__, __LINE__, "GenerateBSplineMesh::BSpline (mismatched dimensions)");
}

void Generate3DRectilinearBSplineMesh(itk::fem::Element::ConstPointer e0, Solver& S, vnl_vector<double>& orig, vnl_vector<double>& size, vnl_vector<double>& Nel, 
                                      unsigned int BSplineOrder, vnl_vector<double>& knots0, vnl_vector<double>& knots1, vnl_vector<double>& knots2, Image<char, 2>::Pointer maskImage)
{
  throw FEMException(__FILE__, __LINE__, "GenerateBSplineMesh::BSpline (mismatched dimensions)");
}





}} // end namespace itk::fem

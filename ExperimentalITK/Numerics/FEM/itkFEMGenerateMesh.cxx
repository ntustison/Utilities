/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMGenerateMesh.cxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:22:49 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itkFEMGenerateMesh.h"
#include "itkFEMElement2DC0LinearQuadrilateral.h"
#include "itkFEMElement3DC0LinearHexahedron.h"
#include <math.h>

namespace itk {
namespace fem {



/*
 * Generate a rectangular mesh of quadrilateral elements
 */
void Generate2DRectilinearMesh(itk::fem::Element::ConstPointer e0, Solver& S, vnl_vector<double>& orig, vnl_vector<double>& size, vnl_vector<double>& Nel)
{

  // Check for correct number of dimensions
  if(orig.size() != Element2DC0LinearQuadrilateral::NumberOfSpatialDimensions ||
     size.size() != Element2DC0LinearQuadrilateral::NumberOfSpatialDimensions ||
     Nel.size()  != Element2DC0LinearQuadrilateral::NumberOfSpatialDimensions)
  {
    throw FEMException(__FILE__, __LINE__, "GenerateMesh<Element2DC0LinearQuadrilateral>::Rectangular");
  }

  // Clear existing elements, loads and nodes in Solver
  S.load.clear();
  S.el.clear();
  S.node.clear();

  // Number of nodes in each dimension
  Nel[0]=floor(Nel[0]);
  Nel[1]=floor(Nel[1]);
  double Ni=static_cast<double>(Nel[0]);
  double Nj=static_cast<double>(Nel[1]);

  // Create nodes
  Node::Pointer n;
  int gn=0; // number of node
  for(double j=0; j<=Nj; j++)
  {
    for(double i=0; i<=Ni; i++)  
    {
      n=new Node(orig[0]+i*size[0]/Nel[0], orig[1]+j*size[1]/Nel[1]);
      n->GN=gn;
      gn++;
      S.node.push_back(FEMP<Node>(n));
    }
  }

  unsigned int order = 2;
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

  // Create elements  
  gn=0; // global number of the element
  
  Element2DC0LinearQuadrilateral::Pointer e;
  for(unsigned int j=0; j<Nj; j++)
  {
    for(unsigned int i=0; i<Ni; i++)
    {
      e=dynamic_cast<Element2DC0LinearQuadrilateral*>(e0->Clone());
      e->SetNode(0,S.node.Find((unsigned int) (i+  (Ni+1)*j)     ));
      e->SetNode(1,S.node.Find((unsigned int) (i+1+(Ni+1)*j)     ));
      e->SetNode(2,S.node.Find((unsigned int) (i+1+(Ni+1)*(j+1)) ));
      e->SetNode(3,S.node.Find((unsigned int) (i+  (Ni+1)*(j+1)) ));      
      e->GN=gn;
      
      gn++;
      S.el.push_back(FEMP<Element>(e));
    }
  }

}


/*
 * Generate a rectangular mesh of hexahedron elements
 */
void Generate3DRectilinearMesh
(itk::fem::Element::ConstPointer e0, Solver& S, vnl_vector<double>& orig, 
 vnl_vector<double>& size, vnl_vector<double>& Nel)
{

  // Check for correct number of dimensions
  if(orig.size() != Element3DC0LinearHexahedron::NumberOfSpatialDimensions ||
     size.size() != Element3DC0LinearHexahedron::NumberOfSpatialDimensions ||
     Nel.size()  != Element3DC0LinearHexahedron::NumberOfSpatialDimensions)
  {
    throw FEMException(__FILE__, __LINE__, "GenerateMesh<Element2DC0LinearQuadrilateral>::Rectangular");
  }

  // Number of nodes in each dimension
  Nel[0]=floor(Nel[0]);
  Nel[1]=floor(Nel[1]);
  Nel[2]=floor(Nel[2]);
  double Ni=static_cast<double>(Nel[0]);
  double Nj=static_cast<double>(Nel[1]);
  double Nk=static_cast<double>(Nel[2]);

  // Create nodes
  Node::Pointer n;
  int gn=0; // number of node
  for(double k=0; k<=Nk; k++)
  {
    for(double j=0; j<=Nj; j++)
    {
      for(double i=0; i<=Ni; i++)
      {
        double xx,yy,zz;
        xx=orig[0]+i*size[0]/Nel[0]; 
        yy=orig[1]+j*size[1]/Nel[1];
        zz=orig[2]+k*size[2]/Nel[2];
//std::cout << " xx " << xx << " yy " << yy << " zz " << zz << std::endl;
        n=new Node(xx,yy,zz);
        n->GN=gn;
        gn++;
        S.node.push_back(FEMP<Node>(n));
      }
    }
  }

  unsigned int order = 2;
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
       
//  S.SetMaximumNumberOfNonZeroElements(nze);
  

  // Create elements  
  gn=0; // global number of the element
  itk::fem::Element3DC0LinearHexahedron::Pointer e;
  for(unsigned int k=0; k<Nk; k++)
  {
    for(unsigned int j=0; j<Nj; j++)
    {
      for(unsigned int i=0; i<Ni; i++)
      {
        e=dynamic_cast<Element3DC0LinearHexahedron*>(e0->Clone());
        e->SetNode(0,S.node.Find((unsigned int) (i+  (Ni+1)*(j  +(Nj+1)*k) )));
        e->SetNode(1,S.node.Find((unsigned int) (i+1+(Ni+1)*(j  +(Nj+1)*k) )));
        e->SetNode(2,S.node.Find((unsigned int) (i+1+(Ni+1)*(j+1+(Nj+1)*k) )));
        e->SetNode(3,S.node.Find((unsigned int) (i+  (Ni+1)*(j+1+(Nj+1)*k) )));
        e->SetNode(4,S.node.Find((unsigned int) (i+  (Ni+1)*(j  +(Nj+1)*(k+1)) )));
        e->SetNode(5,S.node.Find((unsigned int) (i+1+(Ni+1)*(j  +(Nj+1)*(k+1)) )));
        e->SetNode(6,S.node.Find((unsigned int) (i+1+(Ni+1)*(j+1+(Nj+1)*(k+1)) )));
        e->SetNode(7,S.node.Find((unsigned int) (i+  (Ni+1)*(j+1+(Nj+1)*(k+1)) )));

        e->GN=gn;
        gn++;
        S.el.push_back(FEMP<Element>(e));
      }
    }
  }
}


}} // end namespace itk::fem

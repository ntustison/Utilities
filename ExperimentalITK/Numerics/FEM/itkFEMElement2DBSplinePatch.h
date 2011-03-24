/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMElement2DBSplinePatch.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:22:47 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkFEMElement2DBSplinePatch_h
#define __itkFEMElement2DBSplinePatch_h

#include "itkFEMElementBase.h"
#include "itkFEMBSplineShapeFunctions.h"

#include "vcl_vector.h"

namespace itk {
namespace fem {

/**
 * \class Element2DBSplinePatch
 * \brief C^(BSplineOrder-2) continuous finite element in 2D space.
 */
class Element2DBSplinePatch : public Element
{
typedef Element TemplatedParentClass;
FEM_ABSTRACT_CLASS( Element2DBSplinePatch, TemplatedParentClass )

  /*************************************************************************/
  //  The first part of this class contains information related to FEM-
  //  based B-spline image registration when using the entire 2-D image.
  //  Most of the functions are adapted to also employ *extended* B-splines.
  /*************************************************************************/
public:

  Element2DBSplinePatch() : m_BSplineOrder(3), m_NumberOfNodes(9)
   { 
     this->GenerateShapeFunctions(); 
     m_localCoordinateBoundaries.set_size(4);  
     m_localCoordinateBoundaries[0] = 0.0;
     m_localCoordinateBoundaries[1] = 1.0;
     m_localCoordinateBoundaries[2] = 0.0;
     m_localCoordinateBoundaries[3] = 1.0;

     m_ExtendedBSplineStatus = false;  
   }; 

  /**
   * Number of dimensions of space in which element can exist.
   */
  enum { NumberOfSpatialDimensions = 2 };

  typedef BSplineShapeFunctions::PolynomialType PolynomialType;

  /*
   * Methods related to numeric integration
   */
  enum { DefaultIntegrationOrder = 2 };
  virtual void GetIntegrationPointAndWeight(unsigned int i, VectorType& pt, Float& w, unsigned int order) const;
  virtual unsigned int GetNumberOfIntegrationPoints(unsigned int order) const;

  /*
   * Methods related to the geometry of an element
   */
  virtual VectorType ShapeFunctions( const VectorType& pt ) const;
  virtual void ShapeFunctionDerivatives( const VectorType& pt, MatrixType& shapeD ) const;
  virtual bool GetLocalFromGlobalCoordinates( const VectorType& globalPt, VectorType& localPt ) const;
  void SetLocalCoordinateBoundaries(VectorType v) { m_localCoordinateBoundaries = v; };
  void SetNumberOfNodes(unsigned int n) { m_NumberOfNodes = n; m_node.resize(n); };  
  VectorType GetLocalCoordinateBoundaries(void) { return m_localCoordinateBoundaries; };
  virtual unsigned int GetNumberOfNodes( void ) const { return m_NumberOfNodes; };  
  virtual unsigned int GetNumberOfCorners( void ) const { return 4; }  
  virtual NodeIDType GetNode(unsigned int n) const { return ((n >= m_NumberOfNodes) ? 0 : m_node[n]); };
  virtual void SetNode(unsigned int n, NodeIDType node) { if (n < m_NumberOfNodes) m_node[n] = node; };
  virtual const VectorType& GetNodeCoordinates( unsigned int n ) const  { return m_node[n]->GetCoordinates(); };
  virtual const VectorType GetCornerCoordinates( unsigned int n ) const;
  virtual unsigned int GetNumberOfSpatialDimensions() const { return NumberOfSpatialDimensions; };

  /**
   * Added functionality to handle the B-spline aspect of this element
   */
  void SetBSplineOrder(unsigned int order) 
    { 
      m_BSplineOrder = (order >= 2) ? order : 2; 
      this->SetNumberOfNodes(m_BSplineOrder*m_BSplineOrder);
      this->GenerateShapeFunctions();
    };
  unsigned int GetBSplineOrder(void) { return m_BSplineOrder; };
  void SetBSplineShapeFunctions(MatrixType M) { m_BSplineShapeFunctions = M; };
  MatrixType GetBSplineShapeFunctions(void) { return m_BSplineShapeFunctions; };
  
  /**
   * Draw the element on the specified device context
   */
  #ifdef FEM_BUILD_VISUALIZATION
  void Draw(CDC* pDC, Solution::ConstPointer sol) const; 
  #endif
private: 

  /**
   * Generate the shape functions if they weren't produced from the mesh generation.
   */
  void GenerateShapeFunctions();

  /**
   * BSplineOrder - the order of the BSpline object (degree = BSplineOrder-1)
   */
  unsigned int m_BSplineOrder;

  /**
   * BSplineShapeFunctions - 2-D Array of coefficient vectors for the B-spline polynomial 
   *                              shape functions.  Since the m_BSplineOrder^(dimension) element 
   *                              shape functions are separable, we only store m_BSplineOrder*(dimension) 
   *                              functions in each parametric direction.
   */
  MatrixType m_BSplineShapeFunctions;

  /**
   * m_NumberOfNodes - equal to m_BSplineOrder^(dimension) except in the case of 
   *                   extended b-splines
   */
  unsigned int m_NumberOfNodes; 
  
  /**
   * Array of pointers to point objects that define the element
   */
  vcl_vector<NodeIDType> m_node;
  
  /**
   * localCoordinateBoundaries - The local coordinates for the elements are typically in the range
   *                               [-1, 1].  However, the local coordinates of the B-spline finite 
   *                               elements are defined by the corresponding knot elements.  This variable
   *                               contains those limits (length = 4 --- [u0, u1, v0, v1]).
   */
  VectorType m_localCoordinateBoundaries;  


  /*************************************************************************/
  //  The second part of this class contains information related to FEM-
  //  based *extended* B-spline image registration when using an image mask
  //  region.
  /*************************************************************************/

public:
  /**
   * ExtendedBSplineStatus - 
   */     
  void SetExtendedBSplineStatus(bool b) { m_ExtendedBSplineStatus = b; };
  bool GetExtendedBSplineStatus(void) { return m_ExtendedBSplineStatus; };  
  
  /**
   * PartitionStatus - the status of an element ('inner', 'border', 'outer').  
   *             If the entire region of image is used all elements are 
   *             classified as 'inner.'
   */     
  enum PartitionStatus { INNER = -1, BORDER, OUTER };
  void SetPartitionStatus(PartitionStatus p) { m_PartitionStatus = p; };
  PartitionStatus GetPartitionStatus(void) { return m_PartitionStatus; };

  void SetBSplineShapeFunctionsIndexI(IntVectorType I) { m_BSplineShapeFunctionsIndexI = I; }
  IntVectorType GetBSplineShapeFunctionsIndexI(void) { return m_BSplineShapeFunctionsIndexI; };
  void SetBSplineShapeFunctionsIndexJi(vcl_vector<UIntVectorType*> Ji) { m_BSplineShapeFunctionsIndexJi = Ji; };  
  vcl_vector<UIntVectorType*> GetBSplineShapeFunctionsIndexJi(void) { return m_BSplineShapeFunctionsIndexJi; }; 
  void SetBSplineShapeFunctionsJi_e(vcl_vector<VectorType*> Ji_e) { m_BSplineShapeFunctionsJi_e = Ji_e;};  
  vcl_vector<VectorType*> GetBSplineShapeFunctionsJi_e(void) { return m_BSplineShapeFunctionsJi_e; };
 
private:  
   
  PartitionStatus m_PartitionStatus;     
  bool m_ExtendedBSplineStatus;

  IntVectorType m_BSplineShapeFunctionsIndexI;
  vcl_vector<UIntVectorType*> m_BSplineShapeFunctionsIndexJi;
  vcl_vector<VectorType*> m_BSplineShapeFunctionsJi_e;

};



}} // end namespace itk::fem

#endif  // #ifndef __itkFEMElement2DBSplinePatch_h

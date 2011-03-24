/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSparseMatrixEigenAnalysis.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:20:04 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkSparseMatrixEigenAnalysis_h
#define __itkSparseMatrixEigenAnalysis_h

#include "itkObject.h"

#include "itkMacro.h"

#include "vnl/vnl_sparse_matrix.h"
#include "vnl/vnl_vector.h"

#include <vector>

namespace itk
{
  
/** \class SparseMatrixEigenAnalysis
 */  

template < class TRealType = double, 
           class TSparseMatrix = vnl_sparse_matrix<TRealType>, 
           class TVector = vnl_vector<TRealType> >
class ITK_EXPORT SparseMatrixEigenAnalysis : public Object
{
public:
  /** Standard "Self" typedef. */
  typedef SparseMatrixEigenAnalysis                        Self;
  typedef Object                                           Superclass;
  typedef SmartPointer<Self>                               Pointer;
  typedef SmartPointer<const Self>                         ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro( SparseMatrixEigenAnalysis, Object );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );
  
  typedef TRealType                                        RealType;
  typedef TSparseMatrix                                    SparseMatrixType;
  typedef TVector                                          VectorType;
  typedef std::vector<VectorType>                          VectorContainerType;

  itkSetObjectMacro( Matrix1, SparseMatrixType );
//  itkGetObjectMacro( Matrix1, SparseMatrixType );

  itkSetObjectMacro( Matrix2, SparseMatrixType );
//  itkGetObjectMacro( Matrix2, SparseMatrixType );

  itkSetMacro( IsSymmetric, bool );
  itkGetConstMacro( IsSymmetric, bool );
  itkBooleanMacro( IsSymmetric );

  itkSetMacro( SolveForSmallestEigenValues, bool );
  itkGetConstMacro( SolveForSmallestEigenValues, bool );
  itkBooleanMacro( SolveForSmallestEigenValues );

  itkSetMacro( Tolerance, RealType );
  itkGetConstMacro( Tolerance, RealType );

  itkSetMacro( MaximumNumberOfIterations, unsigned int );
  itkGetConstMacro( MaximumNumberOfIterations, unsigned int );

  itkSetMacro( NumberOfEigenPairs, unsigned int );
  itkGetConstMacro( NumberOfEigenPairs, unsigned int );

  itkSetMacro( NumberOfLanczosVectors, unsigned int );
  itkGetConstMacro( NumberOfLanczosVectors, unsigned int );

  VectorType GetEigenVector( unsigned int n )
    {
    if ( n < this->m_NumberOfEigenPairs )
      {
      return this->m_EigenVectors[n];
      }  
    };

  RealType GetEigenValue( unsigned int n )
    {
    if ( n < this->m_NumberOfEigenPairs )
      {
      return this->m_EigenValues[n];
      }  
    };

  /** dummy method that calls the GenerateData method to 
   * produce the eigen vectors and values. */
  void Update()
  { this->GenerateData(); }

protected:
  SparseMatrixEigenAnalysis();
  virtual ~SparseMatrixEigenAnalysis();
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Produces the eigen vectors and values. */
  void GenerateData();

private:
  SparseMatrixEigenAnalysis(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  void CalculateEigenPairsFloat();
  void CalculateEigenPairsDouble();


  VectorType                                                m_EigenValues;
  VectorContainerType                                       m_EigenVectors;

  bool                                                      m_IsSymmetric;
//bool                                                      m_IsReal;
  bool                                                      m_SolveForSmallestEigenValues;
  
  RealType                                                  m_Tolerance;
  unsigned int                                              m_MaximumNumberOfIterations;
  unsigned int                                              m_NumberOfEigenPairs;
  unsigned int                                              m_NumberOfLanczosVectors;

  SparseMatrixType                                          *m_Matrix1;
  SparseMatrixType                                          *m_Matrix2;

};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSparseMatrixEigenAnalysis.txx"
#endif

#endif


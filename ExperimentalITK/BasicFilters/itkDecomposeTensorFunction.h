/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDecomposeTensorFunction.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:51 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDecomposeTensorFunction_h
#define __itkDecomposeTensorFunction_h

#include "itkProcessObject.h"

#include "itkVariableSizeMatrix.h"

namespace itk
{
/** \class DecomposeTensorFunction
 * \brief Wrapper for vnl_matrix decomposition routines.
 *
 * Several matrix decomposition routines are available in vnl.  Using th
 * these routines, the following decompositions are made available:
 *
 *   1. eigen-decomposition
 *   2. symmetric eigen-decomposition
 *   3. QR decomposition
 *   4. SVD decomposition
 *   5. Left and right polar decomposition
 *   6. Cholesky decomposition
 *
 * This class also provides an interface to evaluating the determinant
 * from an NxN matrix using vnl routines.
 *
 * \author Nicholas J. Tustison
 *
 */
template <typename TInput,
          typename TRealType = float,
          typename TOutput = itk::VariableSizeMatrix<TRealType> >
class ITK_EXPORT DecomposeTensorFunction : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef DecomposeTensorFunction                         Self;
  typedef ProcessObject                                   Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( DecomposeTensorFunction, ProcessObject );

  /** Extract some information from the input types.  Dimensionality
   * of the two inputs is assumed to be the same. */
  typedef TInput                                          InputMatrixType;
  typedef TOutput                                         OutputMatrixType;

  /** Define the data type and the vector of data type used in calculations. */
  typedef TRealType                                       RealType;

  // Wrappers for vnl routines
  void EvaluateEigenDecomposition( InputMatrixType&,
    OutputMatrixType&, OutputMatrixType& );

  void EvaluateSymmetricEigenDecomposition( InputMatrixType&,
    OutputMatrixType&, OutputMatrixType& );

  void EvaluateQRDecomposition( InputMatrixType&,
    OutputMatrixType&, OutputMatrixType& );

  void EvaluateSVDDecomposition( InputMatrixType&,
    OutputMatrixType&, OutputMatrixType&, OutputMatrixType& );

  void EvaluateSVDEconomyDecomposition( InputMatrixType&,
    OutputMatrixType&, OutputMatrixType& );

  void EvaluateLeftPolarDecomposition( InputMatrixType&,
    OutputMatrixType&, OutputMatrixType& );

  void EvaluateRightPolarDecomposition( InputMatrixType&,
    OutputMatrixType&, OutputMatrixType& );

  void EvaluateCholeskyDecomposition( InputMatrixType&,
    OutputMatrixType& );

  RealType EvaluateDeterminant( InputMatrixType& );

  DecomposeTensorFunction();
  virtual ~DecomposeTensorFunction() {}

protected:

  void PrintSelf ( std::ostream& os, Indent indent ) const;

private:

  DecomposeTensorFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDecomposeTensorFunction.txx"
#endif

#endif

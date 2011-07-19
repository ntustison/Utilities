/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGaussianClassDensityFunction.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGaussianClassDensityFunction_h
#define __itkGaussianClassDensityFunction_h

#include "itkArray.h"
#include "itkVariableSizeMatrix.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include "vnl/algo/vnl_determinant.h"
#include "vnl/vnl_math.h"

#include "itkMatrix.h"
#include "itkDensityFunction.h"

namespace itk{
namespace Statistics{

/** \class GaussianClassDensityFunction
 * \brief GaussianClassDensityFunction class represents Gaussian
 *  Density Function.
 *
 * This class keeps parameter to define Gaussian Density Function  and has
 * method to return the probability density of an instance (pattern).
 * If the all element of the covariance matrix is zero the "usual" density
 * calculations ignored. if the measurement vector to be evaluated is equal to
 * the mean, then the Evaluate method will return maximum value of
 * RealType and return 0 for others.
 *
 */

template<class TMeasurementVector>
class ITK_EXPORT GaussianClassDensityFunction :
    public DensityFunction<TMeasurementVector>
{
public:
  /** Standard class typedefs */
  typedef GaussianClassDensityFunction                    Self;
  typedef DensityFunction<TMeasurementVector>             Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;

  /** Strandard macros */
  itkTypeMacro( GaussianClassDensityFunction, DensityFunction );
  itkNewMacro( Self );

  /** Typedef alias for the measurement vectors */
  typedef TMeasurementVector MeasurementVectorType;

  /** Length of each measurement vector */
  typedef typename
    Superclass::MeasurementVectorSizeType            MeasurementVectorSizeType;

  typedef double                                     RealType;

  /** Type of the mean and measurement type */
  typedef Array<RealType>                            MeanType;
  typedef Array<RealType>                            MeasurementType;

  /** Type of the covariance matrix */
  typedef VariableSizeMatrix<RealType>               MatrixType;

  /** Sets the mean */
  void SetMean( const MeanType mean )
    {
    if ( this->m_Mean != mean )
      {
      this->SetMeasurementVectorSize( mean.Size() );
      this->m_Mean = mean;
      this->Modified();
      }
    }

  /** Gets the mean */
  const MeanType GetMean() const
    {
    return this->m_Mean;
    }

  /** Sets the covariance matrix.
   * Also, this function calculates inverse covariance and pre factor of
   * Gaussian Distribution to speed up GetProbability */
  void SetCovariance( MatrixType cov );

  /** Gets the covariance matrix */
  MatrixType GetCovariance();
  MatrixType GetInverseCovariance();

  /** Get/Set the prefactor */
  itkSetMacro( PreFactor, RealType );
  itkGetConstMacro( PreFactor, RealType );

  /** Gets the probability density of a measurement vector. */
  RealType Evaluate( const MeasurementVectorType &measurement ) const;
  RealType Evaluate( MeasurementType &measurement );

protected:
  GaussianClassDensityFunction( void );
  virtual ~GaussianClassDensityFunction( void ) {}
  void PrintSelf( std::ostream& os, Indent indent ) const;

private:
  MeanType        m_Mean;           // mean
  MatrixType      m_Covariance;     // covariance matrix

  // inverse covariance matrix which is automatically calculated
  // when covariance matrix is set.
  MatrixType  m_InverseCovariance;

  // Factor which is automatically calculated when covariace matrix is set.
  RealType m_NormalizationFactor;

  // Additional factor for multiplication.
  RealType m_PreFactor;

  /** if the all element of the given covarinace is zero, then this
   * value set to true */
  bool m_IsCovarianceZero;

  bool m_IsInverseCalculated;

};

} // end of namespace Statistics
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGaussianClassDensityFunction.hxx"
#endif

#endif

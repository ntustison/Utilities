/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGaussianClassDensityFunction.hxx,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGaussianClassDensityFunction_hxx
#define __itkGaussianClassDensityFunction_hxx

#include "itkGaussianClassDensityFunction.h"

#include "itkDecomposeTensorFunction.h"

namespace itk {

namespace Statistics {

template <class TMeasurementVector>
GaussianClassDensityFunction<TMeasurementVector>
::GaussianClassDensityFunction()
{
  this->m_Covariance.SetSize( this->GetMeasurementVectorSize(),
    this->GetMeasurementVectorSize() );
  this->m_Covariance.SetIdentity();
  this->m_InverseCovariance.SetSize( this->GetMeasurementVectorSize(),
    this->GetMeasurementVectorSize() );
  this->m_InverseCovariance.SetIdentity();

  this->m_IsInverseCalculated = false;

  this->m_Mean.SetSize( this->GetMeasurementVectorSize() );
  this->m_Mean.Fill( 0.0 );

  this->m_PreFactor = 1.0;
  this->m_NormalizationFactor = 1.0;
}

template <class TMeasurementVector>
void
GaussianClassDensityFunction<TMeasurementVector>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Mean: ";
  os << this->m_Mean << std::endl;

  os << indent << "Covariance: " << std::endl;
  os << this->m_Covariance.GetVnlMatrix();
  if( this->m_IsInverseCalculated )
    {
    os << indent << "InverseCovariance: " << std::endl;
    os << indent << this->m_InverseCovariance.GetVnlMatrix();
    }

  os << indent << "Prefactor: " << this->m_PreFactor << std::endl;
}

template <class TMeasurementVector>
void
GaussianClassDensityFunction<TMeasurementVector>
::SetCovariance( MatrixType cov )
{
  if( this->m_Covariance != cov )
    {
    this->m_Covariance = cov;
    this->m_IsInverseCalculated = false;
    this->SetMeasurementVectorSize( cov.GetVnlMatrix().rows() );
    this->Modified();
    }
}

template <class TMeasurementVector>
typename GaussianClassDensityFunction<TMeasurementVector>::MatrixType
GaussianClassDensityFunction<TMeasurementVector>
::GetCovariance()
{
  return this->m_Covariance;
}

template <class TMeasurementVector>
typename GaussianClassDensityFunction<TMeasurementVector>::MatrixType
GaussianClassDensityFunction<TMeasurementVector>
::GetInverseCovariance()
{
  return this->m_InverseCovariance;
}

template <class TMeasurementVector>
inline typename GaussianClassDensityFunction<TMeasurementVector>::RealType
GaussianClassDensityFunction<TMeasurementVector>
::Evaluate( const MeasurementVectorType &measurement ) const
{
  itkExceptionMacro( "Not implemented." );
}

template <class TMeasurementVector>
inline typename GaussianClassDensityFunction<TMeasurementVector>::RealType
GaussianClassDensityFunction<TMeasurementVector>
::Evaluate( MeasurementType &measurement )
{
  if( this->m_IsInverseCalculated == false )
    {
    this->m_IsCovarianceZero = this->m_Covariance.GetVnlMatrix().is_zero();

    if( !this->m_IsCovarianceZero )
      {
      this->m_InverseCovariance.GetVnlMatrix() =
        vnl_matrix_inverse<RealType>( this->m_Covariance.GetVnlMatrix() );

      // the determinant of the covaraince matrix
      typedef DecomposeTensorFunction<MatrixType, RealType> DecomposerType;
      typename DecomposerType::Pointer decomposer = DecomposerType::New();
      RealType det = decomposer->EvaluateDeterminant( this->m_Covariance );

      // calculate coefficient C of multivariate gaussian
      this->m_NormalizationFactor = 1.0 / ( vcl_sqrt( det ) *
        vcl_pow( sqrt( 2.0 * vnl_math::pi ),
        RealType( this->GetMeasurementVectorSize() ) ) );

      this->m_IsInverseCalculated = true;
      }
    }

  const MeasurementVectorSizeType measurementVectorSize
    = this->GetMeasurementVectorSize();

  MeanType tempVector;
  tempVector.SetSize( measurementVectorSize );

  // Compute |y - mean |
  for( unsigned int i = 0; i < measurementVectorSize; i++ )
    {
    tempVector[i] = measurement[i] - this->m_Mean[i];
    }

  MeanType tempVector2;
  tempVector2.SetSize( measurementVectorSize );

  // Compute |y - mean | * inverse(cov)
  for( unsigned int i = 0; i < measurementVectorSize; i++ )
    {
    RealType temp = 0.0;
    for( unsigned int j = 0; j < measurementVectorSize; j++ )
      {
      temp += tempVector[j]
        * this->m_InverseCovariance.GetVnlMatrix().get(j, i);
      }
    tempVector2[i] = temp;
    }

  // Compute |y - mean | * inverse(cov) * |y - mean|^T
  RealType temp = 0.0;
  for( unsigned int i = 0; i < measurementVectorSize; i++ )
    {
    temp += tempVector2[i] * tempVector[i];
    }
  return this->m_PreFactor * this->m_NormalizationFactor * vcl_exp( -0.5 * temp );
}

} // end namespace Statistics
} // end of namespace itk

#endif

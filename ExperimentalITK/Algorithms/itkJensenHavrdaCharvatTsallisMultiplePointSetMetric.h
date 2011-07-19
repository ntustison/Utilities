/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkJensenHavrdaCharvatTsallisMultiplePointSetMetric.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:13:42 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkJensenHavrdaCharvatTsallisMultiplePointSetMetric_h
#define __itkJensenHavrdaCharvatTsallisMultiplePointSetMetric_h

#include "itkMultiplePointSetMetric.h"

#include "itkIdentityTransform.h"
#include "itkManifoldParzenWindowsPointSetFunction.h"

#include <vector>

namespace itk {

/** \class JensenHavrdaCharvatTsallisMultiplePointSetMetric
 *
 *
 *
 *
 */
template<class TPointSet>
class ITK_EXPORT JensenHavrdaCharvatTsallisMultiplePointSetMetric :
    public MultiplePointSetMetric<TPointSet>
{
public:
  /** Standard class typedefs. */
  typedef JensenHavrdaCharvatTsallisMultiplePointSetMetric       Self;
  typedef MultiplePointSetMetric<TPointSet>                      Superclass;
  typedef SmartPointer<Self>                                     Pointer;
  typedef SmartPointer<const Self>                               ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods) */
  itkTypeMacro( JensenHavrdaCharvatTsallisMultiplePointSetMetric,
    MultiplePointSetMetric );

  itkStaticConstMacro( PointDimension, unsigned int,
                       TPointSet::PointDimension );

  /** Types transferred from the base class */
  typedef typename Superclass::RealType                 RealType;
  typedef typename Superclass::MeasureType              MeasureType;
  typedef typename Superclass::DerivativeType           DerivativeType;

  typedef TPointSet                                     PointSetType;
  typedef typename PointSetType::Pointer                PointSetPointer;
  typedef typename PointSetType::PointType              PointType;

  /**
   * Other typedefs
   */
  typedef ManifoldParzenWindowsPointSetFunction
    <PointSetType, RealType>                            DensityFunctionType;
  typedef typename DensityFunctionType::Pointer         DensityFunctionPointer;
  typedef typename DensityFunctionType::GaussianType    GaussianType;


  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void Initialize( void ) throw ( ExceptionObject );

  /** Get the number of values */
  unsigned int GetNumberOfValues() const;

  /** Get the derivatives of the match measure. */
  virtual void GetDerivative( DerivativeType & Derivative ) const;

  /**  Get the value for single valued optimizers. */
  virtual MeasureType GetValue() const;

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative( MeasureType& Value,
    DerivativeType& Derivative ) const;

  itkSetClampMacro( Alpha, RealType, 1.0, 2.0 );
  itkGetConstMacro( Alpha, RealType );

  itkSetMacro( UseRegularizationTerm, bool );
  itkGetConstMacro( UseRegularizationTerm, bool );
  itkBooleanMacro( UseRegularizationTerm );

  itkSetMacro( UseInputAsSamples, bool );
  itkGetConstMacro( UseInputAsSamples, bool );
  itkBooleanMacro( UseInputAsSamples );

  void SetPointSetSigma( unsigned int i, RealType sigma )
    {
    if( i >= this->m_PointSetSigma.size() )
      {
      this->m_PointSetSigma.resize( i+1 );
      this->m_PointSetSigma[i] = sigma;
      this->Modified();
      }
    else if( this->m_PointSetSigma[i] != sigma )
      {
      this->m_PointSetSigma[i] = sigma;
      this->Modified();
      }
    }
  RealType GetPointSetSigma( unsigned int i )
    {
    if( i < this->m_PointSetSigma.size() )
      {
      return this->m_PointSetSigma[i];
      }
    return 0;
    }

  void SetEvaluationKNeighborhood( unsigned int i, unsigned int n )
    {
    if( i >= this->m_EvaluationKNeighborhood.size() )
      {
      this->m_EvaluationKNeighborhood.resize( i+1 );
      this->m_EvaluationKNeighborhood[i] = n;
      this->Modified();
      }
    else if( this->m_EvaluationKNeighborhood[i] != n )
      {
      this->m_EvaluationKNeighborhood[i] = n;
      this->Modified();
      }
    }
  unsigned int GetEvaluationKNeighborhood( unsigned int i )
    {
    if( i < this->m_EvaluationKNeighborhood.size() )
      {
      return this->m_EvaluationKNeighborhood[i];
      }
    }

  /**
   * If this->m_UseInputAsSamples = true, the following
   * two variables are not used.
   */

  void SetNumberOfSamples( unsigned int i, unsigned long n )
    {
    if( i >= this->m_NumberOfSamples.size() )
      {
      this->m_NumberOfSamples.resize( i+1 );
      this->m_NumberOfSamples[i] = n;
      this->Modified();
      }
    else if( this->m_NumberOfSamples[i] != n )
      {
      this->m_NumberOfSamples[i] = n;
      this->Modified();
      }
    }
  unsigned long GetNumberOfSamples( unsigned int i )
    {
    if( i < this->m_NumberOfSamples.size() )
      {
      return this->m_NumberOfSamples[i];
      }
    }

  itkSetMacro( UseAnisotropicCovariances, bool );
  itkGetConstMacro( UseAnisotropicCovariances, bool );
  itkBooleanMacro( UseAnisotropicCovariances );

  /**
   * If this->m_UseAnisotropicCovariances = false, the
   * following four variables are not used.
   */

  void SetCovarianceKNeighborhood( unsigned int i, unsigned int n )
    {
    if( i >= this->m_CovarianceKNeighborhood.size() )
      {
      this->m_CovarianceKNeighborhood.resize( i+1 );
      this->m_CovarianceKNeighborhood[i] = n;
      this->Modified();
      }
    else if( this->m_CovarianceKNeighborhood[i] != n )
      {
      this->m_CovarianceKNeighborhood[i] = n;
      this->Modified();
      }
    }
  unsigned int GetCovarianceKNeighborhood( unsigned int i )
    {
    if( i < this->m_CovarianceKNeighborhood.size() )
      {
      return this->m_CovarianceKNeighborhood[i];
      }
    }

  void SetKernelSigma( unsigned int i, RealType sigma )
    {
    if( i >= this->m_KernelSigma.size() )
      {
      this->m_KernelSigma.resize( i+1 );
      this->m_KernelSigma[i] = sigma;
      this->Modified();
      }
    else if( this->m_KernelSigma[i] != sigma )
      {
      this->m_KernelSigma[i] = sigma;
      this->Modified();
      }
    }
  RealType GetKernelSigma( unsigned int i )
    {
    if( i < this->m_KernelSigma.size() )
      {
      return this->m_KernelSigma[i];
      }
    }

protected:
  JensenHavrdaCharvatTsallisMultiplePointSetMetric();
  ~JensenHavrdaCharvatTsallisMultiplePointSetMetric() {}

  void PrintSelf( std::ostream& os, Indent indent ) const;

private:
  //purposely not implemented
  JensenHavrdaCharvatTsallisMultiplePointSetMetric(const Self&);
  void operator=(const Self&);

  bool                                     m_UseRegularizationTerm;
  bool                                     m_UseInputAsSamples;
  bool                                     m_UseAnisotropicCovariances;

  std::vector<DensityFunctionPointer>      m_DensityFunction;
  std::vector<PointSetPointer>             m_SamplePoints;
  std::vector<RealType>                    m_PointSetSigma;
  std::vector<RealType>                    m_KernelSigma;
  std::vector<unsigned int>                m_CovarianceKNeighborhood;
  std::vector<unsigned int>                m_EvaluationKNeighborhood;
  std::vector<unsigned long>               m_NumberOfSamples;

  RealType                                 m_Alpha;

};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkJensenHavrdaCharvatTsallisMultiplePointSetMetric.hxx"
#endif

#endif

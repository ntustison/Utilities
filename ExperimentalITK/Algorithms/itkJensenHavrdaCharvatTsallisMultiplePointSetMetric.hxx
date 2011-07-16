/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkJensenHavrdaCharvatTsallisMultiplePointSetMetric.hxx,v $
  Language:  C++
  Date:      $Date: 2008/11/03 15:19:47 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkJensenHavrdaCharvatTsallisMultiplePointSetMetric_hxx
#define __itkJensenHavrdaCharvatTsallisMultiplePointSetMetric_hxx

#include "itkJensenHavrdaCharvatTsallisMultiplePointSetMetric.h"

namespace itk {

template <class TPointSet>
JensenHavrdaCharvatTsallisMultiplePointSetMetric<TPointSet>
::JensenHavrdaCharvatTsallisMultiplePointSetMetric()
{
  this->m_UseRegularizationTerm = false;
  this->m_UseInputAsSamples = true;
  this->m_UseAnisotropicCovariances = false;
  this->m_Alpha = 2.0;

  this->SetNumberOfRequiredInputs( 2 );
}

/** Initialize the metric */
template <class TPointSet>
void
JensenHavrdaCharvatTsallisMultiplePointSetMetric<TPointSet>
::Initialize( void ) throw ( ExceptionObject )
{
  Superclass::Initialize();

  /**
   * Ensure that each input point-set has a corresponding group of
   * parameters.
   */
  std::vector<RealType>                    m_PointSetSigma;
  std::vector<RealType>                    m_KernelSigma;
  std::vector<unsigned int>                m_CovarianceKNeighborhood;
  std::vector<unsigned int>                m_EvaluationKNeighborhood;
  std::vector<unsigned long>               m_NumberOfSamples;

  if( this->GetUseAnisotropicCovariances() &&
      this->m_KernelSigma.size() < this->GetNumberOfInputs() )
    {
    itkExceptionMacro( "m_KernelSigma only has " << this->m_KernelSigma.size()
      << "elements." );
    }
  if( this->m_PointSetSigma.size() < this->GetNumberOfInputs() )
    {
    itkExceptionMacro( "m_PointSetSigma only has "
      << this->m_PointSetSigma.size() << "elements." );
    }
  if( this->GetUseAnisotropicCovariances() &&
      this->m_CovarianceKNeighborhood.size() < this->GetNumberOfInputs() )
    {
    itkExceptionMacro( "m_CovarianceKNeighborhood only has "
      << this->m_CovarianceKNeighborhood.size() << "elements." );
    }
  if( this->m_EvaluationKNeighborhood.size() < this->GetNumberOfInputs() )
    {
    itkExceptionMacro( "m_EvaluationKNeighborhood only has "
      << this->m_EvaluationKNeighborhood.size() << "elements." );
    }
  if( !this->m_UseInputAsSamples &&
      this->m_NumberOfSamples.size() < this->GetNumberOfInputs() )
    {
    itkExceptionMacro( "m_NumberOfSamples only has "
      << this->m_NumberOfSamples.size() << "elements." );
    }

  /**
   * Initialize the points
   */
  this->m_DensityFunction.clear();
  this->m_SamplePoints.clear();
  for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
    {
    DensityFunctionPointer densityFunction = DensityFunctionType::New();
    densityFunction->SetBucketSize( 4 );
    densityFunction->SetRegularizationSigma( this->m_PointSetSigma[i] );
    densityFunction->SetNormalize( true );
    densityFunction->SetUseAnisotropicCovariances(
      this->m_UseAnisotropicCovariances );
    if( this->m_UseAnisotropicCovariances )
      {
      densityFunction->SetCovarianceKNeighborhood(
        this->m_CovarianceKNeighborhood[i] );
      densityFunction->SetKernelSigma( this->m_KernelSigma[i] );
      }
    densityFunction->SetEvaluationKNeighborhood(
      this->m_EvaluationKNeighborhood[i] );
    densityFunction->SetInputPointSet( this->GetInput( i ) );
    this->m_DensityFunction.push_back( densityFunction );

    PointSetPointer samplePoints = PointSetType::New();
    samplePoints->Initialize();
    
    if( this->m_UseInputAsSamples )
      {
      typename PointSetType::PointsContainerConstIterator It
        = this->GetInput( i )->GetPoints()->Begin();
      while( It != this->GetInput( i )->GetPoints()->End() )
        {
        samplePoints->SetPoint( It.Index(), It.Value() );
        ++It;
        }
      }
    else
      {
      for( unsigned long j = 0; j < this->m_NumberOfSamples[i]; j++ )
        {
        samplePoints->SetPoint( j,
          this->m_DensityFunction[i]->GenerateRandomSample() );
        }
      }  
    this->m_SamplePoints.push_back( samplePoints );
    }
}

/** Return the number of values, i.e the total number of points */
template <class TPointSet>
unsigned int
JensenHavrdaCharvatTsallisMultiplePointSetMetric<TPointSet>
::GetNumberOfValues() const
{
  unsigned int numberOfValues = 0;
  for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
    {
    PointSetPointer points = const_cast<PointSetType*>(
      static_cast<const PointSetType*>( this->ProcessObject::GetInput( i ) ) );
    numberOfValues += points->GetNumberOfPoints();
    }
  return numberOfValues;
}


/** Get the match Measure */
template <class TPointSet>
typename JensenHavrdaCharvatTsallisMultiplePointSetMetric
  <TPointSet>::MeasureType
JensenHavrdaCharvatTsallisMultiplePointSetMetric<TPointSet>
::GetValue() const
{

  RealType totalNumberOfPoints = 0;
  for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
    {
    PointSetPointer points = const_cast<PointSetType*>(
      static_cast<const PointSetType*>( this->ProcessObject::GetInput( i ) ) );
    totalNumberOfPoints += points->GetNumberOfPoints();
    }
  RealType totalNumberOfSamples = 0;
  for( unsigned int i = 0; i < this->m_SamplePoints.size(); i++ )
    {
    totalNumberOfSamples += static_cast<RealType>(
      this->m_SamplePoints[i]->GetNumberOfPoints() );
    }

  PointSetPointer points[this->GetNumberOfInputs()];
  PointSetPointer samples[this->GetNumberOfInputs()];

  for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
    {
    points[i] = const_cast<PointSetType *>(
      static_cast<const PointSetType *>( this->ProcessObject::GetInput( i ) ) );
    if( this->m_UseInputAsSamples )
      {
      samples[i] = points[i];
      }
    else
      {
      samples[i] = this->m_SamplePoints[i];
      }
    }

  MeasureType measure;
  measure.SetSize( 1 );
  measure.Fill( 0 );

  RealType energyTerm1 = 0.0;
  RealType energyTerm2 = 0.0;

  /**
    * first term
    */
  for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
    {
    RealType prefactor = -1.0 / totalNumberOfSamples;
    if( this->m_Alpha != 1.0 )
      {
      prefactor /= ( this->m_Alpha - 1.0 );
      }

    typename PointSetType::PointsContainerConstIterator It
      = samples[i]->GetPoints()->Begin();
    while( It != samples[i]->GetPoints()->End() )
      {
      PointType samplePoint = It.Value();

      RealType probabilityStar = 0.0;
      for( unsigned int j = 0; j < this->m_DensityFunction.size(); j++ )
        {
        probabilityStar += this->m_DensityFunction[j]->Evaluate( samplePoint )
          * static_cast<RealType>( points[j]->GetNumberOfPoints() );
        }
      probabilityStar /= totalNumberOfPoints;

      if( probabilityStar == 0 )
        {
        ++It;
        continue;
        }

      if( this->m_Alpha == 1.0 )
        {
        energyTerm1 += ( prefactor *
          vcl_log( probabilityStar )  / vnl_math::ln2 );
        }
      else
        {
        energyTerm1 += ( prefactor * vcl_pow( probabilityStar,
          static_cast<RealType>( this->m_Alpha - 1.0 ) ) );
        }
      ++It;
      }
    }

  /**
    * second term, i.e. regularization term
    */
  if( this->m_UseRegularizationTerm )
    {
    for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
      {
      RealType prefactor = -static_cast<RealType>(
        points[i]->GetNumberOfPoints() ) / ( totalNumberOfPoints *
        static_cast<RealType>( samples[i]->GetNumberOfPoints() ) );
      if( this->m_Alpha != 1.0 )
        {
        prefactor /= ( this->m_Alpha - 1.0 );
        }

      typename PointSetType::PointsContainerConstIterator It
        = samples[i]->GetPoints()->Begin();
      while( It != samples[i]->GetPoints()->End() )
        {
        PointType samplePoint = It.Value();

        RealType probability
          = this->m_DensityFunction[i]->Evaluate( samplePoint );

        if( probability == 0 )
          {
          ++It;
          continue;
          }

        if( this->m_Alpha == 1.0 )
          {
          energyTerm2 += ( prefactor * vcl_log( probability ) / vnl_math::ln2 );
          }
        else
          {
          energyTerm2 += ( prefactor * vcl_pow( probability,
            static_cast<RealType>( this->m_Alpha - 1.0 ) ) );
          }
        ++It;
        }
      }
    }

  measure[0] = energyTerm1 - energyTerm2;

  return measure;
}

/** Get the Derivative */
template <class TPointSet>
void
JensenHavrdaCharvatTsallisMultiplePointSetMetric<TPointSet>
::GetDerivative( DerivativeType & derivative ) const
{

  RealType totalNumberOfPoints = 0;
  for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
    {
    PointSetPointer points = const_cast<PointSetType*>(
      static_cast<const PointSetType*>( this->ProcessObject::GetInput( i ) ) );
    totalNumberOfPoints += points->GetNumberOfPoints();
    }
  RealType totalNumberOfSamples = 0;
  for( unsigned int i = 0; i < this->m_SamplePoints.size(); i++ )
    {
    totalNumberOfSamples += static_cast<RealType>(
      this->m_SamplePoints[i]->GetNumberOfPoints() );
    }

  PointSetPointer points[this->GetNumberOfInputs()];
  PointSetPointer samples[this->GetNumberOfInputs()];

  for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
    {
    points[i] = const_cast<PointSetType *>(
      static_cast<const PointSetType *>( this->ProcessObject::GetInput( i ) ) );
    if( this->m_UseInputAsSamples )
      {
      samples[i] = points[i];
      }
    else
      {
      samples[i] = this->m_SamplePoints[i];
      }
    }

  derivative.SetSize( this->GetNumberOfValues(), PointDimension );
  derivative.Fill( 0 );

  /**
   * first term
   */
  for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
    {
    unsigned int index = 0;
    for ( unsigned int j = 0; j < i; j++ )
      {
      PointSetPointer pts = const_cast<PointSetType*>(
        static_cast<const PointSetType*>( this->ProcessObject::GetInput( j ) ) );
      index += pts->GetNumberOfPoints();
      }
    RealType prefactor = 1.0 / ( totalNumberOfSamples * totalNumberOfPoints );

    for( unsigned int k = 0; k < this->GetNumberOfInputs(); k++ )
      {
      typename PointSetType::PointsContainerConstIterator It
        = samples[k]->GetPoints()->Begin();
      while( It != samples[k]->GetPoints()->End() )
        {
        PointType samplePoint = It.Value();
  
        RealType probabilityStar = 0.0;
        for( unsigned int j = 0; j < this->m_DensityFunction.size(); j++ )
          {
          probabilityStar += this->m_DensityFunction[j]->Evaluate( samplePoint )
            * static_cast<RealType>( points[j]->GetNumberOfPoints() );
          }
        probabilityStar /= totalNumberOfPoints;
  
        if( probabilityStar == 0 )
          {
          ++It;
          continue;
          }
  
        RealType probabilityStarFactor = vcl_pow( probabilityStar,
          static_cast<RealType>( 2.0 - this->m_Alpha ) );
  
        typename GaussianType::MeasurementVectorType sampleMeasurement;
        for( unsigned int d = 0; d < PointDimension; d++ )
          {
          sampleMeasurement[d] = samplePoint[d];
          }
  
        typename DensityFunctionType::NeighborhoodIdentifierType neighbors
          = this->m_DensityFunction[i]->GetNeighborhoodIdentifiers(
            sampleMeasurement, this->m_EvaluationKNeighborhood[i] );
  
        for( unsigned int n = 0; n < neighbors.size(); n++ )
          {
          RealType gaussian = this->m_DensityFunction[i]->
            GetGaussian( neighbors[n] )->Evaluate( sampleMeasurement );
  
          if( gaussian == 0 )
            {
            continue;
            }
  
          typename GaussianType::MeanType mean
            = this->m_DensityFunction[i]->GetGaussian( neighbors[n] )->GetMean();
  
          for( unsigned int d = 0; d < PointDimension; d++ )
            {
            mean[d] -= samplePoint[d];
            }
  
          if( this->m_UseAnisotropicCovariances )
            {
            typename GaussianType::MatrixType Ci
              = this->m_DensityFunction[i]->GetGaussian( neighbors[n] )
                ->GetInverseCovariance();
            mean = Ci * mean;
            }
          else
            {
            mean /= vnl_math_sqr( this->m_DensityFunction[i]->
              GetGaussian( neighbors[n] )->GetSigma() );
            }
  
          mean *= ( prefactor * gaussian / probabilityStarFactor );
          for( unsigned int d = 0; d < PointDimension; d++ )
            {
            derivative(index+neighbors[n], d) += mean[d];
            }
          }
        ++It;
        }
      }
    }
  /**
   * second term, i.e. regularization term
   */
  if( this->m_UseRegularizationTerm )
    {
    for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
      {
      unsigned int index = 0;
      for ( unsigned int j = 0; j < i; j++ )
        {
        PointSetPointer pts = const_cast<PointSetType*>(
          static_cast<const PointSetType*>( this->ProcessObject::GetInput( j ) ) );
        index += pts->GetNumberOfPoints();
        }
      RealType prefactor = -1.0 / ( static_cast<RealType>(
        samples[i]->GetNumberOfPoints() ) * totalNumberOfPoints );

      typename PointSetType::PointsContainerConstIterator It
        = samples[i]->GetPoints()->Begin();
      while( It != samples[i]->GetPoints()->End() )
        {
        PointType samplePoint = It.Value();

        RealType probability
          = this->m_DensityFunction[i]->Evaluate( samplePoint );

        if( probability == 0 )
          {
          ++It;
          continue;
          }

        RealType probabilityFactor = vcl_pow( probability,
          static_cast<RealType>( 2.0 - this->m_Alpha ) );

        typename GaussianType::MeasurementVectorType sampleMeasurement;
        for( unsigned int d = 0; d < PointDimension; d++ )
          {
          sampleMeasurement[d] = samplePoint[d];
          }

        typename DensityFunctionType::NeighborhoodIdentifierType neighbors
          = this->m_DensityFunction[i]->GetNeighborhoodIdentifiers(
            sampleMeasurement, this->m_EvaluationKNeighborhood[i] );

        for( unsigned int n = 0; n < neighbors.size(); n++ )
          {
          RealType gaussian = this->m_DensityFunction[i]->
            GetGaussian( neighbors[n] )->Evaluate( sampleMeasurement );
          if( gaussian == 0 )
            {
            continue;
            }

          typename GaussianType::MeanType mean
            = this->m_DensityFunction[i]->GetGaussian( neighbors[n] )->GetMean();
          for( unsigned int d = 0; d < PointDimension; d++ )
            {
            mean[d] -= samplePoint[d];
            }

          if( this->m_UseAnisotropicCovariances )
            {
            typename GaussianType::MatrixType Ci
              = this->m_DensityFunction[i]->GetGaussian( neighbors[n] )
                ->GetInverseCovariance();
            mean = Ci * mean;
            }
          else
            {
            mean /= vnl_math_sqr( this->m_DensityFunction[i]->
              GetGaussian( neighbors[n] )->GetSigma() );
            }

          mean *= ( prefactor * gaussian / probabilityFactor );
          for( unsigned int d = 0; d < PointDimension; d++ )
            {
            derivative(index+neighbors[n], d) += mean[d];
            }
          }
        ++It;
        }
      }
    }
}

/** Get both the match Measure and theDerivative Measure  */
template <class TPointSet>
void
JensenHavrdaCharvatTsallisMultiplePointSetMetric<TPointSet>
::GetValueAndDerivative( MeasureType & value,
  DerivativeType  & derivative ) const
{

  RealType totalNumberOfPoints = 0;
  for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
    {
    PointSetPointer points = const_cast<PointSetType*>(
      static_cast<const PointSetType*>( this->ProcessObject::GetInput( i ) ) );
    totalNumberOfPoints += points->GetNumberOfPoints();
    }
  RealType totalNumberOfSamples = 0;
  for( unsigned int i = 0; i < this->m_SamplePoints.size(); i++ )
    {
    totalNumberOfSamples += static_cast<RealType>(
      this->m_SamplePoints[i]->GetNumberOfPoints() );
    }

  PointSetPointer points[this->GetNumberOfInputs()];
  PointSetPointer samples[this->GetNumberOfInputs()];

  for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
    {
    points[i] = const_cast<PointSetType *>(
      static_cast<const PointSetType *>( this->ProcessObject::GetInput( i ) ) );
    if( this->m_UseInputAsSamples )
      {
      samples[i] = points[i];
      }
    else
      {
      samples[i] = this->m_SamplePoints[i];
      }
    }

  derivative.SetSize( this->GetNumberOfValues(), PointDimension );
  derivative.Fill( 0 );

  value.SetSize( 1 );
  value.Fill( 0 );

  RealType energyTerm1 = 0.0;
  RealType energyTerm2 = 0.0;

  /**
   * first term
   */
  for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
    {
    unsigned int index = 0;
    for ( unsigned int j = 0; j < i; j++ )
      {
      PointSetPointer pts = const_cast<PointSetType*>(
        static_cast<const PointSetType*>( this->ProcessObject::GetInput( j ) ) );
      index += pts->GetNumberOfPoints();
      }
    RealType prefactor[2];
    prefactor[0] = -1.0 / totalNumberOfSamples;
    if( this->m_Alpha != 1.0 )
      {
      prefactor[0] /= ( this->m_Alpha - 1.0 );
      }
    prefactor[1] = 1.0 / ( totalNumberOfSamples * totalNumberOfPoints );

    for( unsigned int k = 0; k < this->GetNumberOfInputs(); k++ )
      {
      typename PointSetType::PointsContainerConstIterator It
        = samples[k]->GetPoints()->Begin();
      while( It != samples[k]->GetPoints()->End() )
        {
        PointType samplePoint = It.Value();
  
        RealType probabilityStar = 0.0;
        for( unsigned int j = 0; j < this->m_DensityFunction.size(); j++ )
          {
          probabilityStar += this->m_DensityFunction[j]->Evaluate( samplePoint )
            * static_cast<RealType>( points[j]->GetNumberOfPoints() );
          }
        probabilityStar /= totalNumberOfPoints;
  
        if( probabilityStar == 0 )
          {
          ++It;
          continue;
          }
  
        if( this->m_Alpha == 1.0 )
          {
          energyTerm1 += ( prefactor[0] *
            vcl_log( probabilityStar )  / vnl_math::ln2 );
          }
        else
          {
          energyTerm1 += ( prefactor[0] * vcl_pow( probabilityStar,
            static_cast<RealType>( this->m_Alpha - 1.0 ) ) );
          }
  
        RealType probabilityStarFactor = vcl_pow( probabilityStar,
          static_cast<RealType>( 2.0 - this->m_Alpha ) );
  
        typename GaussianType::MeasurementVectorType sampleMeasurement;
        for( unsigned int d = 0; d < PointDimension; d++ )
          {
          sampleMeasurement[d] = samplePoint[d];
          }
  
        typename DensityFunctionType::NeighborhoodIdentifierType neighbors
          = this->m_DensityFunction[i]->GetNeighborhoodIdentifiers(
            sampleMeasurement, this->m_EvaluationKNeighborhood[i] );
  
        for( unsigned int n = 0; n < neighbors.size(); n++ )
          {
          RealType gaussian = this->m_DensityFunction[i]->
            GetGaussian( neighbors[n] )->Evaluate( sampleMeasurement );
  
          if( gaussian == 0 )
            {
            continue;
            }
  
          typename GaussianType::MeanType mean
            = this->m_DensityFunction[i]->GetGaussian( neighbors[n] )->GetMean();
  
          for( unsigned int d = 0; d < PointDimension; d++ )
            {
            mean[d] -= samplePoint[d];
            }
  
          if( this->m_UseAnisotropicCovariances )
            {
            typename GaussianType::MatrixType Ci
              = this->m_DensityFunction[i]->GetGaussian( neighbors[n] )
                ->GetInverseCovariance();
            mean = Ci * mean;
            }
          else
            {
            mean /= vnl_math_sqr( this->m_DensityFunction[i]->
              GetGaussian( neighbors[n] )->GetSigma() );
            }
  
          mean *= ( prefactor[1] * gaussian / probabilityStarFactor );
          for( unsigned int d = 0; d < PointDimension; d++ )
            {
            derivative(index+neighbors[n], d) += mean[d];
            }
          }
        ++It;
        }
      }  
    }

  /**
   * second term, i.e. regularization term
   */
  if( this->m_UseRegularizationTerm )
    {
    for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
      {
      unsigned int index = 0;
      for ( unsigned int j = 0; j < i; j++ )
        {
        PointSetPointer pts = const_cast<PointSetType*>(
          static_cast<const PointSetType*>( this->ProcessObject::GetInput( j ) ) );
        index += pts->GetNumberOfPoints();
        }
      RealType prefactor[2];
      prefactor[0] = -static_cast<RealType>(
        points[i]->GetNumberOfPoints() ) / ( totalNumberOfPoints *
        static_cast<RealType>( samples[i]->GetNumberOfPoints() ) );
      if( this->m_Alpha != 1.0 )
        {
        prefactor[0] /= ( this->m_Alpha - 1.0 );
        }
      prefactor[1] = -1.0 / ( static_cast<RealType>(
        samples[i]->GetNumberOfPoints() ) * totalNumberOfPoints );

      typename PointSetType::PointsContainerConstIterator It
        = samples[i]->GetPoints()->Begin();
      while( It != samples[i]->GetPoints()->End() )
        {
        PointType samplePoint = It.Value();

        RealType probability
          = this->m_DensityFunction[i]->Evaluate( samplePoint );

        if( probability == 0 )
          {
          ++It;
          continue;
          }

        if( this->m_Alpha == 1.0 )
          {
          energyTerm2 += ( prefactor[0]
            * vcl_log( probability ) / vnl_math::ln2 );
          }
        else
          {
          energyTerm2 += ( prefactor[0] * vcl_pow( probability,
            static_cast<RealType>( this->m_Alpha - 1.0 ) ) );
          }

        RealType probabilityFactor = vcl_pow( probability,
          static_cast<RealType>( 2.0 - this->m_Alpha ) );

        typename GaussianType::MeasurementVectorType sampleMeasurement;
        for( unsigned int d = 0; d < PointDimension; d++ )
          {
          sampleMeasurement[d] = samplePoint[d];
          }

        typename DensityFunctionType::NeighborhoodIdentifierType neighbors
          = this->m_DensityFunction[i]->GetNeighborhoodIdentifiers(
            sampleMeasurement, this->m_EvaluationKNeighborhood[i] );

        for( unsigned int n = 0; n < neighbors.size(); n++ )
          {
          RealType gaussian = this->m_DensityFunction[i]->
            GetGaussian( neighbors[n] )->Evaluate( sampleMeasurement );
          if( gaussian == 0 )
            {
            continue;
            }

          typename GaussianType::MeanType mean
            = this->m_DensityFunction[i]->GetGaussian( neighbors[n] )->GetMean();
          for( unsigned int d = 0; d < PointDimension; d++ )
            {
            mean[d] -= samplePoint[d];
            }

          if( this->m_UseAnisotropicCovariances )
            {
            typename GaussianType::MatrixType Ci
              = this->m_DensityFunction[i]->GetGaussian( neighbors[n] )
                ->GetInverseCovariance();
            mean = Ci * mean;
            }
          else
            {
            mean /= vnl_math_sqr( this->m_DensityFunction[i]->
              GetGaussian( neighbors[n] )->GetSigma() );
            }

          mean *= ( prefactor[1] * gaussian / probabilityFactor );
          for( unsigned int d = 0; d < PointDimension; d++ )
            {
            derivative(index+neighbors[n], d) += mean[d];
            }
          }
        ++It;
        }
      }
    }

  value[0] = energyTerm1 - energyTerm2;
}

template <class TPointSet>
void
JensenHavrdaCharvatTsallisMultiplePointSetMetric<TPointSet>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Use regularization term: "
     << this->m_UseRegularizationTerm << std::endl;
  os << indent << "Alpha: "
     << this->m_Alpha << std::endl;
  if( !this->m_UseAnisotropicCovariances )
    {
    os << indent << "Isotropic covariances are used." << std::endl;
    }
  if( this->m_UseInputAsSamples )
    {
    os << indent << "Use input points as samples." << std::endl;
    }

  for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
    {
    os << indent << " Point-set " << i << ":" << std::endl;
    PointSetPointer pts = const_cast<PointSetType*>(
      static_cast<const PointSetType*>( this->ProcessObject::GetInput( i ) ) );
    os << indent << "   Number of points = " 
       << pts->GetNumberOfPoints() << std::endl;
    os << indent << "   Point-set sigma = "
       << this->m_PointSetSigma[i] << std::endl;
    os << indent << "   Evaluation k-neighborhood: "
        << this->m_EvaluationKNeighborhood[i] << std::endl;
    if( !this->m_UseInputAsSamples )
      {
      os << indent << "   Number of samples = "
         << this->m_NumberOfSamples[i] << std::endl;
      }
    if( this->m_UseAnisotropicCovariances )
      {
      os << indent << "   Kernel sigma = "
         << this->m_KernelSigma[i] << std::endl;
      os << indent << "   Covariance k-neighborhood: "
         << this->m_CovarianceKNeighborhood[i] << std::endl;
      }
    }
}

} // end namespace itk


#endif

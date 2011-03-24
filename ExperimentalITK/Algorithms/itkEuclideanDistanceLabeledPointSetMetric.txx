/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkEuclideanDistanceLabeledPointSetMetric.txx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:13:42 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkEuclideanDistanceLabeledPointSetMetric_txx
#define _itkEuclideanDistanceLabeledPointSetMetric_txx

#include "itkEuclideanDistanceLabeledPointSetMetric.h"

#include "itkEuclideanDistancePointSetMetric.h"
namespace itk {

template <class TPointSet>
EuclideanDistanceLabeledPointSetMetric<TPointSet>
::EuclideanDistanceLabeledPointSetMetric()
{
  this->m_UseWithRespectToTheMovingPointSet = true;

  typename DefaultTransformType::Pointer transform
    = DefaultTransformType::New();
  transform->SetIdentity();

  Superclass::SetTransform( transform );
}

/** Initialize the metric */
template <class TPointSet>
void
EuclideanDistanceLabeledPointSetMetric<TPointSet>
::Initialize( void ) throw ( ExceptionObject )
{
  Superclass::Initialize();

  if ( this->m_FixedLabelSet.size() <= 0
       && this->m_FixedPointSet->GetNumberOfPoints() > 0 )
    {
    typename PointSetType::PointDataContainerIterator It
      = this->m_FixedPointSet->GetPointData()->Begin();
    while ( It != this->m_FixedPointSet->GetPointData()->End() )
      {
      if ( find( this->m_FixedLabelSet.begin(), this->m_FixedLabelSet.end(),
             It.Value() ) == this->m_FixedLabelSet.end() )
        {
        this->m_FixedLabelSet.push_back( It.Value() );
        }
      ++It;
      }
    }
  if ( this->m_MovingLabelSet.size() <= 0
       && this->m_MovingPointSet->GetNumberOfPoints() > 0 )
    {
    typename PointSetType::PointDataContainerIterator It
      = this->m_MovingPointSet->GetPointData()->Begin();
    while ( It != this->m_MovingPointSet->GetPointData()->End() )
      {
      if ( find( this->m_MovingLabelSet.begin(), this->m_MovingLabelSet.end(),
             It.Value() ) == this->m_MovingLabelSet.end() )
        {
        this->m_MovingLabelSet.push_back( It.Value() );
        }
      ++It;
      }
    }
}

/** Return the number of values, i.e the number of points in the moving set */
template <class TPointSet>
unsigned int
EuclideanDistanceLabeledPointSetMetric<TPointSet>
::GetNumberOfValues() const
{
  if ( this->m_UseWithRespectToTheMovingPointSet )
    {
    if ( this->m_MovingPointSet )
      {
      return  this->m_MovingPointSet->GetPoints()->Size();
      }
    }
  else
    {
    if ( this->m_FixedPointSet )
      {
      return  this->m_FixedPointSet->GetPoints()->Size();
      }
    }

  return 0;
}


/** Get the match Measure */
template <class TPointSet>
typename EuclideanDistanceLabeledPointSetMetric<TPointSet>::MeasureType
EuclideanDistanceLabeledPointSetMetric<TPointSet>
::GetValue( const TransformParametersType & parameters ) const
{
  MeasureType measure;
  measure.SetSize( 1 );
  measure.Fill( 0 );

  typename LabelSetType::const_iterator iter;

  for ( iter = this->m_FixedLabelSet.begin();
        iter != this->m_FixedLabelSet.end(); ++iter )
    {
    PixelType currentLabel = *iter;

    /**
     * check to see if the moving label set contains the same label
     */
    if ( find( this->m_MovingLabelSet.begin(), this->m_MovingLabelSet.end(),
           currentLabel ) == this->m_MovingLabelSet.end() )
      {
      continue;
      }

    /**
     * Collect all the fixed and moving points with the currentLabel
     */
    typename PointSetType::Pointer fixedLabelPoints
      = PointSetType::New();
    fixedLabelPoints->Initialize();
    unsigned long fixedCount = 0;

    typename PointSetType::PointsContainerConstIterator ItF =
      this->m_FixedPointSet->GetPoints()->Begin();
    typename PointSetType::PointDataContainerIterator ItFD =
      this->m_FixedPointSet->GetPointData()->Begin();

    while ( ItF != this->m_FixedPointSet->GetPoints()->End() )
      {
      if ( ItFD.Value() == currentLabel )
        {
        fixedLabelPoints->SetPoint( fixedCount++, ItF.Value() );
        }
      ++ItF;
      ++ItFD;
      }

    typename PointSetType::Pointer movingLabelPoints
      = PointSetType::New();
    movingLabelPoints->Initialize();
    unsigned long movingCount = 0;

    typename PointSetType::PointsContainerConstIterator ItM =
      this->m_MovingPointSet->GetPoints()->Begin();
    typename PointSetType::PointDataContainerIterator ItMD =
      this->m_MovingPointSet->GetPointData()->Begin();

    while ( ItM != this->m_MovingPointSet->GetPoints()->End() )
      {
      if ( ItMD.Value() == currentLabel )
        {
        movingLabelPoints->SetPoint( movingCount++, ItM.Value() );
        }
      ++ItM;
      ++ItMD;
      }

    /**
     * Invoke the single label JensenTsallis measure
     */
    typedef EuclideanDistancePointSetMetric
      <PointSetType> MetricType;
    typename MetricType::Pointer metric = MetricType::New();

    metric->SetFixedPointSet( fixedLabelPoints );
    metric->SetMovingPointSet( movingLabelPoints );
    metric->SetUseWithRespectToTheMovingPointSet(
      this->m_UseWithRespectToTheMovingPointSet );
    metric->Initialize();

    MeasureType value = metric->GetValue( parameters );
    measure[0] += value[0];
    }

  return measure;
}

/** Get the Derivative Measure */
template <class TPointSet>
void
EuclideanDistanceLabeledPointSetMetric<TPointSet>
::GetDerivative( const TransformParametersType & parameters,
  DerivativeType & derivative ) const
{
  unsigned long numberOfDerivatives = this->m_FixedPointSet->GetNumberOfPoints();
  if ( this->m_UseWithRespectToTheMovingPointSet )
    {
    numberOfDerivatives = this->m_MovingPointSet->GetNumberOfPoints();
    }
  derivative.SetSize( numberOfDerivatives, PointDimension );
  derivative.Fill( 0 );

  typename LabelSetType::const_iterator iter;

  for ( iter = this->m_FixedLabelSet.begin();
        iter != this->m_FixedLabelSet.end(); ++iter )
    {
    PixelType currentLabel = *iter;

    /**
     * check to see if the moving label set contains the same label
     */
    if ( find( this->m_MovingLabelSet.begin(), this->m_MovingLabelSet.end(),
           currentLabel ) == this->m_MovingLabelSet.end() )
      {
      continue;
      }

    /**
     * Collect all the fixed and moving points with the currentLabel
     */
    typename PointSetType::Pointer fixedLabelPoints
      = PointSetType::New();
    fixedLabelPoints->Initialize();
    unsigned long fixedCount = 0;

    std::vector<long> fixedIndices;
    fixedIndices.clear();

    typename PointSetType::PointsContainerConstIterator ItF =
      this->m_FixedPointSet->GetPoints()->Begin();
    typename PointSetType::PointDataContainerIterator ItFD =
      this->m_FixedPointSet->GetPointData()->Begin();

    while ( ItF != this->m_FixedPointSet->GetPoints()->End() )
      {
      if ( ItFD.Value() == currentLabel )
        {
        fixedLabelPoints->SetPoint( fixedCount++, ItF.Value() );
        fixedIndices.push_back( ItF.Index() );
        }
      ++ItF;
      ++ItFD;
      }

    typename PointSetType::Pointer movingLabelPoints
      = PointSetType::New();
    movingLabelPoints->Initialize();
    unsigned long movingCount = 0;

    std::vector<long> movingIndices;
    movingIndices.clear();

    typename PointSetType::PointsContainerConstIterator ItM =
      this->m_MovingPointSet->GetPoints()->Begin();
    typename PointSetType::PointDataContainerIterator ItMD =
      this->m_MovingPointSet->GetPointData()->Begin();

    while ( ItM != this->m_MovingPointSet->GetPoints()->End() )
      {
      if ( ItMD.Value() == currentLabel )
        {
        movingLabelPoints->SetPoint( movingCount++, ItM.Value() );
        movingIndices.push_back( ItM.Index() );
        }
      ++ItM;
      ++ItMD;
      }

    /**
     * Invoke the single label JensenTsallis measure
     */
    typedef EuclideanDistancePointSetMetric<PointSetType> MetricType;
    typename MetricType::Pointer metric = MetricType::New();

    metric->SetFixedPointSet( fixedLabelPoints );
    metric->SetMovingPointSet( movingLabelPoints );
    metric->SetUseWithRespectToTheMovingPointSet(
      this->m_UseWithRespectToTheMovingPointSet );

    metric->Initialize();

    DerivativeType labelDerivative;
    metric->GetDerivative( parameters, labelDerivative );

    RealType avgNorm = 0.0;
    for ( unsigned int i = 0; i < metric->GetNumberOfValues(); i++ )
      {
      RealType norm = 0.0;
      for ( unsigned int j = 0; j < PointDimension; j++ )
        {
        norm += ( labelDerivative(i, j) * labelDerivative(i, j) );
        }
      avgNorm += vcl_sqrt( norm );
      }
    avgNorm /= static_cast<RealType>( metric->GetNumberOfValues() );
    labelDerivative /= avgNorm;

    if ( this->m_UseWithRespectToTheMovingPointSet )
      {
      std::vector<long>::const_iterator it;
      unsigned long index = 0;
      for ( it = movingIndices.begin(); it != movingIndices.end(); ++it )
        {
        for ( unsigned int d = 0; d < PointDimension; d++ )
          {
          derivative( *it, d ) = labelDerivative( index, d );
          }
        index++;
        }
      }
    else
      {
      std::vector<long>::const_iterator it;
      unsigned long index = 0;
      for ( it = fixedIndices.begin(); it != fixedIndices.end(); ++it )
        {
        for ( unsigned int d = 0; d < PointDimension; d++ )
          {
          derivative( *it, d ) = labelDerivative( index, d );
          }
        index++;
        }
      }
    }
}

/** Get both the match Measure and theDerivative Measure  */
template <class TPointSet>
void
EuclideanDistanceLabeledPointSetMetric<TPointSet>
::GetValueAndDerivative( const TransformParametersType & parameters,
  MeasureType & value, DerivativeType  & derivative ) const
{
  unsigned long numberOfDerivatives = this->m_FixedPointSet->GetNumberOfPoints();
  if ( this->m_UseWithRespectToTheMovingPointSet )
    {
    numberOfDerivatives = this->m_MovingPointSet->GetNumberOfPoints();
    }
  derivative.SetSize( numberOfDerivatives, PointDimension );
  derivative.Fill( 0 );

  value.SetSize( 1 );
  value.Fill( 0 );

  typename LabelSetType::const_iterator iter;

  for ( iter = this->m_FixedLabelSet.begin();
        iter != this->m_FixedLabelSet.end(); ++iter )
    {
    PixelType currentLabel = *iter;

    /**
     * check to see if the moving label set contains the same label
     */
    if ( find( this->m_MovingLabelSet.begin(), this->m_MovingLabelSet.end(),
           currentLabel ) == this->m_MovingLabelSet.end() )
      {
      continue;
      }

    /**
     * Collect all the fixed and moving points with the currentLabel
     */
    typename PointSetType::Pointer fixedLabelPoints
      = PointSetType::New();
    fixedLabelPoints->Initialize();
    unsigned long fixedCount = 0;

    std::vector<long> fixedIndices;
    fixedIndices.clear();

    typename PointSetType::PointsContainerConstIterator ItF =
      this->m_FixedPointSet->GetPoints()->Begin();
    typename PointSetType::PointDataContainerIterator ItFD =
      this->m_FixedPointSet->GetPointData()->Begin();

    while ( ItF != this->m_FixedPointSet->GetPoints()->End() )
      {
      if ( ItFD.Value() == currentLabel )
        {
        fixedLabelPoints->SetPoint( fixedCount++, ItF.Value() );
        fixedIndices.push_back( ItF.Index() );
        }
      ++ItF;
      ++ItFD;
      }

    typename PointSetType::Pointer movingLabelPoints
      = PointSetType::New();
    movingLabelPoints->Initialize();
    unsigned long movingCount = 0;

    std::vector<long> movingIndices;
    movingIndices.clear();

    typename PointSetType::PointsContainerConstIterator ItM =
      this->m_MovingPointSet->GetPoints()->Begin();
    typename PointSetType::PointDataContainerIterator ItMD =
      this->m_MovingPointSet->GetPointData()->Begin();

    while ( ItM != this->m_MovingPointSet->GetPoints()->End() )
      {
      if ( ItMD.Value() == currentLabel )
        {
        movingLabelPoints->SetPoint( movingCount++, ItM.Value() );
        movingIndices.push_back( ItM.Index() );
        }
      ++ItM;
      ++ItMD;
      }

    /**
     * Invoke the single label Euclidean measure
     */
    typedef EuclideanDistancePointSetMetric<PointSetType> MetricType;
    typename MetricType::Pointer metric = MetricType::New();

    metric->SetFixedPointSet( fixedLabelPoints );
    metric->SetMovingPointSet( movingLabelPoints );
    metric->SetUseWithRespectToTheMovingPointSet(
      this->m_UseWithRespectToTheMovingPointSet );

    metric->Initialize();

    DerivativeType labelDerivative;
    MeasureType labelValue;

    metric->GetValueAndDerivative( parameters, labelValue, labelDerivative );

    value[0] += labelValue[0];

    RealType avgNorm = 0.0;
    for ( unsigned int i = 0; i < metric->GetNumberOfValues(); i++ )
      {
      RealType norm = 0.0;
      for ( unsigned int j = 0; j < PointDimension; j++ )
        {
        norm += ( labelDerivative(i, j) * labelDerivative(i, j) );
        }
      avgNorm += vcl_sqrt( norm );
      }
    avgNorm /= static_cast<RealType>( metric->GetNumberOfValues() );
    labelDerivative /= avgNorm;

    if ( this->m_UseWithRespectToTheMovingPointSet )
      {
      std::vector<long>::const_iterator it;
      unsigned long index = 0;
      for ( it = movingIndices.begin(); it != movingIndices.end(); ++it )
        {
        for ( unsigned int d = 0; d < PointDimension; d++ )
          {
          derivative( *it, d ) = labelDerivative( index, d );
          }
        index++;
        }
      }
    else
      {
      std::vector<long>::const_iterator it;
      unsigned long index = 0;
      for ( it = fixedIndices.begin(); it != fixedIndices.end(); ++it )
        {
        for ( unsigned int d = 0; d < PointDimension; d++ )
          {
          derivative( *it, d ) = labelDerivative( index, d );
          }
        index++;
        }
      }
    }
}

template <class TPointSet>
void
EuclideanDistanceLabeledPointSetMetric<TPointSet>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Use with respect to the moving point set: "
     << this->m_UseWithRespectToTheMovingPointSet << std::endl;
}

} // end namespace itk


#endif

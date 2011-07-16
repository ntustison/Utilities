/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkEuclideanDistancePointSetMetric.hxx,v $
  Language:  C++
  Date:      $Date: 2008/12/06 05:19:54 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkEuclideanDistancePointSetMetric_hxx
#define __itkEuclideanDistancePointSetMetric_hxx

#include "itkEuclideanDistancePointSetMetric.h"

namespace itk {

template <class TPointSet>
EuclideanDistancePointSetMetric<TPointSet>
::EuclideanDistancePointSetMetric()
{
  this->m_UseWithRespectToTheMovingPointSet = true;

  this->m_KdTreeGenerator = NULL;
  this->m_SamplePoints = NULL;

  typename DefaultTransformType::Pointer transform
    = DefaultTransformType::New();
  transform->SetIdentity();

  Superclass::SetTransform( transform );
}

/** Initialize the metric */
template <class TPointSet>
void
EuclideanDistancePointSetMetric<TPointSet>
::Initialize( void ) throw ( ExceptionObject )
{
  Superclass::Initialize();

  /**
   * Generate KdTrees for the opposite point set
   */
  this->m_SamplePoints = SampleType::New();
  this->m_SamplePoints->SetMeasurementVectorSize( PointDimension );

  if( this->m_UseWithRespectToTheMovingPointSet )
    {
    typename PointSetType::PointsContainerConstIterator It
      = this->m_FixedPointSet->GetPoints()->Begin();
    while( It != this->m_FixedPointSet->GetPoints()->End() )
      {
      MeasurementVectorType mv;
      PointType point = It.Value();
      for( unsigned int d = 0; d < PointDimension; d++ )
        {
        mv[d] = point[d];
        }
      this->m_SamplePoints->PushBack( mv );
      ++It;
      }
    }
  else
    {
    typename PointSetType::PointsContainerConstIterator It
      = this->m_MovingPointSet->GetPoints()->Begin();
    while( It != this->m_MovingPointSet->GetPoints()->End() )
      {
      MeasurementVectorType mv;
      PointType point = It.Value();
      for( unsigned int d = 0; d < PointDimension; d++ )
        {
        mv[d] = point[d];
        }
      this->m_SamplePoints->PushBack( mv );
      ++It;
      }
    }

  this->m_KdTreeGenerator = TreeGeneratorType::New();
  this->m_KdTreeGenerator->SetSample( this->m_SamplePoints );
  this->m_KdTreeGenerator->SetBucketSize( 4 );
  this->m_KdTreeGenerator->Update();
}

/** Return the number of values, i.e the number of points in the moving set */
template <class TPointSet>
unsigned int
EuclideanDistancePointSetMetric<TPointSet>
::GetNumberOfValues() const
{
  if( this->m_UseWithRespectToTheMovingPointSet )
    {
    if( this->m_MovingPointSet )
      {
      return this->m_MovingPointSet->GetPoints()->Size();
      }
    }
  else
    {
    if( this->m_FixedPointSet )
      {
      return this->m_FixedPointSet->GetPoints()->Size();
      }
    }

  return 0;
}


/** Get the match Measure */
template <class TPointSet>
typename EuclideanDistancePointSetMetric
  <TPointSet>::MeasureType
EuclideanDistancePointSetMetric<TPointSet>
::GetValue( const TransformParametersType & parameters ) const
{
  /**
   * Only identity transform is valid
   */
  //  this->SetTransformParameters( parameters );

  PointSetPointer points[2];

  if( this->m_UseWithRespectToTheMovingPointSet )
    {
    points[0] = const_cast<PointSetType *>(
      static_cast<const PointSetType *>( this->m_FixedPointSet ) );
    points[1] = const_cast<PointSetType *>(
      static_cast<const PointSetType *>( this->m_MovingPointSet ) );
    }
  else
    {
    points[1] = const_cast<PointSetType *>(
      static_cast<const PointSetType *>( this->m_FixedPointSet ) );
    points[0] = const_cast<PointSetType *>(
      static_cast<const PointSetType *>( this->m_MovingPointSet ) );
    }

  MeasureType measure;
  measure.SetSize( 1 );
  measure.Fill( 0 );

  typename PointSetType::PointsContainerConstIterator It
    = points[1]->GetPoints()->Begin();
  while( It != points[1]->GetPoints()->End() )
    {
    PointType point = It.Value();

    MeasurementVectorType queryPoint;
    for( unsigned int d = 0; d < PointDimension; d++ )
      {
      queryPoint[d] = point[d];
      }
    typename TreeGeneratorType::KdTreeType
      ::InstanceIdentifierVectorType neighbors;
    this->m_KdTreeGenerator->GetOutput()->Search( queryPoint, 1u, neighbors );

    MeasurementVectorType neighbor =
      this->m_KdTreeGenerator->GetOutput()->GetMeasurementVector( neighbors[0] );

    RealType sum = 0.0;
    for( unsigned int d = 0; d < PointDimension; d++ )
      {
      sum += vnl_math_sqr( neighbor[d] - queryPoint[d] );
      }
    measure[0] += vcl_sqrt( sum );

    ++It;
    }

  measure[0] /= static_cast<RealType>( points[1]->GetNumberOfPoints() );
  return measure;
}

/** Get the Derivative Measure */
template <class TPointSet>
void
EuclideanDistancePointSetMetric<TPointSet>
::GetDerivative( const TransformParametersType & parameters,
                 DerivativeType & derivative ) const
{
  /**
   * Only identity transform is valid
   */
  //  this->SetTransformParameters( parameters );

  PointSetPointer points[2];

  if( this->m_UseWithRespectToTheMovingPointSet )
    {
    points[0] = const_cast<PointSetType *>(
      static_cast<const PointSetType *>( this->m_FixedPointSet ) );
    points[1] = const_cast<PointSetType *>(
      static_cast<const PointSetType *>( this->m_MovingPointSet ) );
    }
  else
    {
    points[1] = const_cast<PointSetType *>(
      static_cast<const PointSetType *>( this->m_FixedPointSet ) );
    points[0] = const_cast<PointSetType *>(
      static_cast<const PointSetType *>( this->m_MovingPointSet ) );
    }
  derivative.SetSize( points[1]->GetPoints()->Size(), PointDimension );
  derivative.Fill( 0 );

  unsigned long count = 0;
  typename PointSetType::PointsContainerConstIterator It
    = points[1]->GetPoints()->Begin();
  while( It != points[1]->GetPoints()->End() )
    {
    PointType point = It.Value();

    MeasurementVectorType queryPoint;
    for( unsigned int d = 0; d < PointDimension; d++ )
      {
      queryPoint[d] = point[d];
      }
    typename TreeGeneratorType::KdTreeType
      ::InstanceIdentifierVectorType neighbors;
    this->m_KdTreeGenerator->GetOutput()->Search( queryPoint, 1u, neighbors );

    MeasurementVectorType neighbor =
      this->m_KdTreeGenerator->GetOutput()->GetMeasurementVector( neighbors[0] );

    for( unsigned int d = 0; d < PointDimension; d++ )
      {
      derivative( count, d ) = -( neighbor[d] - queryPoint[d] );
      }
    count++;
    ++It;
    }
}

/** Get both the match Measure and theDerivative Measure  */
template <class TPointSet>
void
EuclideanDistancePointSetMetric<TPointSet>
::GetValueAndDerivative( const TransformParametersType & parameters,
  MeasureType & value, DerivativeType  & derivative ) const
{
  /**
   * Only identity transform is valid
   */
  //  this->SetTransformParameters( parameters );

  PointSetPointer points[2];

  if( this->m_UseWithRespectToTheMovingPointSet )
    {
    points[0] = const_cast<PointSetType *>(
      static_cast<const PointSetType *>( this->m_FixedPointSet ) );
    points[1] = const_cast<PointSetType *>(
      static_cast<const PointSetType *>( this->m_MovingPointSet ) );
    }
  else
    {
    points[1] = const_cast<PointSetType *>(
      static_cast<const PointSetType *>( this->m_FixedPointSet ) );
    points[0] = const_cast<PointSetType *>(
      static_cast<const PointSetType *>( this->m_MovingPointSet ) );
    }
  derivative.SetSize( points[1]->GetPoints()->Size(), PointDimension );
  derivative.Fill( 0.0 );

  value.SetSize( 1 );
  value.Fill( 0.0 );

  unsigned long count = 0;
  typename PointSetType::PointsContainerConstIterator It
    = points[1]->GetPoints()->Begin();
  while( It != points[1]->GetPoints()->End() )
    {
    PointType point = It.Value();

    MeasurementVectorType queryPoint;
    for( unsigned int d = 0; d < PointDimension; d++ )
      {
      queryPoint[d] = point[d];
      }
    typename TreeGeneratorType::KdTreeType
      ::InstanceIdentifierVectorType neighbors;
    this->m_KdTreeGenerator->GetOutput()->Search( queryPoint, 1u, neighbors );

    MeasurementVectorType neighbor =
      this->m_KdTreeGenerator->GetOutput()->GetMeasurementVector( neighbors[0] );

    RealType sum = 0.0;
    for( unsigned int d = 0; d < PointDimension; d++ )
      {
      derivative( count, d ) = -( neighbor[d] - queryPoint[d] );
      sum += vnl_math_sqr( neighbor[d] - queryPoint[d] );
      }
    value[0] += vcl_sqrt( sum );

    count++;
    ++It;
    }

  value[0] /= static_cast<RealType>( points[1]->GetNumberOfPoints() );
}

template <class TPointSet>
void
EuclideanDistancePointSetMetric<TPointSet>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Use with respect to the moving point set: "
     << this->m_UseWithRespectToTheMovingPointSet << std::endl;
}

} // end namespace itk


#endif

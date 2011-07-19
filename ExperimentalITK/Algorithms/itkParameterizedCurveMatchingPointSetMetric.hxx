/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkParameterizedCurveMatchingPointSetMetric_hxx
#define __itkParameterizedCurveMatchingPointSetMetric_hxx

#include "itkParameterizedCurveMatchingPointSetMetric.h"

#include "itkBSplineControlPointImageFilter.h"
#include "itkBSplineControlPointImageFunction.h"
#include "itkVariableSizeMatrix.h"

namespace itk {

template<class TPointSet>
ParameterizedCurveMatchingPointSetMetric<TPointSet>
::ParameterizedCurveMatchingPointSetMetric() :
  m_SplineOrder( 3 ),
  m_NumberOfFittingLevels( 4 ),
  m_SamplingRate( 0.0001 ),
  m_UseClosedCurves( false ),
  m_ParametricWindow( 0.1 )
{

  this->m_FixedBSplineCurveControlPointImage = NULL;
  this->m_MovingBSplineCurveControlPointImage = NULL;
}

/** Initialize the metric */
template<class TPointSet>
void
ParameterizedCurveMatchingPointSetMetric<TPointSet>
::Initialize( void ) throw ( ExceptionObject )
{
  if ( !this->m_FixedTransform )
    {
    itkExceptionMacro( "Fixed transform is not present" );
    }

  if ( !this->m_MovingTransform )
    {
    itkExceptionMacro( "Moving transform is not present" );
    }

  if ( !this->m_FixedPointSet )
    {
    itkExceptionMacro( "Fixed point set is not present" );
    }

  if ( !this->m_MovingPointSet )
    {
    itkExceptionMacro( "Moving point set is not present" );
    }

  // If the PointSet is provided by a source, update the source.
  if( this->m_MovingPointSet->GetSource() )
    {
    this->m_MovingPointSet->GetSource()->Update();
    }
  this->TransformMovingPointSet();

  // If the point set is provided by a source, update the source.
  if( this->m_FixedPointSet->GetSource() )
    {
    this->m_FixedPointSet->GetSource()->Update();
    }
  this->TransformFixedPointSet();

  if( this->m_UseClosedCurves )
    {
    this->InitializePointsLocators();
    }

  // Construct the fixed curve
  this->m_FixedBSplineCurveControlPointImage =
    this->GenerateBSplineCurveControlPoints(
    this->m_FixedTransformedPointSet.GetPointer() );

  // Construct the moving curve
  this->m_MovingBSplineCurveControlPointImage =
    this->GenerateBSplineCurveControlPoints(
    this->m_MovingTransformedPointSet.GetPointer() );
}

/** Get the match Measure */
template<class TPointSet>
typename ParameterizedCurveMatchingPointSetMetric<TPointSet>
::MeasureType
ParameterizedCurveMatchingPointSetMetric<TPointSet>
::GetValue() const
{
  MeasureType value;
  DerivativeType derivative;
  this->GetValueAndDerivative( value, derivative );
  return value;
}

/** Get the Derivative Measure */
template<class TPointSet>
void
ParameterizedCurveMatchingPointSetMetric<TPointSet>
::GetDerivative( DerivativeType &derivative ) const
{
  MeasureType value;
  this->GetValueAndDerivative( value, derivative );
}

/** Get both the match Measure and the Derivative Measure  */
template<class TPointSet>
void
ParameterizedCurveMatchingPointSetMetric<TPointSet>
::GetValueAndDerivative( MeasureType &value, DerivativeType &derivative ) const
{
  CurveImagePointer fixedCurveImage = this->ReconstructBSplineCurve(
    this->m_FixedBSplineCurveControlPointImage.GetPointer() );

  CurveImagePointer movingCurveImage = this->ReconstructBSplineCurve(
    this->m_MovingBSplineCurveControlPointImage.GetPointer() );

  typedef BSplineControlPointImageFunction<CurveImageType> BSplineFunctionType;
  typename BSplineFunctionType::ArrayType isClosed( this->m_UseClosedCurves );
  typename CurveImageType::PointType origin( 0.0 );
  typename CurveImageType::SpacingType spacing( this->m_SamplingRate );
  typename CurveImageType::SizeType size(
    static_cast<unsigned int>( 1.0/spacing[0] ) + 1 );

  typename BSplineFunctionType::Pointer fixedBSplineFunction =
    BSplineFunctionType::New();
  fixedBSplineFunction->SetCloseDimension( isClosed );
  fixedBSplineFunction->SetSplineOrder( this->m_SplineOrder );
  fixedBSplineFunction->SetOrigin( origin );
  fixedBSplineFunction->SetSpacing( spacing );
  fixedBSplineFunction->SetSize( size );
  fixedBSplineFunction->SetInputImage(
    this->m_FixedBSplineCurveControlPointImage );

  typename BSplineFunctionType::Pointer movingBSplineFunction =
    BSplineFunctionType::New();
  movingBSplineFunction->SetCloseDimension( isClosed );
  movingBSplineFunction->SetSplineOrder( this->m_SplineOrder );
  movingBSplineFunction->SetOrigin( origin );
  movingBSplineFunction->SetSpacing( spacing );
  movingBSplineFunction->SetSize( size );
  movingBSplineFunction->SetInputImage(
    this->m_FixedBSplineCurveControlPointImage );

  VariableSizeMatrix<MeasureType> costMatrix(
    fixedCurveImage->GetLargestPossibleRegion().GetSize()[0],
    movingCurveImage->GetLargestPossibleRegion().GetSize()[0] );

}

template<class TPointSet>
typename ParameterizedCurveMatchingPointSetMetric<TPointSet>::CurveImagePointer
ParameterizedCurveMatchingPointSetMetric<TPointSet>
::GenerateBSplineCurveControlPoints( PointSetType *sampledPoints )
{
  RealType totalDistance = 0.0;

  PointType pointEnd;
  pointEnd.Fill( 0.0 );

  PointsContainerConstIterator It = sampledPoints->GetPoints()->Begin();
  PointType pointBegin = It.Value();
  ++It;
  while( It != sampledPoints->GetPoints()->End() )
    {
    pointEnd = It.Value();
    totalDistance += pointBegin.EuclideanDistanceTo( pointEnd );
    pointBegin = pointEnd;
    ++It;
    }

  typename SampledCurveType::Pointer curvePoints = SampledCurveType::New();
  curvePoints->Initialize();

  RealType distance = 0.0;

  pointEnd.Fill( 0.0 );
  typename SampledCurveType::PointType parametricValue;
  curvePoints->SetPoint( 0, parametricValue );
  curvePoints->SetPointData( 0, pointBegin.GetVectorFromOrigin() );

  unsigned int count = 1;

  It = sampledPoints->GetPoints()->Begin();
  pointBegin = It.Value();
  ++It;
  while( It != sampledPoints->GetPoints()->End() )
    {
    pointEnd = It.Value();
    distance += pointBegin.EuclideanDistanceTo( pointEnd );
    pointBegin = pointEnd;

    parametricValue[0] = distance / totalDistance;
    curvePoints->SetPoint( count, parametricValue );
    curvePoints->SetPointData( count, pointBegin.GetVectorFromOrigin() );

    ++count;
    ++It;
    }

  typename BSplineFilterType::ArrayType numberOfControlPoints(
    this->m_SplineOrder + 1 );
  typename BSplineFilterType::ArrayType isClosed( this->m_UseClosedCurves );
  typename CurveImageType::PointType origin( 0.0 );
  typename CurveImageType::SpacingType spacing( this->m_SamplingRate );
  typename CurveImageType::SizeType size(
    static_cast<unsigned int>( 1.0/spacing[0] ) + 1 );

  typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
  bspliner->SetInput( curvePoints );
  bspliner->SetNumberOfControlPoints( numberOfControlPoints );
  bspliner->SetCloseDimension( isClosed );
  bspliner->SetNumberOfLevels( this->m_NumberOfLevels );
  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetOrigin( origin );
  bspliner->SetSpacing( spacing );
  bspliner->SetSize( size );
  bspliner->GenerateOutputImage( false );
  bspliner->Update();

  CurveImageType::Pointer controlPointImage = bspliner->GetControlPointPhiLattice();
  controlPointImage->DisconnectPipeline();

  return controlPointImage;
}

template<class TPointSet>
typename ParameterizedCurveMatchingPointSetMetric<TPointSet>::CurveImagePointer
ParameterizedCurveMatchingPointSetMetric<TPointSet>
::ReconstructBSplineCurve( CurveImageType *controlPoints )
{
  typedef BSplineControlPointImageFilter<CurveImageType> ReconstructorType;
  typename ReconstructorType::Pointer bspliner = ReconstructorType::New();

  typename ReconstructorType::ArrayType isClosed( this->m_UseClosedCurves );
  typename CurveImageType::PointType origin( 0.0 );
  typename CurveImageType::SpacingType spacing( this->m_SamplingRate );
  typename CurveImageType::SizeType size(
    static_cast<unsigned int>( 1.0/spacing[0] ) + 1 );

  bspliner->SetInput( curvePoints );
  bspliner->SetNumberOfControlPoints( numberOfControlPoints );
  bspliner->SetCloseDimension( isClosed );
  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetOrigin( origin );
  bspliner->SetSpacing( spacing );
  bspliner->SetSize( size );

  CurveImageType::Pointer curveImage = bspliner->GetOutput();
  curveImage->Update();
  curveImage->DisconnectPipeline();

  return curveImage;
}

template<class TPointSet>
void
ParameterizedCurveMatchingPointSetMetric<TPointSet>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Use regularization term: "
     << this->m_UseRegularizationTerm << std::endl;
  os << indent << "Alpha: "
     << this->m_Alpha << std::endl;

  os << indent << "Point set sigma: "
     << this->m_PointSetSigma << std::endl;

  if( this->m_UseAnisotropicCovariances )
    {
    os << indent << "Kernel sigma: "
       << this->m_KernelSigma << std::endl;
    os << indent << "Covariance k-neighborhood: "
       << this->m_CovarianceKNeighborhood << std::endl;
    }
  else
    {
    os << indent << "Isotropic covariances are used." << std::endl;
    }

}
} // end namespace itk


#endif

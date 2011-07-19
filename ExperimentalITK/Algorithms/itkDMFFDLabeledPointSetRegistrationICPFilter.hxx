/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: itkDMFFDLabeledPointSetRegistrationICPFilter.hxx,v $
Language:  C++

Date:      $Date: 2009/06/02 17:32:46 $
Version:   $Revision: 1.6 $

Copyright (c) Insight Software Consortium. All rights reser
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for detail.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkDMFFDLabeledPointSetRegistrationICPFilter_hxx
#define __itkDMFFDLabeledPointSetRegistrationICPFilter_hxx

#include "itkDMFFDLabeledPointSetRegistrationICPFilter.h"

#include "itkAddImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkEuclideanDistanceLabeledPointSetMetric.h"
#include "itkMultiplyImageFilter.h"

#include "vnl/vnl_math.h"

namespace itk {

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
DMFFDLabeledPointSetRegistrationICPFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::DMFFDLabeledPointSetRegistrationICPFilter()
{
  this->SetNumberOfRequiredInputs( 2 );

  // Default values
  this->m_SplineOrder = 3;
  this->m_Size.Fill( 0 );
  this->m_Direction.SetIdentity();
  this->m_Directionality.Fill( 1 );
  this->m_CalculateInitialSimilarityTransform = true;

  this->m_MaximumNumberOfIterations.SetSize( 3 );
  this->m_MaximumNumberOfIterations.Fill( 10 );
  this->m_GradientScalingFactor.Fill( NumericTraits<RealType>::One );
  this->m_LineSearchMaximumIterations = 0;
  this->m_LineSearchMaximumStepSize = NumericTraits<RealType>::max();

  this->m_Verbose = true;

  this->m_InitialDeformationFieldControlPoints
    = ControlPointLatticeType::New();
  this->m_InitialDeformationFieldControlPoints = NULL;
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
DMFFDLabeledPointSetRegistrationICPFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::~DMFFDLabeledPointSetRegistrationICPFilter()
{
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
void
DMFFDLabeledPointSetRegistrationICPFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::GenerateData()
{
  for( this->m_CurrentLevel = 0; this->m_CurrentLevel
    < this->m_MaximumNumberOfIterations.Size(); this->m_CurrentLevel++ )
    {
    this->Initialize();

    if( this->m_Verbose )
      {
      std::cout << "Current level = " << this->m_CurrentLevel+1 << " ("
        << this->m_MaximumNumberOfIterations.Size() << " total levels).  "
        << "Control point grid size = "
        << this->m_TotalDeformationFieldControlPoints
          ->GetLargestPossibleRegion().GetSize() << ". " << std::endl;
      }

    itkDebugMacro( "Current level = " << this->m_CurrentLevel+1
      << " (" << this->m_MaximumNumberOfIterations.Size() << " total levels).  "
      << "Control point grid size = "
      << this->m_TotalDeformationFieldControlPoints
        ->GetLargestPossibleRegion().GetSize() << ". " );

    if( this->m_MaximumNumberOfIterations[this->m_CurrentLevel] > 0 )
      {
      itkDebugMacro( "Begin iterative solve. " );
      this->IterativeSolve();
      itkDebugMacro( "Finish iterative solve. " );
      }
    }

  // Warp the input points

  typename WarpedPointSetType::Pointer warpedPoints
    = WarpedPointSetType::New();
  warpedPoints->Initialize();

  typename ControlPointFilterType::Pointer bspliner
    = ControlPointFilterType::New();
  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetInput( this->m_TotalDeformationFieldControlPoints );
  bspliner->SetOrigin( this->m_Origin );
  bspliner->SetSize( this->m_Size );
  bspliner->SetSpacing( this->m_Spacing );
  bspliner->SetDirection( this->m_Direction );

  typename MovingPointSetType::PointsContainerConstIterator ItM =
    this->m_SimilarityTransformMovingPointSet->GetPoints()->Begin();
  typename MovingPointSetType::PointDataContainer::ConstIterator ItD =
    this->m_SimilarityTransformMovingPointSet->GetPointData()->Begin();
  while( ItM != this->m_SimilarityTransformMovingPointSet->GetPoints()->End() )
    {
    MovingPointType inputPoint = ItM.Value();

    typename ControlPointFilterType::PointType point;
    point.CastFrom( inputPoint );

    VectorType vector;
    vector.Fill( 0 );
    try
      {
      bspliner->EvaluateAtPoint( point, vector );
      }
    catch(...)
      {
      vector.Fill( 0 );
      }
    point += vector;
    for( unsigned int d = 0; d < Dimension; d++ )
      {
      inputPoint[d] = point[d];
      }

    warpedPoints->SetPoint( ItM.Index(), inputPoint );
    warpedPoints->SetPointData( ItM.Index(), ItD.Value() );
    ++ItM;
    ++ItD;
    }
  this->SetNthOutput( 0, warpedPoints );

}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
void
DMFFDLabeledPointSetRegistrationICPFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::Initialize()
{
  if( this->m_CurrentLevel == 0 )
    {
    this->CalculateInitialCentroidsAndNorms();

    if( this->m_Size[0] == 0 )
      {

      /**
       * Define the transformation domain based on the bounding box of the
       * two point-sets.
       */
      FixedPointType minPoint;
      minPoint.Fill( NumericTraits<RealType>::max() );
      FixedPointType maxPoint;
      maxPoint.Fill( NumericTraits<RealType>::min() );
      for( unsigned int i = 0; i < this->GetNumberOfInputs(); i++ )
        {
        for( unsigned int d = 0; d < Dimension; d++ )
          {
          minPoint[d] = vnl_math_min( minPoint[d],
            this->GetInput( i )->GetBoundingBox()->GetMinimum()[d] );
          maxPoint[d] = vnl_math_max( maxPoint[d],
            this->GetInput( i )->GetBoundingBox()->GetMaximum()[d] );
          }
        }

      this->m_Size.Fill( 100 );

      for( unsigned int d = 0; d < Dimension; d++ )
        {
        this->m_Origin[d] = minPoint[d];
        this->m_Spacing[d] = ( maxPoint[d] - minPoint[d] )
          / static_cast<RealType>( this->m_Size[d] - 1 );
        }
      }

    if( this->m_InitialDeformationFieldControlPoints )
      {
      for( unsigned int d = 0; d < Dimension; d++ )
        {
        this->m_InitialMeshResolution[d] =
          this->m_InitialDeformationFieldControlPoints->
            GetLargestPossibleRegion().GetSize()[d] - this->m_SplineOrder;

        if( this->m_InitialMeshResolution[d] < 1 )
          {
          itkExceptionMacro( "Invalid size for initial deformation field "
            << "control point lattice." );
          }
        }
      }
    else
      {
      typename ControlPointLatticeType::RegionType::SizeType size;
      for( unsigned int d = 0; d < Dimension; d++ )
        {
        size[d] = this->m_InitialMeshResolution[d] + this->m_SplineOrder;
        }

      this->m_InitialDeformationFieldControlPoints
        = ControlPointLatticeType::New();
      this->m_InitialDeformationFieldControlPoints->SetRegions( size );
      this->m_InitialDeformationFieldControlPoints->Allocate();

      VectorType V;
      V.Fill( 0.0 );
      this->m_InitialDeformationFieldControlPoints->Allocate();
      this->m_InitialDeformationFieldControlPoints->FillBuffer( V );
      }
    typedef ImageDuplicator<ControlPointLatticeType> DuplicatorType;
    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage( this->m_InitialDeformationFieldControlPoints );
    duplicator->Update();
    this->m_TotalDeformationFieldControlPoints = duplicator->GetOutput();

    VectorType V;
    V.Fill( 0 );
    this->m_GradientFieldControlPoints = ControlPointLatticeType::New();
    this->m_GradientFieldControlPoints->SetRegions(
      this->m_InitialDeformationFieldControlPoints->GetLargestPossibleRegion() );
    this->m_GradientFieldControlPoints->Allocate();
    this->m_GradientFieldControlPoints->FillBuffer( V );
    }
  else
    {

    typename ControlPointFilterType::Pointer bspliner
      = ControlPointFilterType::New();

    bspliner->SetSplineOrder( this->m_SplineOrder );
    bspliner->SetInput( this->m_TotalDeformationFieldControlPoints );

    typename BSplineFilterType::ArrayType nlevels;
    nlevels.Fill( 2 );
    this->m_TotalDeformationFieldControlPoints
      = bspliner->RefineControlPointLattice( nlevels );

    VectorType V;
    V.Fill( 0 );
    this->m_GradientFieldControlPoints = ControlPointLatticeType::New();
    this->m_GradientFieldControlPoints->SetRegions(
      this->m_TotalDeformationFieldControlPoints
        ->GetLargestPossibleRegion().GetSize() );
    this->m_GradientFieldControlPoints->Allocate();
    this->m_GradientFieldControlPoints->FillBuffer( V );
    }
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
void
DMFFDLabeledPointSetRegistrationICPFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::CalculateInitialCentroidsAndNorms()
{
  /**
   * Calculate the initial transform by finding the appropriate
   * translation and scale.
   */
  this->m_FixedPointSetCentroid.Fill( 0.0 );
  this->m_MovingPointSetCentroid.Fill( 0.0 );
  this->m_ScaleFactor = 1.0;

  bool calculateNonIdentitySimilarityTransform
    = this->m_CalculateInitialSimilarityTransform;
  if ( this->m_InitialDeformationFieldControlPoints )
    {
    calculateNonIdentitySimilarityTransform = false;
    }
  for( unsigned int d = 0; d < Dimension; d++ )
    {
    if( !this->m_Directionality[d] )
      {
      calculateNonIdentitySimilarityTransform = false;
      break;
      }
    }

  typename FixedPointSetType::PointsContainerConstIterator ItF =
    this->GetInput( 0 )->GetPoints()->Begin();
  while( ItF != this->GetInput( 0 )->GetPoints()->End() )
    {
    for( unsigned int d = 0; d < Dimension; d++ )
      {
      this->m_FixedPointSetCentroid[d] += ItF.Value()[d];
      }
    ++ItF;
    }
  for( unsigned int d = 0; d < Dimension; d++ )
    {
    this->m_FixedPointSetCentroid[d] /= static_cast<RealType>(
      this->GetInput( 0 )->GetNumberOfPoints() );
    }
  RealType fixedPointSetNorm = 0.0;
  ItF = this->GetInput( 0 )->GetPoints()->Begin();
  while( ItF != this->GetInput( 0 )->GetPoints()->End() )
    {
    fixedPointSetNorm +=
      ( ItF.Value() - this->m_FixedPointSetCentroid ).GetSquaredNorm();
    ++ItF;
    }
  fixedPointSetNorm = vcl_sqrt( fixedPointSetNorm
    / static_cast<RealType>( this->GetInput( 0 )->GetNumberOfPoints() ) );

  typename MovingPointSetType::PointsContainerConstIterator ItM =
    this->GetInput( 1 )->GetPoints()->Begin();
  while( ItM != this->GetInput( 1 )->GetPoints()->End() )
    {
    for( unsigned int d = 0; d < Dimension; d++ )
      {
      this->m_MovingPointSetCentroid[d] += ItM.Value()[d];
      }
    ++ItM;
    }
  for( unsigned int d = 0; d < Dimension; d++ )
    {
    this->m_MovingPointSetCentroid[d] /= static_cast<RealType>(
      this->GetInput( 1 )->GetNumberOfPoints() );
    }
  RealType movingPointSetNorm = 0.0;
  ItM = this->GetInput( 1 )->GetPoints()->Begin();
  while( ItM != this->GetInput( 1 )->GetPoints()->End() )
    {
    movingPointSetNorm +=
      ( ItM.Value() - this->m_MovingPointSetCentroid ).GetSquaredNorm();
    ++ItM;
    }
  movingPointSetNorm = vcl_sqrt( movingPointSetNorm
    / static_cast<RealType>( this->GetInput( 1 )->GetNumberOfPoints() ) );

  this->m_ScaleFactor = fixedPointSetNorm / movingPointSetNorm;

  if( !calculateNonIdentitySimilarityTransform )
    {
    this->m_FixedPointSetCentroid.Fill( 0.0 );
    this->m_MovingPointSetCentroid.Fill( 0.0 );
    this->m_ScaleFactor = 1.0;
    }


  /**
   * Create a new moving point set that has been shifted and scaled
   */
  this->m_SimilarityTransformMovingPointSet = MovingPointSetType::New();
  this->m_SimilarityTransformMovingPointSet->Initialize();

  typename MovingPointSetType::PointsContainerConstIterator ItP =
    this->GetInput( 1 )->GetPoints()->Begin();
  typename MovingPointSetType::PointDataContainer::ConstIterator ItD =
    this->GetInput( 1 )->GetPointData()->Begin();

  while( ItP != this->GetInput( 1 )->GetPoints()->End() )
    {
    MovingPointType inputPoint = ItP.Value();

    inputPoint -= ( this->m_MovingPointSetCentroid.GetVectorFromOrigin() );
    for( unsigned int d = 0; d < Dimension; d++ )
      {
      inputPoint[d] *= this->m_ScaleFactor;
      }
    inputPoint += ( this->m_FixedPointSetCentroid.GetVectorFromOrigin() );

    this->m_SimilarityTransformMovingPointSet->SetPoint(
      ItP.Index(), inputPoint );
    this->m_SimilarityTransformMovingPointSet->SetPointData(
      ItP.Index(), ItD.Value() );

    ++ItP;
    ++ItD;
    }
}


template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
void
DMFFDLabeledPointSetRegistrationICPFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::IterativeSolve()
{
  const RealType ftol = 1.0e-50;
  const RealType eps = 1.0e-50;

  this->m_CurrentIteration = -1;

  RealType fp = this->EvaluateMetricAndGradient();

  ControlPointLatticePointer G = ControlPointLatticeType::New();
  ControlPointLatticePointer H = ControlPointLatticeType::New();

  typedef ImageDuplicator<ControlPointLatticeType> DuplicatorType;
  typename DuplicatorType::Pointer duplicatorH = DuplicatorType::New();
  duplicatorH->SetInputImage( this->m_GradientFieldControlPoints );
  duplicatorH->Update();
  H = duplicatorH->GetOutput();

  typename DuplicatorType::Pointer duplicatorG = DuplicatorType::New();
  duplicatorG->SetInputImage( this->m_GradientFieldControlPoints );
  duplicatorG->Update();
  G = duplicatorG->GetOutput();

  ImageRegionIterator<ControlPointLatticeType>
    ItG( G, G->GetLargestPossibleRegion() );
  ImageRegionIterator<ControlPointLatticeType>
    ItH( H, H->GetLargestPossibleRegion() );
  ImageRegionIterator<ControlPointLatticeType>
    ItXi( this->m_GradientFieldControlPoints,
    this->m_GradientFieldControlPoints->GetLargestPossibleRegion() );

  ItG.GoToBegin();
  ItH.GoToBegin();
  ItXi.GoToBegin();
  while( !ItG.IsAtEnd() )
    {
    ItG.Set( -ItG.Get() );
    ItH.Set( -ItH.Get() );
    ItXi.Set( -ItXi.Get() );
    ++ItG;
    ++ItH;
    ++ItXi;
    }

  for( unsigned its = 1;
    its <= this->m_MaximumNumberOfIterations[this->m_CurrentLevel]; its++ )
    {
    this->m_CurrentIteration = its-1;

    if( this->m_Verbose )
      {
      std::cout << "  Iteration " << its
        << " (of " << this->m_MaximumNumberOfIterations[this->m_CurrentLevel]
        << "): Energy = " << fp << std::endl;
      }
    itkDebugMacro( "  Iteration " << its
      << " (of " << this->m_MaximumNumberOfIterations[this->m_CurrentLevel]
      << "): Energy = " << fp );

    RealType gradientStep = 1.0;
    RealType fret;
    if( this->m_LineSearchMaximumIterations > 0 )
      {
      this->LineMinimization( &gradientStep, &fret );
      if( gradientStep == 0.0 )
        {
        std::cout << "  Exit condition (gradientStep == 0)" << std::endl;
        return;
        }
      }
    else
      {
      fret = this->EvaluateMetric( gradientStep );
      }

    if( 2.0*vnl_math_abs( fret - fp )
      <= ftol*( vnl_math_abs( fret ) + vnl_math_abs( fp ) + eps ) )
      {
      std::cout << "  Exit condition (normal):  |fret - fp| = "
        << vnl_math_abs( fret - fp ) << std::endl;
      return;
      }
    fp = fret;

    typename RealImageType::Pointer lambda = RealImageType::New();
    lambda->SetOrigin( this->m_GradientFieldControlPoints->GetOrigin() );
    lambda->SetSpacing( this->m_GradientFieldControlPoints->GetSpacing() );
    lambda->SetRegions( this->m_GradientFieldControlPoints
      ->GetLargestPossibleRegion().GetSize() );
    lambda->Allocate();
    lambda->FillBuffer( gradientStep );

    typedef MultiplyImageFilter<ControlPointLatticeType, RealImageType,
      ControlPointLatticeType> MultiplierType;
    typename MultiplierType::Pointer multiplier = MultiplierType::New();
    multiplier->SetInput1( this->m_GradientFieldControlPoints );
    multiplier->SetInput2( lambda );
    multiplier->Update();

    typedef AddImageFilter<ControlPointLatticeType, ControlPointLatticeType,
      ControlPointLatticeType> AdderType;
    typename AdderType::Pointer adder = AdderType::New();
    adder->SetInput1( multiplier->GetOutput() );
    adder->SetInput2( this->m_TotalDeformationFieldControlPoints );
    adder->Update();

    this->m_TotalDeformationFieldControlPoints = adder->GetOutput();

    if( its == this->m_MaximumNumberOfIterations[this->m_CurrentLevel] )
      {
      return;
      }

    this->EvaluateMetricAndGradient();

    RealType gg = 0.0;
    RealType dgg = 0.0;

    ImageRegionIterator<ControlPointLatticeType>
      ItG( G, G->GetLargestPossibleRegion() );
    ImageRegionIterator<ControlPointLatticeType>
      ItH( H, H->GetLargestPossibleRegion() );
    ImageRegionIterator<ControlPointLatticeType>
      ItXi( this->m_GradientFieldControlPoints,
      this->m_GradientFieldControlPoints->GetLargestPossibleRegion() );

    for( ItG.GoToBegin(), ItXi.GoToBegin(); !ItG.IsAtEnd(); ++ItG, ++ItXi )
      {
      VectorType xi = ItXi.Get();
      VectorType g = ItG.Get();
      gg += g * g;
      // Polak-Ribiere
      dgg += ( xi + g ) * xi;
      }
    if( gg == 0.0 )
      {
      continue;
      }

    RealType gamma = vnl_math_max( static_cast<RealType>( 0.0 ),
      static_cast<RealType>( dgg / gg ) );

    ItG.GoToBegin();
    ItH.GoToBegin();
    ItXi.GoToBegin();
    while( !ItG.IsAtEnd() )
      {
      ItG.Set( -ItXi.Get() );
      VectorType tmp = ItG.Get() + gamma*ItH.Get();
      ItH.Set( tmp );
      ItXi.Set( tmp );
      ++ItG;
      ++ItH;
      ++ItXi;
      }
    }
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
typename DMFFDLabeledPointSetRegistrationICPFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>::RealType
DMFFDLabeledPointSetRegistrationICPFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::EvaluateMetric( RealType t = 0 )
{
  /**
   * Warp the moving points
   */
  typename MovingPointSetType::Pointer warpedPoints
    = MovingPointSetType::New();
  warpedPoints->Initialize();

  typename RealImageType::Pointer lambda = RealImageType::New();
  lambda->SetRegions( this->m_TotalDeformationFieldControlPoints
    ->GetLargestPossibleRegion().GetSize() );
  lambda->Allocate();
  lambda->FillBuffer( t );

  typedef MultiplyImageFilter<ControlPointLatticeType, RealImageType,
    ControlPointLatticeType> MultiplierType;
  typename MultiplierType::Pointer multiplier = MultiplierType::New();
  multiplier->SetInput1( this->m_GradientFieldControlPoints );
  multiplier->SetInput2( lambda );
  multiplier->Update();

  typedef AddImageFilter<ControlPointLatticeType, ControlPointLatticeType,
    ControlPointLatticeType> AdderType;
  typename AdderType::Pointer adder = AdderType::New();
  adder->SetInput1( multiplier->GetOutput() );
  adder->SetInput2( this->m_TotalDeformationFieldControlPoints );
  adder->Update();

  typename ControlPointFilterType::Pointer bspliner
    = ControlPointFilterType::New();
  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetInput( adder->GetOutput() );
  bspliner->SetOrigin( this->m_Origin );
  bspliner->SetSize( this->m_Size );
  bspliner->SetSpacing( this->m_Spacing );
  bspliner->SetDirection( this->m_Direction );

  unsigned long count = 0;
  typename MovingPointSetType::PointsContainerConstIterator ItM =
    this->m_SimilarityTransformMovingPointSet->GetPoints()->Begin();
  typename MovingPointSetType::PointDataContainer::ConstIterator ItD =
    this->m_SimilarityTransformMovingPointSet->GetPointData()->Begin();
  while( ItM != this->m_SimilarityTransformMovingPointSet->GetPoints()->End() )
    {
    MovingPointType inputPoint = ItM.Value();

    typename ControlPointFilterType::PointType point;
    point.CastFrom( inputPoint );

    VectorType vector;
    try
      {
      bspliner->EvaluateAtPoint( point, vector );
      }
    catch(...)
      {
      ++ItM;
      ++ItD;
      continue;
      }
    inputPoint += vector;

    warpedPoints->SetPoint( count, inputPoint );
    warpedPoints->SetPointData( count, ItD.Value() );
    count++;
    ++ItM;
    ++ItD;
    }

  /**
   * Set up the point-set function
   */
  typedef itk::EuclideanDistanceLabeledPointSetMetric
    <FixedPointSetType> PointSetFunctionType;
  typename PointSetFunctionType::Pointer pointSetFunction
    = PointSetFunctionType::New();

  pointSetFunction->SetUseWithRespectToTheMovingPointSet( true );
  pointSetFunction->SetFixedPointSet( this->GetInput( 0 ) );
  pointSetFunction->SetMovingPointSet( warpedPoints );

  pointSetFunction->Initialize();

  typename PointSetFunctionType
    ::DefaultTransformType::ParametersType parameters;
  parameters.Fill( 0 );

  typename PointSetFunctionType::MeasureType value
    = pointSetFunction->GetValue( parameters );

  RealType sumValue = 0.0;
  for( unsigned int n = 0; n < value.Size(); n++ )
    {
    sumValue += value[n];
    }

  return sumValue;
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
typename DMFFDLabeledPointSetRegistrationICPFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>::RealType
DMFFDLabeledPointSetRegistrationICPFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::EvaluateMetricAndGradient()
{
  /**
   * Warp the moving points
   */
  typename MovingPointSetType::Pointer warpedPoints
    = MovingPointSetType::New();
  warpedPoints->Initialize();

  typename ControlPointFilterType::Pointer bspliner
    = ControlPointFilterType::New();

  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetInput( this->m_TotalDeformationFieldControlPoints );
  bspliner->SetOrigin( this->m_Origin );
  bspliner->SetSize( this->m_Size );
  bspliner->SetSpacing( this->m_Spacing );
  bspliner->SetDirection( this->m_Direction );

  unsigned long count = 0;
  typename MovingPointSetType::PointsContainerConstIterator ItM =
    this->m_SimilarityTransformMovingPointSet->GetPoints()->Begin();
  typename MovingPointSetType::PointDataContainer::ConstIterator ItD =
    this->m_SimilarityTransformMovingPointSet->GetPointData()->Begin();
  while( ItM != this->m_SimilarityTransformMovingPointSet->GetPoints()->End() )
    {
    MovingPointType inputPoint = ItM.Value();

    typename ControlPointFilterType::PointType point;
    point.CastFrom( inputPoint );

    VectorType vector;
    try
      {
      bspliner->EvaluateAtPoint( point, vector );
      }
    catch(...)
      {
      ++ItM;
      ++ItD;
      continue;
      }
    inputPoint += vector;

    warpedPoints->SetPoint( count, inputPoint );
    warpedPoints->SetPointData( count, ItD.Value() );
    count++;
    ++ItM;
    ++ItD;
    }

  /**
   * Set up the point-set function
   */
  typedef itk::EuclideanDistanceLabeledPointSetMetric
    <FixedPointSetType> PointSetFunctionType;
  typename PointSetFunctionType::Pointer pointSetFunction
    = PointSetFunctionType::New();

  pointSetFunction->SetUseWithRespectToTheMovingPointSet( true );
  pointSetFunction->SetFixedPointSet( this->GetInput( 0 ) );
  pointSetFunction->SetMovingPointSet( warpedPoints );

  pointSetFunction->Initialize();

  typename PointSetFunctionType::DefaultTransformType
    ::ParametersType parameters;
  parameters.Fill( 0 );

  typename PointSetFunctionType::DerivativeType pointSetGradient;
  typename PointSetFunctionType::MeasureType value;
  pointSetFunction->GetValueAndDerivative(
    parameters, value, pointSetGradient );

  RealType sumValue = 0.0;
  for( unsigned int n = 0; n < value.Size(); n++ )
    {
    sumValue += value[n];
    }

  typename BSplinePointSetType::Pointer fieldPoints
    = BSplinePointSetType::New();
  fieldPoints->Initialize();
  typename BSplineFilterType::WeightsContainerType::Pointer weights
    = WeightsContainerType::New();
  weights->Initialize();

  count = 0;
  for( unsigned int n = 0; n < pointSetFunction->GetNumberOfValues(); n++ )
    {
    typename BSplinePointSetType::PointType fieldPoint;
    VectorType gradient;

    MovingPointType warpedPoint;
    warpedPoints->GetPoint( n, &warpedPoint );

    WarpedPixelType label = 1;
    warpedPoints->GetPointData( n, &label );

    bool isInside = true;
    for( unsigned d = 0; d < Dimension; d++ )
      {
      gradient[d] = pointSetGradient(n, d);
      if( !this->m_Directionality[d] )
        {
        gradient[d] = 0.0;
        }
      fieldPoint[d] = warpedPoint[d];
      if( fieldPoint[d] <= this->m_Origin[d] ||
          fieldPoint[d] >= this->m_Origin[d] + this->m_Spacing[d]
        * static_cast<RealType>( this->m_Size[d] - 1 ) )
        {
        isInside = false;
        }
      }
    if( isInside )
      {
      fieldPoints->SetPoint( count, fieldPoint );
      fieldPoints->SetPointData( count, gradient );
      if( this->m_LabelWeights.size() == 0 )
        {
        weights->InsertElement( count, 1.0 );
        }
      else
        {
        weights->InsertElement( count, this->m_LabelWeights[label] );
        }
      count++;
      }
    }

  /**
    * Since the different metrics have different gradient scales,
    * we must scale the gradient so that an effective line search
    * is carried out.  We take the limiting case of gradients with
    * normal i.i.d. components which is a Rayleigh distribution for
    * 2-D images and a Maxwell distribution for 3-D images.  The
    * sigma parameter which governs both distributions is chosen
    * such that the mode of the distribution is equal to the spacing
    * of one pixel.
    */
  itkDebugMacro( "Normalizing gradient values." );

  RealType sumSquaredNorm = 0.0;
  for( unsigned int i = 0; i < fieldPoints->GetNumberOfPoints(); i++ )
    {
    VectorType gradient;
    fieldPoints->GetPointData( i, &gradient );
    sumSquaredNorm += gradient.GetSquaredNorm();
    }
  VectorType V;
  RealType sigma;
  for( unsigned int i = 0; i < Dimension; i++ )
    {
    V[i] = this->GetSpacing()[i];
    }
  if( Dimension == 2 )
    {
    sigma = V.GetNorm();
    }
  else if( Dimension == 3 )
    {
    sigma = V.GetNorm()/vcl_sqrt( 2.0 );
    }
  RealType gradientNormalizationFactor = sigma*vcl_sqrt
    ( static_cast<RealType>( Dimension *
    fieldPoints->GetNumberOfPoints() ) / sumSquaredNorm );

  for( unsigned int i = 0; i < fieldPoints->GetNumberOfPoints(); i++ )
    {
    VectorType gradient;
    fieldPoints->GetPointData( i, &gradient );
    fieldPoints->SetPointData( i,
      gradient*gradientNormalizationFactor
      *this->m_GradientScalingFactor[this->m_CurrentLevel] );
    }

  typename BSplineFilterType::ArrayType nlevels;
  typename BSplineFilterType::ArrayType ncps;

  nlevels.Fill( 1 );
  for( unsigned int d = 0; d < Dimension; d++ )
    {
    ncps[d] = this->m_TotalDeformationFieldControlPoints
      ->GetLargestPossibleRegion().GetSize()[d];
    }

  typename BSplineFilterType::Pointer bsplineFilter = BSplineFilterType::New();

  bsplineFilter->SetInput( fieldPoints );
  bsplineFilter->SetOrigin( this->m_Origin );
  bsplineFilter->SetSpacing( this->m_Spacing );
  bsplineFilter->SetSize( this->m_Size );
  bsplineFilter->SetDirection( this->m_Direction );
  bsplineFilter->SetNumberOfLevels( nlevels );
  bsplineFilter->SetSplineOrder( this->m_SplineOrder );
  bsplineFilter->SetNumberOfControlPoints( ncps );
  bsplineFilter->SetGenerateOutputImage( false );
  bsplineFilter->SetPointWeights( weights.GetPointer() );
  bsplineFilter->Update();

  this->m_GradientFieldControlPoints = bsplineFilter->GetPhiLattice();

  return sumValue;
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
void
DMFFDLabeledPointSetRegistrationICPFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Transformation variables" << std::endl;
  os << indent << "   " << "Spline order: "
    << this->m_SplineOrder << std::endl;
  os << indent << "  " << "Initial mesh resolution: "
    << this->m_InitialMeshResolution << std::endl;
  os << indent << "  " << "Origin: "
    << this->m_Origin << std::endl;
  os << indent << "  " << "Spacing: "
    << this->m_Spacing << std::endl;
  os << indent << "  " << "Size: "
    << this->m_Size << std::endl;

  os << indent << "Optimization variables" << std::endl;
  os << indent << "  " << "Maximum number of iterations: "
    << this->m_MaximumNumberOfIterations << std::endl;
  os << indent << "  " << "Gradient scaling factor: "
    << this->m_GradientScalingFactor << std::endl;
  os << indent << "  " << "Line search maximum iterations: "
    << this->m_LineSearchMaximumIterations << std::endl;
  os << indent << "  " << "Line search maximum step size: "
    << this->m_LineSearchMaximumStepSize << std::endl;
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
void
DMFFDLabeledPointSetRegistrationICPFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::LineMinimization( RealType *step, RealType *fret )
{
  const RealType Gold = 1.618034;

  if( this->m_Verbose )
    {
    std::cout << "    Begin line search..." << std::endl;
    }
  RealType ax, bx, fbx, cx;
  fbx = this->FindBracketingTriplet( &ax, &bx, &cx );

  if( ax > this->m_LineSearchMaximumStepSize )
    {
    *step = this->m_LineSearchMaximumStepSize;
    *fret = this->EvaluateEnergyForLineSearch( *step );
    if( this->m_Verbose )
      {
      std::cout << "      ***** exceeded line search maximum step size. ***** " << std::endl;
      }
    itkDebugMacro( "      ***** exceeded line search maximum step size. ***** " );

    ax = 0;
    cx = this->m_LineSearchMaximumStepSize;
    bx = ( cx + Gold*ax )/( 1.0 + Gold );
    fbx = this->EvaluateEnergyForLineSearch( bx );
    }
  this->BrentSearch( ax, bx, fbx, cx, step, fret );
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
typename DMFFDLabeledPointSetRegistrationICPFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>::RealType
DMFFDLabeledPointSetRegistrationICPFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::EvaluateEnergyForLineSearch( RealType lambda )
{
  return this->EvaluateMetric( lambda );
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
typename DMFFDLabeledPointSetRegistrationICPFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>::RealType
DMFFDLabeledPointSetRegistrationICPFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::FindBracketingTriplet( RealType *ax, RealType *bx, RealType *cx )
{
  const RealType Gold = 1.618034;
  const RealType Glimit = 100;
  const RealType Tiny = 1e-20;
  *ax = 0.0;
  *bx = 1.0;

  RealType fa = this->EvaluateEnergyForLineSearch( *ax );
  RealType fb = this->EvaluateEnergyForLineSearch( *bx );

  RealType dum;
  if( fb > fa )
    {
    dum = *ax;
    *ax = *bx;
    *bx = dum;
    dum = fb;
    fb = fa;
    fa = dum;
    }
  *cx = *bx + Gold*( *bx-*ax );
  RealType fc = this->EvaluateEnergyForLineSearch( *cx );
  if( this->m_Verbose )
    {
    if( *cx < *ax )
      {
      std::cout << "      Bracket triple: "
                << "f(" << *cx << ") = " << fc << ", "
                << "f(" << *bx << ") = " << fb << ", "
                << "f(" << *ax << ") = " << fa << std::endl;
      if( *cx > this->m_LineSearchMaximumStepSize )
        {
        return fb;
        }
      }
    else
      {
      std::cout << "      Bracket triple: "
                << "f(" << *ax << ") = " << fa << ", "
                << "f(" << *bx << ") = " << fb << ", "
                << "f(" << *cx << ") = " << fc << std::endl;
      if( *ax > this->m_LineSearchMaximumStepSize )
        {
        RealType u = *ax;
        *ax = *cx;
        *cx = u;
        RealType fu = fa;
        fa = fc;
        fc = fu;
        return fb;
        }
      }
    }
  RealType ulim, u, r, q, fu;

  while( fb > fc )
    {
    r = ( *bx-*ax )*( fb-fc );
    q = ( *bx-*cx )*( fb-fa );
    RealType denom = 2.0*vnl_math_max( vnl_math_abs( q-r ), Tiny );
    if( q-r < 0.0 )
      {
      denom = -denom;
      }
    u = *bx - ( ( *bx-*cx )*q - ( *bx-*ax )*r )/denom;
    ulim = *bx + Glimit*( *cx-*bx );
    if( ( *bx-u )*( u-*cx ) > 0.0 )
      {
      fu = this->EvaluateEnergyForLineSearch( u );
      if( fu < fc )
        {
        *ax = *bx;
        *bx = u;
        fa = fb;
        fb = fu;
        if( this->m_Verbose )
          {
          if( *cx < *ax )
            {
            std::cout << "      Bracket triple: "
                      << "f(" << *cx << ") = " << fc << ", "
                      << "f(" << *bx << ") = " << fb << ", "
                      << "f(" << *ax << ") = " << fa << std::endl;
            if( *cx > this->m_LineSearchMaximumStepSize )
              {
              return fb;
              }
            }
          else
            {
            std::cout << "      Bracket triple: "
                      << "f(" << *ax << ") = " << fa << ", "
                      << "f(" << *bx << ") = " << fb << ", "
                      << "f(" << *cx << ") = " << fc << std::endl;
            if( *ax > this->m_LineSearchMaximumStepSize )
              {
              u = *ax;
              *ax = *cx;
              *cx = u;
              fu = fa;
              fa = fc;
              fc = fu;
              return fb;
              }
            }
          }
        return fb;
        }
      else if( fu > fb )
        {
        *cx = u;
        fc = fu;
        if( this->m_Verbose )
          {
          if( *cx < *ax )
            {
            std::cout << "      Bracket triple: "
                      << "f(" << *cx << ") = " << fc << ", "
                      << "f(" << *bx << ") = " << fb << ", "
                      << "f(" << *ax << ") = " << fa << std::endl;
            if( *cx > this->m_LineSearchMaximumStepSize )
              {
              return fb;
              }
            }
          else
            {
            std::cout << "      Bracket triple: "
                      << "f(" << *ax << ") = " << fa << ", "
                      << "f(" << *bx << ") = " << fb << ", "
                      << "f(" << *cx << ") = " << fc << std::endl;
            if( *ax > this->m_LineSearchMaximumStepSize )
              {
              u = *ax;
              *ax = *cx;
              *cx = u;
              fu = fa;
              fa = fc;
              fc = fu;
              return fb;
              }
            }
          }
        return fb;
        }
      u = *cx + Gold*( *cx-*bx );
      fu = this->EvaluateEnergyForLineSearch( u );
      }
    else if( ( *cx-u )*( u-ulim ) > 0.0 )
      {
      fu = this->EvaluateEnergyForLineSearch( u );
      if( fu < fc )
        {
        *bx = *cx;
        *cx = u;
        u = *cx + Gold*( *cx-*bx );
        fb = fc;
        fc = fu;
        fu = this->EvaluateEnergyForLineSearch( u );
        }
      }
    else if( ( u-ulim )*( ulim-*cx ) >= 0.0 )
      {
      u = ulim;
      fu = this->EvaluateEnergyForLineSearch( u );
      }
    else
      {
      u = *cx + Gold*( *cx-*bx );
      fu = this->EvaluateEnergyForLineSearch( u );
      }

    *ax = *bx;
    *bx = *cx;
    *cx = u;
    fa = fb;
    fb = fc;
    fc = fu;
    if( this->m_Verbose )
      {
      if( *cx < *ax )
        {
        std::cout << "      Bracket triple: "
                  << "f(" << *cx << ") = " << fc << ", "
                  << "f(" << *bx << ") = " << fb << ", "
                  << "f(" << *ax << ") = " << fa << std::endl;
        if( *cx > this->m_LineSearchMaximumStepSize )
          {
          return fb;
          }
        }
      else
        {
        std::cout << "      Bracket triple: "
                  << "f(" << *ax << ") = " << fa << ", "
                  << "f(" << *bx << ") = " << fb << ", "
                  << "f(" << *cx << ") = " << fc << std::endl;
        if( *ax > this->m_LineSearchMaximumStepSize )
          {
          u = *ax;
          *ax = *cx;
          *cx = u;
          fu = fa;
          fa = fc;
          fc = fu;
          return fb;
          }
        }
      }
    }
  return fb;
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
void
DMFFDLabeledPointSetRegistrationICPFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::BrentSearch( RealType ax, RealType bx, RealType fbx,
  RealType cx, RealType *step, RealType *fret )
{
  const RealType R = 0.6180339;
  const RealType CGOLD = 1.0 - R;
  const RealType tol = 1e-20;
  const RealType ZEPS = 1.0e-10;

  RealType a;
  RealType b;
  RealType d = 0;
  RealType p;
  RealType q;
  RealType r;
  RealType etemp;
  RealType e = 0.0;

  if( ax < cx )
    {
    a = ax;
    b = cx;
    }
  else
    {
    a = cx;
    b = ax;
    }
  RealType x, fx;
  RealType w, fw;
  RealType v, fv;
  RealType u, fu;
  x = w = v = bx;
  fx = fw = fv = fbx; //this->EvaluateEnergyForLineSearch( x );

  for( unsigned int iter = 1;
    iter <= this->m_LineSearchMaximumIterations; iter++ )
    {
    if( this->m_Verbose )
      {
      std::cout << "        Line iteration " << iter
                << ": f(x = " << x << ") = " << fx << ", x in ["
                << a << ", " << b << "] " << std::endl;
      }
    if( a > this->m_LineSearchMaximumStepSize )
      {
      *step = this->m_LineSearchMaximumStepSize;
      *fret = this->EvaluateEnergyForLineSearch( *step );
      if( this->m_Verbose )
        {
        std::cout <<  "    Results of line search: E(" << *step << ") = "
        << *fret << " (exceeded line search maximum step size)." << std::endl;
        }
      itkDebugMacro( "    Results of line search: E(" << *step << ") = "
        << *fret << " (exceeded line search maximum step size)." );
      return;
      }

    RealType xm = 0.5 * ( a + b );
    RealType tol1 = tol * vnl_math_abs( x ) + ZEPS;
    RealType tol2 = 2.0 * tol1;
    if( vnl_math_abs( x - xm ) <= ( tol2 - 0.5 * ( b - a ) ) )
      {
      if( x > this->m_LineSearchMaximumStepSize )
        {
        *step = this->m_LineSearchMaximumStepSize;
        *fret = this->EvaluateEnergyForLineSearch( *step );
        if( this->m_Verbose )
          {
          std::cout <<  "    Results of line search: E(" << *step << ") = "
          << *fret << " (exceeded line search maximum step size)." << std::endl;
          }
        itkDebugMacro( "    Results of line search: E(" << *step << ") = "
          << *fret << " (exceeded line search maximum step size)." );
        return;
        }
      else
        {
        if( this->m_Verbose )
          {
          std::cout <<  "    Results of line search: E(" << x << ") = "
            << fx << "." << std::endl;
          }
        itkDebugMacro( << "Results of line search: E(" << x << ") = "
          << fx << "." );
        *step = x;
        *fret = fx;
        return;
        }
      }
    if( vnl_math_abs( e ) > tol1 )
      {
      r = ( x - w ) * ( fx - fv );
      q = ( x - v ) * ( fx - fw );
      p = ( x - v ) * q - ( x - w ) * r;
      q = 2.0 * ( q - r );
      if( q > 0.0 )
        {
        p = -p;
        }
      q = vnl_math_abs( q );
      etemp = e;
      e = d;
      if( vnl_math_abs( p ) >= vnl_math_abs( 0.5 * q * etemp )
           || p <= q * ( a - x ) || p >= q * ( b - x ) )
        {
        if( x >= xm )
          {
          e = a - x;
          }
        else
          {
          e = b - x;
          }
        d = CGOLD * e;
        }
      else
        {
        d = p / q;
        u = x + d;
        if( u - a < tol2 || b - u < tol2 )
          {
          d = vnl_math_abs( tol1 );
          if( xm - x <= 0 )
            {
            d = -d;
            }
          }
        }
      }
    else
      {
      if( x >= xm )
        {
        e = a - x;
        }
      else
        {
        e = b - x;
        }
      d = CGOLD * e;
      }
    if( vnl_math_abs( d ) >= tol1 )
      {
      u = x + d;
      }
    else
      {
      u = x + vnl_math_abs( tol1 );
      if( d <= 0 )
        {
        u = x - vnl_math_abs( tol1 );
        }
      }
    fu = this->EvaluateEnergyForLineSearch( u );
    if( fu <= fx )
      {
      if( u >= x )
        {
        a = x;
        }
      else
        {
        b = x;
        }
      v = w;
      w = x;
      x = u;
      fv = fw;
      fw = fx;
      fx = fu;
      }
    else
      {
      if( u < x )
        {
        a = u;
        }
      else
        {
        b = u;
        }
      if( fu <= fw || w == x )
        {
        v = w;
        w = u;
        fv = fw;
        fw = fu;
        }
      else if( fu <= fv || v == x || v == w )
        {
        v = u;
        fv = fu;
        }
      }
    }

  if( x > this->m_LineSearchMaximumStepSize )
    {
    *step = this->m_LineSearchMaximumStepSize;
    *fret = this->EvaluateEnergyForLineSearch( *step );
    if( this->m_Verbose )
      {
      std::cout <<  "    Results of line search: E(" << *step << ") = "
      << *fret << " (exceeded line search maximum step size)." << std::endl;
      }
    itkDebugMacro( "    Results of line search: E(" << *step << ") = "
      << *fret << " (exceeded line search maximum step size)." );
    return;
    }
  else
    {
    if( this->m_Verbose )
      {
      std::cout <<  "    Results of line search: E(" << x << ") = "
        << fx << "." << std::endl;
      }
    itkDebugMacro( << "Results of line search: E(" << x << ") = "
      << fx << "." );
    *step = x;
    *fret = fx;
    return;
    }
}

}
#endif

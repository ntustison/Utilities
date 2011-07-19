/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: itkFFDPointSetRegistrationFilter.hxx,v $
Language:  C++

Date:      $Date: 2009/01/14 02:56:30 $
Version:   $Revision: 1.14 $

Copyright (c) Insight Software Consortium. All rights reser
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for detail.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef _itkFFDPointSetRegistrationFilter_hxx_
#define _itkFFDPointSetRegistrationFilter_hxx_

// disable debug warnings in MS compiler
#ifdef _MSC_VER
#pragma warning(disable: 4786)
#endif

#include "itkFFDPointSetRegistrationFilter.h"

#include "itkAddImageFilter.h"
#include "itkJensenHavrdaCharvatTsallisLabeledPointSetMetric.h"
#include "itkMultiplyImageFilter.h"
#include "itkTimeProbe.h"

#include "itkImageDuplicator.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionExclusionConstIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVectorImageFileWriter.h"

#include "itkTimeProbe.h"

#include <fstream.h>

#include "vnl/vnl_math.h"

namespace itk {

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
FFDPointSetRegistrationFilter<TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::FFDPointSetRegistrationFilter()
{

//  if( FixedPointSetType::Dimension !=
//         MovingPointSetType::Dimension )
//    {
//    itkExceptionMacro( "Point set dimensions must be equal." );
//    }

  this->SetNumberOfRequiredInputs( 2 );

  // Default values
  this->SetNumberOfLevels( 3 );
  this->SetLineSearchMaximumIterations( 10 );
  this->m_MaximumNumberOfIterations.Fill( 10 );
  this->SetSplineOrder( 3 );
  this->SetNumberOfFixedSamples( 1000 );
  this->SetNumberOfMovingSamples( 1000 );

  this->SetEmploySteepestDescent( false );
  this->SetUseInputAsSamples( true );
  this->SetUseAnisotropicCovariances( false );
  this->SetProlificacy( true );
  this->SetWhichGradient( 0 );
  this->m_Directionality.Fill( 1 );
  this->SetEmployTerm2( true );

  this->m_RegularizationSigma = 1.0;
  this->m_CovarianceKNeighborhood = 4;
  this->m_EvaluationKNeighborhood = 50;
  this->m_ExpansionFactor = 0.5;

  this->m_Size.Fill( 0 );

  this->m_AnnealingRate = 0.93;
  this->m_FilePrefix = "output";
  
  this->SetBindBoundary( true );

  this->SetAlpha( 1.0 );

  this->m_InitialDeformationFieldControlPoints = ControlPointLatticeType::New();
  this->m_InitialDeformationFieldControlPoints = NULL;

  this->m_Randomizer = RandomizerType::New();
  this->m_Randomizer->SetSeed(static_cast <ITK_UINT32> (0));
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
FFDPointSetRegistrationFilter<TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::~FFDPointSetRegistrationFilter()
{
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
void
FFDPointSetRegistrationFilter<TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::SetNumberOfLevels( unsigned int n )
{

  if( n != this->m_NumberOfLevels )
    {
    this->m_NumberOfLevels = n;

    unsigned int tmp_MI = 1;

    if( this->m_MaximumNumberOfIterations.Size() > 0 )
      {
      tmp_MI = this->m_MaximumNumberOfIterations[0];
      }

    this->m_MaximumNumberOfIterations.SetSize( n );

    if( this->m_NumberOfLevels != 0 )
      {
      this->m_MaximumNumberOfIterations.Fill( tmp_MI );
      }
    this->Modified();
  }
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
void
FFDPointSetRegistrationFilter<TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::GenerateData()
{
  if( this->GetNumberOfInputs() < 2 )
    {
    itkExceptionMacro( << "Images are not specified." );
    }

  this->m_CurrentAnnealing = 1.0;

  for( this->m_CurrentLevel = 0;
        this->m_CurrentLevel < this->m_NumberOfLevels; this->m_CurrentLevel++ )
    {
    this->Initialize();

    std::cout << "Current level = " << this->m_CurrentLevel+1 << " ("
      << this->m_NumberOfLevels << " total levels).  " << "Grid size = "
      << this->m_TotalDeformationFieldControlPoints
        ->GetLargestPossibleRegion().GetSize() << ". " << std::endl;

    itkDebugMacro( "Current level = " << this->m_CurrentLevel+1 
      << " (" << this->m_NumberOfLevels << " total levels).  "
      << "Grid size = " << this->m_TotalDeformationFieldControlPoints
        ->GetLargestPossibleRegion().GetSize() << ". " );

    if( this->m_MaximumNumberOfIterations[this->m_CurrentLevel] > 0 )
      {
      itkDebugMacro( "Begin iterative solve. " );
      this->IterativeSolve();
      itkDebugMacro( "Finish iterative solve. " );
      }

    typedef VectorImageFileWriter
      <ControlPointLatticeType, RealImageType> WriterType;

    itk::OStringStream buf;
    buf << "_" << "_" << this->m_CurrentLevel << "_" << this->GetSplineOrder();

    std::string filename = std::string( this->m_FilePrefix )
      + std::string( "ControlPointLattice" ) + buf.str()
      + std::string( ".nii.gz" );

    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( filename.c_str() );
    writer->SetInput( this->m_TotalDeformationFieldControlPoints );
    writer->Update();
    }
  // Warp the input points

  typename WarpedPointSetType::Pointer warpedPoints = WarpedPointSetType::New();
  warpedPoints->Initialize();

  typename BSplineControlPointFilterType::Pointer bspliner
    = BSplineControlPointFilterType::New();
  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetInput( this->m_TotalDeformationFieldControlPoints );
  bspliner->SetOrigin( this->m_Origin );
  bspliner->SetSize( this->m_Size );
  bspliner->SetSpacing( this->m_Spacing );

  unsigned long count = 0;
  for( unsigned int j = 0; j < this->GetInput( 1 )->GetNumberOfPoints(); j++ )
    {

    MovingPointType inputPoint;
    this->GetInput( 1 )->GetPoint( j, &inputPoint );

    typename BSplineControlPointFilterType::PointType point;
    for( unsigned int d = 0; d < Dimension; d++ )
      {
      point[d] = inputPoint[d];
      }
    VectorType vector;
    try 
      {
      bspliner->EvaluateAtPoint( point, vector );
      }
    catch(...)
      {
      continue; 
      }  
    point += vector;
    for( unsigned int d = 0; d < Dimension; d++ )
      {
      inputPoint[d] = point[d];
      }
    typename MovingPointSetType::PixelType inputLabel = 1;
    this->GetInput( 1 )->GetPointData( j, &inputLabel );
    
    warpedPoints->SetPoint( j, inputPoint );
    warpedPoints->SetPointData( count, inputLabel );
    count++;
    }
  this->SetNthOutput( 0, warpedPoints );
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
void
FFDPointSetRegistrationFilter<TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::Initialize()
{

  if( this->m_CurrentLevel == 0 )
    {

    if( this->m_InitialDeformationFieldControlPoints )
      {
      for( unsigned int d = 0; d < Dimension; d++ )
        {
        this->m_InitialMeshResolution[d] =
          this->m_InitialDeformationFieldControlPoints->
            GetLargestPossibleRegion().GetSize()[d] - this->m_SplineOrder;

        if( this->m_InitialMeshResolution[d] < 1 )
          {
          std::cerr << "Invalid size for initial deformation field control "
            << "point lattice." << std::endl;
          exit( 0 );
          }
        }
      }

    typename ControlPointLatticeType::RegionType::SizeType size;
    for( unsigned int d = 0; d < Dimension; d++ )
      {
      size[d] = this->m_InitialMeshResolution[d] + this->m_SplineOrder;
      }

    VectorType V;
    V.Fill( 0 );
    if( this->m_InitialDeformationFieldControlPoints )
      {
      typedef ImageDuplicator<ControlPointLatticeType> DuplicatorType;
      typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
      duplicator->SetInputImage( this->m_InitialDeformationFieldControlPoints );
      duplicator->Update();
      this->m_TotalDeformationFieldControlPoints = duplicator->GetOutput();
      }
    else
      {
      this->m_TotalDeformationFieldControlPoints = ControlPointLatticeType::New();
      this->m_TotalDeformationFieldControlPoints->SetRegions( size );
      this->m_TotalDeformationFieldControlPoints->Allocate();
      this->m_TotalDeformationFieldControlPoints->FillBuffer( V );
      }

    this->m_GradientFieldControlPoints = ControlPointLatticeType::New();
    this->m_GradientFieldControlPoints->SetRegions( size );
    this->m_GradientFieldControlPoints->Allocate();
    this->m_GradientFieldControlPoints->FillBuffer( V );


    if( this->m_Size[0] == 0 )
      {

      /**
       * define the transformation domain
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

      /**
       * The size and spacing are irrelevant other 
       * than they define the parametric domain
       */
      this->m_Size.Fill( 100 );

      for( unsigned int d = 0; d < Dimension; d++ )
        {
        /**
         * Expand the bounding box to ensure coverage
         */
        RealType expandedMinimum = minPoint[d] 
          - this->m_ExpansionFactor * ( maxPoint[d] - minPoint[d] );
        RealType expandedMaximum = maxPoint[d] 
          + this->m_ExpansionFactor * ( maxPoint[d] - minPoint[d] );
        this->m_Origin[d] = expandedMinimum;
        this->m_Spacing[d] = ( expandedMaximum - expandedMinimum ) 
          / static_cast<RealType>( this->m_Size[d] - 1 );
        }
      }

    for( unsigned int d = 0; d < Dimension; d++ )
      {
      this->m_MinBoundary[d] = this->m_Origin[d] + this->m_Spacing[d];
      this->m_MaxBoundary[d] = this->m_Origin[d] + this->m_Spacing[d]
        * static_cast<RealType>( this->m_Size[d] - 1 ) - this->m_Spacing[d];
      }

    std::cout << "Transformation domain: " << std::endl;
    std::cout << "   Origin:      " << this->m_Origin << std::endl;
    std::cout << "   Size:        " << this->m_Size << std::endl;
    std::cout << "   Spacing:     " << this->m_Spacing << std::endl;
    std::cout << "   MinBoundary: " << this->m_MinBoundary << std::endl;
    std::cout << "   MaxBoundary: " << this->m_MaxBoundary << std::endl;
    }
  else
    {
    typedef BSplineControlPointImageFilter
      <ControlPointLatticeType, DeformationFieldType> 
      BSplineControlPointFilterType;
    typename BSplineControlPointFilterType::Pointer bspliner
      = BSplineControlPointFilterType::New();

    bspliner->SetSplineOrder( this->m_SplineOrder );
    bspliner->SetInput( this->m_TotalDeformationFieldControlPoints );

    typename BSplineFilterType::ArrayType nlevels;
    nlevels.Fill( 2 );
    this->m_TotalDeformationFieldControlPoints 
      = bspliner->RefineControlLattice( nlevels );

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
FFDPointSetRegistrationFilter<TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::IterativeSolve()
{
  const RealType ftol = 1.0e-50;
  const RealType eps = 1.0e-50;

  this->m_CurrentIteration = -1;

  RealType fp = this->EvaluateMetricAndGradient();

  typename ControlPointLatticeType::Pointer G = ControlPointLatticeType::New();
  typename ControlPointLatticeType::Pointer H = ControlPointLatticeType::New();

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
    this->m_CurrentAnnealing *= this->m_AnnealingRate;

    std::cout << "  Iteration " << its
      << " (of " << this->m_MaximumNumberOfIterations[this->m_CurrentLevel]
      << "): Energy = " << fp << ", " << " Annealing temperature = "
      << vcl_sqrt( this->m_CurrentAnnealing / this->m_AnnealingRate )
         * this->m_RegularizationSigma << std::endl;
    itkDebugMacro( "  Iteration " << its
      << " (of " << this->m_MaximumNumberOfIterations[this->m_CurrentLevel]
      << "): Energy = " << fp << ", " << " Annealing temperature = "
      << vcl_sqrt( this->m_CurrentAnnealing / this->m_AnnealingRate )
         * this->m_RegularizationSigma );

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

    RealType gamma = 0.0;
    if( !this->m_EmploySteepestDescent )
      {
      gamma = vnl_math_max( static_cast<RealType>( 0.0 ),
        static_cast<RealType>( dgg / gg ) );
      }

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
typename FFDPointSetRegistrationFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>::RealType
FFDPointSetRegistrationFilter<TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::EvaluateMetric( RealType t = 0 )
{

  /**
   * Warp the moving points
   */
  typename MovingPointSetType::Pointer warpedPoints = MovingPointSetType::New();
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

  typename BSplineControlPointFilterType::Pointer bspliner
    = BSplineControlPointFilterType::New();

  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetInput( adder->GetOutput() );
  bspliner->SetOrigin( this->m_Origin );
  bspliner->SetSize( this->m_Size );
  bspliner->SetSpacing( this->m_Spacing );

  unsigned long count = 0;
  for( unsigned int n = 0; n < this->GetInput( 1 )->GetNumberOfPoints(); n++ )
    {
    MovingPointType inputPoint;
    this->GetInput( 1 )->GetPoint( n, &inputPoint );

    typename BSplineControlPointFilterType::PointType point;
    point.CastFrom( inputPoint );

    VectorType vector;
    try
      {
      bspliner->EvaluateAtPoint( point, vector );
      }
    catch(...)
      {
      continue;
      }
    inputPoint += vector;  

//    for( unsigned int d = 0; d < Dimension; d++ )
//      {
//      if( inputPoint[d] < this->m_MinBoundary[d] ||
//             inputPoint[d] > this->m_MaxBoundary[d] )
//        {
//        return NumericTraits<RealType>::max();
//        }
//      }
    typename MovingPointSetType::PixelType inputLabel = 1;
    this->GetInput( 1 )->GetPointData( n, &inputLabel );

    warpedPoints->SetPoint( count, inputPoint );
    warpedPoints->SetPointData( count, inputLabel );
    count++;
    }

  /**
   * Set up the point-set function
   */
  typedef itk::JensenHavrdaCharvatTsallisLabeledPointSetMetric
    <FixedPointSetType> PointSetFunctionType;
  typename PointSetFunctionType::Pointer pointSetFunction
    = PointSetFunctionType::New();

  pointSetFunction->SetUseWithRespectToTheMovingPointSet( true );
  pointSetFunction->SetUseAnisotropicCovariances( 
    this->m_UseAnisotropicCovariances );
  pointSetFunction->SetUseInputAsSamples( 
    this->m_UseInputAsSamples );
  pointSetFunction->SetUseRegularizationTerm( this->m_EmployTerm2 );
  pointSetFunction->SetAlpha( this->m_Alpha );

  pointSetFunction->SetFixedPointSet( this->GetInput( 0 ) );
  pointSetFunction->SetFixedPointSetSigma( this->m_RegularizationSigma
    * vcl_sqrt( this->m_CurrentAnnealing ) );
  pointSetFunction->SetFixedKernelSigma( this->m_KernelSigma );
  pointSetFunction->SetFixedCovarianceKNeighborhood( 
    this->m_CovarianceKNeighborhood );
  pointSetFunction->SetFixedEvaluationKNeighborhood( 
    this->m_EvaluationKNeighborhood );
  pointSetFunction->SetNumberOfFixedSamples( this->m_NumberOfFixedSamples );

  pointSetFunction->SetMovingPointSet( warpedPoints );
  pointSetFunction->SetMovingPointSetSigma( this->m_RegularizationSigma
    * vcl_sqrt( this->m_CurrentAnnealing ) );
  pointSetFunction->SetMovingKernelSigma( this->m_KernelSigma );
  pointSetFunction->SetMovingCovarianceKNeighborhood( 
    this->m_CovarianceKNeighborhood );
  pointSetFunction->SetMovingEvaluationKNeighborhood( 
    this->m_EvaluationKNeighborhood );
  pointSetFunction->SetNumberOfFixedSamples( this->m_NumberOfMovingSamples );

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
typename FFDPointSetRegistrationFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>::RealType
FFDPointSetRegistrationFilter<TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::EvaluateMetricAndGradient()
{
  /**
   * Warp the moving points
   */
  typename MovingPointSetType::Pointer warpedPoints
    = MovingPointSetType::New();
  warpedPoints->Initialize();

  typename BSplineControlPointFilterType::Pointer bspliner
    = BSplineControlPointFilterType::New();

  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetInput( this->m_TotalDeformationFieldControlPoints );
  bspliner->SetOrigin( this->m_Origin );
  bspliner->SetSize( this->m_Size );
  bspliner->SetSpacing( this->m_Spacing );

  unsigned long count = 0;
  for( unsigned int n = 0; n < this->GetInput( 1 )->GetNumberOfPoints(); n++ )
    {
    MovingPointType inputPoint;
    this->GetInput( 1 )->GetPoint( n, &inputPoint );

    typename BSplineControlPointFilterType::PointType point;
    point.CastFrom( inputPoint );

    VectorType vector;
    try
      {
      bspliner->EvaluateAtPoint( point, vector );
      }
    catch(...)
      {
      continue;
      }
    inputPoint += vector;

    typename MovingPointSetType::PixelType inputLabel = 1;
    this->GetInput( 1 )->GetPointData( n, &inputLabel );
    warpedPoints->SetPoint( count, inputPoint );
    warpedPoints->SetPointData( count, inputLabel );
    count++;
    }

  /**
   * Set up the point-set function
   */
  typedef itk::JensenHavrdaCharvatTsallisLabeledPointSetMetric
    <FixedPointSetType> PointSetFunctionType;
  typename PointSetFunctionType::Pointer pointSetFunction
    = PointSetFunctionType::New();

  pointSetFunction->SetUseWithRespectToTheMovingPointSet( true );
  pointSetFunction->SetUseAnisotropicCovariances( 
    this->m_UseAnisotropicCovariances );
  pointSetFunction->SetUseInputAsSamples( this->m_UseInputAsSamples );
  pointSetFunction->SetUseRegularizationTerm( this->m_EmployTerm2 );
  pointSetFunction->SetAlpha( this->m_Alpha );

  pointSetFunction->SetFixedPointSet( this->GetInput( 0 ) );
  pointSetFunction->SetFixedPointSetSigma( this->m_RegularizationSigma
      * vcl_sqrt( this->m_CurrentAnnealing ) );
  pointSetFunction->SetFixedKernelSigma( this->m_KernelSigma );
  pointSetFunction->SetFixedCovarianceKNeighborhood( 
    this->m_CovarianceKNeighborhood );
  pointSetFunction->SetFixedEvaluationKNeighborhood( 
    this->m_EvaluationKNeighborhood );
  pointSetFunction->SetNumberOfFixedSamples( this->m_NumberOfFixedSamples );

  pointSetFunction->SetMovingPointSet( warpedPoints );
  pointSetFunction->SetMovingPointSetSigma( this->m_RegularizationSigma
      * vcl_sqrt( this->m_CurrentAnnealing ) );
  pointSetFunction->SetMovingKernelSigma( this->m_KernelSigma );
  pointSetFunction->SetMovingCovarianceKNeighborhood( 
    this->m_CovarianceKNeighborhood );
  pointSetFunction->SetMovingEvaluationKNeighborhood( 
    this->m_EvaluationKNeighborhood );
  pointSetFunction->SetNumberOfMovingSamples( this->m_NumberOfMovingSamples );

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
      fieldPoint[d] = warpedPoint[d];
      if( fieldPoint[d] <= this->m_MinBoundary[d] ||
          fieldPoint[d] >= this->m_MaxBoundary[d] )
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
    * Since the different metrics have different gradient scales, we must scale
    * the gradient so that an effective line search is carried out.  We take the
    * limiting case of gradients with normal i.i.d. components which is a
    * Rayleigh distribution for 2-D images and a Maxwell distribution for 3-D images.
    * The sigma parameter which governs both distributions is chosen such that the
    * mode of the distribution is equal to the spacing of one pixel.
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
  RealType gradientScalingFactor = sigma*vcl_sqrt
    ( static_cast<RealType>( Dimension *
    fieldPoints->GetNumberOfPoints() ) / sumSquaredNorm );

  for( unsigned int i = 0; i < fieldPoints->GetNumberOfPoints(); i++ )
    {
    VectorType gradient;
    fieldPoints->GetPointData( i, &gradient );
    for( unsigned int d = 0; d < Dimension; d++ )
      {
      if( !this->m_Directionality[d] )
        {
        gradient[d] = 0.0;
        }
      }
    fieldPoints->SetPointData( i, gradient*gradientScalingFactor );
    }

  if( this->GetBindBoundary() == true )
    {
    typedef Image<char, Dimension> DummyImageType;
    typename DummyImageType::Pointer dummyImage = DummyImageType::New();
    dummyImage->SetRegions( this->m_Size );
    dummyImage->SetOrigin( this->m_Origin );
    dummyImage->SetSpacing( this->m_Spacing );
    dummyImage->Allocate();

    unsigned long count = fieldPoints->GetNumberOfPoints();
    VectorType zeros;
    zeros.Fill( 0 ); 
    
    ImageRegionExclusionConstIteratorWithIndex<DummyImageType> It( 
      dummyImage, dummyImage->GetLargestPossibleRegion() );
    It.SetExclusionRegionToInsetRegion();
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      typename DummyImageType::PointType imagePoint; 
      dummyImage->TransformIndexToPhysicalPoint( It.GetIndex(), imagePoint );
      
      typename BSplinePointSetType::PointType fieldPoint;
      fieldPoint.CastFrom( imagePoint );
      fieldPoints->SetPoint( count, fieldPoint );
      fieldPoints->SetPointData( count, zeros );
      weights->InsertElement( count, 1.0e6 );
      count++;
      }
    }  
      
  typename BSplineFilterType::ArrayType nlevels;
  typename BSplineFilterType::ArrayType ncps;

  nlevels.Fill( 1 );
  for( unsigned int d = 0; d < Dimension; d++ )
    {
    ncps[d] = this->m_TotalDeformationFieldControlPoints
      ->GetLargestPossibleRegion().GetSize()[d];
    }

  switch ( this->m_WhichGradient )
    {
    case 0: default:
      {
      typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();

      bspliner->SetInput( fieldPoints );
      bspliner->SetOrigin( this->m_Origin );
      bspliner->SetSpacing( this->m_Spacing );
      bspliner->SetSize( this->m_Size );
      bspliner->SetNumberOfLevels( nlevels );
      bspliner->SetSplineOrder( this->m_SplineOrder );
      bspliner->SetNumberOfControlPoints( ncps );
      bspliner->SetGenerateOutputImage( false );
      bspliner->SetPointWeights( weights );
      bspliner->Update();

      this->m_GradientFieldControlPoints = bspliner->GetPhiLattice();

      break;
      }
    case 1:
      {
      typename BSplineControlPointFilterType::Pointer bspliner 
        = BSplineControlPointFilterType::New();

      bspliner->SetSplineOrder( this->m_SplineOrder );
      bspliner->SetInput( this->m_TotalDeformationFieldControlPoints );
      bspliner->SetOrigin( this->m_Origin );
      bspliner->SetSize( this->m_Size );
      bspliner->SetSpacing( this->m_Spacing );

      this->m_GradientFieldControlPoints
        = bspliner->CalculateLatticeGradientFromPoints( fieldPoints );

      break;
      }
    }
  return sumValue;
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
void
FFDPointSetRegistrationFilter<TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Maximum iterations = "
     << this->m_MaximumNumberOfIterations << std::endl;
  os << indent << "Number of levels = "
     << this->m_NumberOfLevels << std::endl;
  os << indent << "Line search maximum iterations = "
     << this->m_LineSearchMaximumIterations << std::endl;
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
void
FFDPointSetRegistrationFilter<TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::LineMinimization( RealType *step, RealType *fret )
{
  std::cout << "    Begin line search..." << std::endl;

  // We should now have a, b and c, as well as f(a), f(b), f(c),
  // where b gives the minimum energy position;
  RealType ax, bx, fbx, cx;
  fbx = this->FindBracketingTriplet( &ax, &bx, &cx );

  this->BrentSearch( ax, bx, fbx, cx, step, fret );
//  this->GoldenSectionSearch( ax, bx, cx, step, fret );
 //this->BruteForceSearch( -0.5, 0.0, 1.5, step, fret );
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
typename FFDPointSetRegistrationFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>::RealType
FFDPointSetRegistrationFilter<TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::EvaluateEnergyForLineSearch( RealType lambda )
{
  return this->EvaluateMetric( lambda );
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
typename FFDPointSetRegistrationFilter
  <TFixedPointSet, TMovingPointSet, TWarpedPointSet>::RealType
FFDPointSetRegistrationFilter<TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::FindBracketingTriplet( RealType *ax, RealType *bx, RealType *cx )
{
  const RealType Gold = 1.618034;
  const RealType Glimit = 100.0;
  const RealType Tiny = 1e-20;
  *ax = 0.0;
  *bx = 1.0;
/*
  for( unsigned int i = 0; i < Dimension; i++ )
    {
    *bx *= static_cast<RealType>( 
    this->m_GradientFieldControlPoints
    ->GetLargestPossibleRegion().GetSize()[i] );
    }
  *bx = sqrt( *bx );
*/
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
  if( *cx < *ax )
    {
    std::cout << "      Bracket triple: "
              << "f(" << *cx << ") = " << fc << ", "
              << "f(" << *bx << ") = " << fb << ", "
              << "f(" << *ax << ") = " << fa << std::endl;
    }
  else
    {
    std::cout << "      Bracket triple: "
              << "f(" << *ax << ") = " << fa << ", "
              << "f(" << *bx << ") = " << fb << ", "
              << "f(" << *cx << ") = " << fc << std::endl;
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
        if( *cx < *ax )
          {
          std::cout << "      Bracket triple: "
                    << "f(" << *cx << ") = " << fc << ", "
                    << "f(" << *bx << ") = " << fb << ", "
                    << "f(" << *ax << ") = " << fa << std::endl;
          }
        else
          {
          std::cout << "      Bracket triple: "
                    << "f(" << *ax << ") = " << fa << ", "
                    << "f(" << *bx << ") = " << fb << ", "
                    << "f(" << *cx << ") = " << fc << std::endl;
          }
        return fb;
        }
      else if( fu > fb )
        {
        *cx = u;
        fc = fu;
        if( *cx < *ax )
          {
          std::cout << "      Bracket triple: "
                    << "f(" << *cx << ") = " << fc << ", "
                    << "f(" << *bx << ") = " << fb << ", "
                    << "f(" << *ax << ") = " << fa << std::endl;
          }
        else
          {
          std::cout << "      Bracket triple: "
                    << "f(" << *ax << ") = " << fa << ", "
                    << "f(" << *bx << ") = " << fb << ", "
                    << "f(" << *cx << ") = " << fc << std::endl;
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
    else if( ( u-ulim )*( ulim-*cx ) >=  0.0 )
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
    if( *cx < *ax )
      {
      std::cout << "      Bracket triple: "
        << "f(" << *cx << ") = " << fc << ", "
        << "f(" << *bx << ") = " << fb << ", "
        << "f(" << *ax << ") = " << fa << std::endl;
      }
    else
      {
      std::cout << "      Bracket triple: "
        << "f(" << *ax << ") = " << fa << ", "
        << "f(" << *bx << ") = " << fb << ", "
        << "f(" << *cx << ") = " << fc << std::endl;
      }
    }
  return fb;
}

template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
void
FFDPointSetRegistrationFilter<TFixedPointSet, TMovingPointSet, TWarpedPointSet>
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

  for( unsigned int iter = 1; iter <= this->m_LineSearchMaximumIterations; iter++ )
    {
    std::cout << "        Iteration (Brent's) " << iter
              << ": f(x = " << x << ") = " << fx << ", x in ["
              << a << ", " << b << "] " << std::endl;

    RealType xm = 0.5 * ( a + b );
    RealType tol1 = tol * vnl_math_abs( x ) + ZEPS;
    RealType tol2 = 2.0 * tol1;
    if( vnl_math_abs( x - xm ) <= ( tol2 - 0.5 * ( b - a ) ) )
      {
      std::cout <<  "    Results of line search: E(" << x << ") = " 
        << fx << "." << std::endl;
      itkDebugMacro( << "Results of line search: E(" << x << ") = " 
        << fx << "." );
      *step = x;
      *fret = fx;
      return;
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
  *step = x;
  *fret = fx;
  std::cout <<  "    Results of line search: E(" << x << ") = " << fx << "." << std::endl;
  itkDebugMacro( "    Results of line search: E(" << x << ") = " << fx << "." );
}

/*
template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
void
FFDPointSetRegistrationFilter<TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::GoldenSectionSearch( RealType ax, RealType bx, RealType cx, RealType *step, RealType *fret )
{
  const RealType R = 0.6180339;
  const RealType C = 1.0 - R;
  const RealType tol = 1.0;

  RealType x0 = ax;
  RealType x1;
  RealType x2;
  RealType x3 = cx;
  if( vnl_math_abs( cx-bx ) > vnl_math_abs( bx-ax ) )
    {
    x1 = bx;
    x2 = bx + C*( cx-bx );
    }
  else
    {
    x2 = bx;
    x1 = bx - C*( bx-ax );
    }

  RealType f1 = this->EvaluateEnergyForLineSearch( x1 );
  RealType f2 = this->EvaluateEnergyForLineSearch( x2 );

  unsigned int iters = 0;
  while( iters++ < this->m_LineSearchMaximumIterations &&
          vnl_math_abs( x3-x0 ) > tol*( vnl_math_abs( x1 ) + vnl_math_abs( x2 ) ) )
    {
    std::cout << "        Iteration (Golden Section) " << iters << ": f(" << x1 << ") = " << f1
              << "    f(" << x2 << ") = " << f2 <<  std::endl;
    if( f2 < f1 )
      {
      x0 = x1;
      x1 = x2;
      x2 = R*x1 + C*x3;
      f1 = f2;
      f2 = this->EvaluateEnergyForLineSearch( x2 );
      }
    else
      {
      x3 = x2;
      x2 = x1;
      x1 = R*x2+C*x0;
      f2 = f1;
      f1 = this->EvaluateEnergyForLineSearch( x1 );
      }
    }
  RealType xmin;
  RealType fmin;
  if( f1 < f2 )
    {
    xmin = x1;
    fmin = f1;
    }
  else
    {
    xmin = x2;
    fmin = f2;
    }

  std::cout <<  "    Results of line search: E(" << xmin << ") = " << fmin << "." << std::endl;
  itkDebugMacro( << "Results of line search: E(" << xmin << ") = " << fmin << "." );
  *step = xmin;
  *fret = fmin;
}
*/
/*
template<class TFixedPointSet, class TMovingPointSet, class TWarpedPointSet>
void
FFDPointSetRegistrationFilter<TFixedPointSet, TMovingPointSet, TWarpedPointSet>
::BruteForceSearch( RealType ax, RealType bx, RealType cx, RealType *step, RealType *fret )
{
  static unsigned int whichSearch = 0;
  const unsigned int numberOfIntervals = 10;


  ofstream str;
  if( whichSearch == 0 )
    {
    str.open( "linesearches.txt", std::ios::out );
    }
  else
    {
    str.open( "linesearches.txt", std::ios::app );
    }


  RealType a;
  RealType b;
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

  RealType dx = ( b - a ) / static_cast<RealType>( numberOfIntervals - 1 );

  RealType minStep = 0.0;
  RealType minValue = this->EvaluateEnergyForLineSearch( minStep );
  str << minStep << " " << minValue << " " << whichSearch << " " << this->m_CurrentLevel << std::endl;
  for( RealType x = a; x <= b; x += dx )
    {
    RealType fx = this->EvaluateEnergyForLineSearch( x );
    str << x << " " << fx << " " << whichSearch << " " << this->m_CurrentLevel << std::endl;
    if( fx < minValue )
      {
      minValue = fx;
      minStep = x;
      }
    }

  *step = minStep;
  *fret = minValue;

  whichSearch++;
  std::cout <<  "    Results of brute force line search: E(" << *step << ") = " << *fret << "." << std::endl;
  itkDebugMacro( "    Results of brute force line search: E(" << *step << ") = " << *fret << "." );
}
*/


//    /**
//     * The following code is strictly for visualization purposes.
//     */
//
//    if( this->m_Prolificacy )
//      {
//      SpacingType gridSpacing;
//      for( unsigned int d = 0; d < Dimension; d++ )
//        {
//        unsigned int numberOfGridLines = this->m_InitialMeshResolution[d]
//        * static_cast<unsigned int>( vcl_pow( static_cast<RealType>( 2 ),
//                   static_cast<RealType>( this->m_CurrentLevel ) ) + 1 );
//
//        gridSpacing[d] = ( this->m_MaxBoundary[d] - this->m_MinBoundary[d] )
//        / static_cast<RealType>( numberOfGridLines - 1 ) - 0.0001;
//        }
//
//      typename BSplineControlPointFilterType::Pointer bspliner
//        = BSplineControlPointFilterType::New();
//
//      bspliner->SetSplineOrder( this->m_SplineOrder );
//      bspliner->SetInput( this->m_TotalDeformationFieldControlPoints );
//      bspliner->SetOrigin( this->m_Origin );
//      bspliner->SetSize( this->m_Size );
//      bspliner->SetSpacing( this->m_Spacing );
//
//      OStringStream buf;
//      buf << "_" << this->m_CurrentLevel << "_" << this->m_CurrentIteration;
//
//      {
//      // draw warped grid
//
//      std::string filename = std::string( this->m_FilePrefix )
//        + std::string( "DeformedGrid" ) + buf.str()
//        + std::string( ".txt" );
//
//      ofstream str( filename.c_str() );
//
//      str << "0 0 0 0" << std::endl;
//
//      unsigned int lineCount = 0;
//      for( unsigned int d = 0; d < Dimension; d++ )
//        {
//        RealType delta = this->m_Spacing[d];
//
//        typename BSplineControlPointFilterType::PointType point;
//        for( unsigned int c = 0; c < Dimension; c++ )
//          {
//          point[c] = this->m_MinBoundary[c];
//          }
//
//        while( true )
//          {
//          VectorType vector;
//          bspliner->EvaluateAtPoint( point, vector );
//
//          str << point[0] + vector[0] << " " << point[1] + vector[1];
//          if( Dimension == 2 )
//            {
//            str << " 0 " << lineCount+1 << std::endl;
//            }
//          else
//            {
//            str << " " << point[2] + vector[2] << " " << lineCount+1 << std::endl;
//            }
//
//          point[d] += delta;
//
//
//          unsigned int doneFlag = 0;
//          if( point[d] > this->m_MaxBoundary[d] )
//            {
//
//            lineCount++;
//            point[d] = this->m_MinBoundary[d];
//            doneFlag++;
//
//            for( unsigned int c = 0; c < Dimension-1; c++ )
//              {
//              point[(d+c+1)%Dimension] += gridSpacing[(d+c+1)%Dimension];
//              if( point[(d+c+1)%Dimension] > this->m_MaxBoundary[(d+c+1)%Dimension] )
//                {
//                point[(d+c+1)%Dimension] = this->m_MinBoundary[(d+c+1)%Dimension];
//                doneFlag++;
//                }
//              else
//                {
//                break;
//                }
//              }
//            }
//          if( doneFlag == Dimension )
//            {
//            break;
//            }
//          }
//        }
//      str << "0 0 0 0" << std::endl;
//      str.close();
//      }
//
//      // warp input points
//      std::string filename = std::string( this->m_FilePrefix )
//        + std::string( "WarpedInputPoints" ) + buf.str()
//        + std::string( ".txt" );
//
//      ofstream str( filename.c_str() );
//      str << "0 0 0 0" << std::endl;
//
//    for( unsigned int n = 0; n < this->GetInput( 1 )->GetNumberOfPoints(); n++ )
//      {
//      MovingPointType inputPoint;
//      this->GetInput( 1 )->GetPoint( n, &inputPoint );
//
//      typename BSplineControlPointFilterType::PointType point;
//      for( unsigned int d = 0; d < Dimension; d++ )
//        {
//        point[d] = inputPoint[d];
//        }
//      VectorType vector;
//      bspliner->EvaluateAtPoint( point, vector );
//
//      point += vector;
//
//      str << point[0] << " " << point[1] << " ";
//      if( Dimension == 3 )
//        {
//        str << point[2] << " " << n+1 << std::endl;
//        }
//      else
//        {
//        str << "0 " << n+1 << std::endl;
//        }
//      }
//    str << "0 0 0 0" << std::endl;
//    str.close();
//
//    // Output sample points
//    {
//    std::string filename = std::string( this->m_FilePrefix )
//      + std::string( "WarpedSamplePoints" ) + buf.str()
//      + std::string( ".txt" );
//
//    ofstream str( filename.c_str() );
//    str << "0 0 0 0" << std::endl;
//    for( unsigned int j = 0; j < this->m_WarpedSamplePoints->GetNumberOfPoints(); j++ )
//      {
//      WarpedPointType point;
//      this->m_WarpedSamplePoints->GetPoint( j, &point );
//      str << point[0] << " " << point[1] << " ";
//      if( Dimension == 3 )
//        {
//        str << point[2] << " " << j+1 << std::endl;
//        }
//      else
//        {
//        str << "0 " << j+1 << std::endl;
//        }
//      }
//    str << "0 0 0 0" << std::endl;
//    str.close();
//    }
//
//    }
//    /**
//     * end visualization code
//     */


}
#endif

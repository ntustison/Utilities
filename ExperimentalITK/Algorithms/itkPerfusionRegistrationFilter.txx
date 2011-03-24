/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: itkPerfusionRegistrationFilter.txx,v $
Language:  C++

Date:      $Date: 2008/10/18 00:13:41 $
Version:   $Revision: 1.1.1.1 $

Copyright (c) Insight Software Consortium. All rights reser
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for detail.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef _itkPerfusionRegistrationFilter_txx_
#define _itkPerfusionRegistrationFilter_txx_

// disable debug warnings in MS compiler
#ifdef _MSC_VER
#pragma warning(disable: 4786)
#endif

#include "itkPerfusionRegistrationFilter.h"

#include "itkAddImageFilter.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkCrossCorrelationRegistrationFunction.h"
#include "itkMultiplyImageFilter.h"
#include "itkVectorImageFileReader.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkVectorFieldGradientImageFunction.h"
#include "itkWarpImageFilter.h"

#include "vnl/vnl_math.h"
#include "vnl/algo/vnl_determinant.h"
#include "vnl/algo/vnl_matrix_inverse.h"

namespace itk {

template<class TImage, class TWarpedImage>
PerfusionRegistrationFilter<TImage, TWarpedImage>
::PerfusionRegistrationFilter()
{
  this->SetNumberOfRequiredInputs( 0 );

  // Default values

  this->m_LineSearchMaximumIterations = 0;
  this->m_MaximumNumberOfIterations = 10;
  this->m_SplineOrder = 3;
  this->m_TemporalSplineOrder = 3;
  this->m_TemporalOrigin = 0.0;
  this->m_TemporalEnd = 1.0;
  this->m_WrapTime = false;

  this->m_MeshResolution.Fill( 4 );

  this->m_TotalDeformationFieldControlPoints = NULL;

  // Set up the default metric type (mean squares)
  typedef CrossCorrelationRegistrationFunction
        <ImageType, ImageType, DeformationFieldType> DefaultMetricType;
  typename DefaultMetricType::Pointer DefaultMetric
          = DefaultMetricType::New();
  this->m_PDEDeformableMetric = DefaultMetric;
  this->m_PDEDeformableMetric->SetNormalizeGradient( false );
  this->SetMaximizeMetric( false );

  this->m_PDEDeformableMetric = NULL;

  // Set up the default image interpolator (linear)
  typedef LinearInterpolateImageFunction
          <ImageType, double> DefaultImageInterpolatorType;
  typename DefaultImageInterpolatorType::Pointer interpolator
          = DefaultImageInterpolatorType::New();
  this->m_ImageInterpolator = interpolator;
}

template<class TImage, class TWarpedImage>
PerfusionRegistrationFilter<TImage, TWarpedImage>
::~PerfusionRegistrationFilter()
{
}

template<class TImage, class TWarpedImage>
void
PerfusionRegistrationFilter<TImage, TWarpedImage>
::GenerateData()
{
  for ( this->m_CurrentLevel = 0;
    this->m_CurrentLevel < this->m_MaximumNumberOfIterations.size(); this->m_CurrentLevel++ )
    {

    this->InitializeImages();

    std::cout << "Current Level = " << this->m_CurrentLevel+1 << " (" << this->m_MaximumNumberOfIterations.size() << " total levels).  "
              << "Image size = " << this->m_PyramidLevelImageSizes[this->m_CurrentLevel] << ".  "
              << "Grid Size = " << this->m_TotalDeformationFieldControlPoints->GetLargestPossibleRegion().GetSize() << ". " << std::endl;

    itkDebugMacro( "Current Level = " << this->m_CurrentLevel+1 << " (" << this->m_MaximumNumberOfIterations.size() << " total levels).  "
              << "Image Size = " << this->m_PyramidLevelImageSizes[this->m_CurrentLevel] << ".  "
              << "Grid Size = " << this->m_TotalDeformationFieldControlPoints->GetLargestPossibleRegion().GetSize() << ". " );

    if ( this->m_MaximumNumberOfIterations[this->m_CurrentLevel] > 0 )
      {
      std::cout << "Begin iterative solve. " << std::endl;
      itkDebugMacro( "Begin iterative solve. " );
      this->IterativeSolve();
      std::cout << "Finish iterative solve. " << std::endl;
      itkDebugMacro( "Finish iterative solve. " );
      }
    }
}

template<class TImage, class TWarpedImage>
void
PerfusionRegistrationFilter<TImage, TWarpedImage>
::InitializeImages()
{
  if ( this->m_CurrentLevel == 0 )
    {
    typename ControlPointLatticeType::RegionType::SizeType ncps;
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      ncps[i] = this->m_MeshResolution[i] + this->m_SplineOrder;
      }
    ncps[ImageDimension] = this->m_TemporalMeshResolution + this->m_TemporalSplineOrder;
    if ( this->m_WrapTime )
      {
      ncps[ImageDimension] += this->m_TemporalSplineOrder;
      }

    VectorType V;
    V.Fill( 0.0 );
    this->m_TotalDeformationFieldControlPoints = ControlPointLatticeType::New();
    this->m_TotalDeformationFieldControlPoints->SetRegions( ncps );
    this->m_TotalDeformationFieldControlPoints->Allocate();
    this->m_TotalDeformationFieldControlPoints->FillBuffer( V );

    typedef ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( this->m_ImageFileNames[0].c_str() );
    reader->Update();

    typename ImagePyramidType::Pointer pyramid = ImagePyramidType::New();
    pyramid->SetInput( reader->GetOutput() );
    pyramid->SetNumberOfLevels( this->m_MaximumNumberOfIterations.size() );
    pyramid->SetStartingShrinkFactors( this->m_ImageShrinkFactors.GetDataPointer() );
    pyramid->Update();

    for( unsigned int i = 0; i < pyramid->GetNumberOfOutputs(); i++ )
      {
      this->m_PyramidLevelImageSizes.push_back(
        pyramid->GetOutput( i )->GetLargestPossibleRegion().GetSize() );
      }
    }
  else
    {
    typename ImageType::Pointer referenceImage =
      this->EvaluateImageAtPyramidLevel( 0 );

    // Refine the control lattice

    typedef BSplineControlPointImageFilter
      <ControlPointLatticeType, DeformationFieldType> BSplineControlPointsFilterType;
    typename BSplineControlPointsFilterType::Pointer bspliner = BSplineControlPointsFilterType::New();
    typename BSplineControlPointsFilterType::ArrayType close;

    close.Fill( false );
    bspliner->SetSplineOrder( this->m_SplineOrder );
    bspliner->SetCloseDimension( close );
    bspliner->SetInput( this->m_TotalDeformationFieldControlPoints );
    bspliner->SetOrigin( referenceImage->GetOrigin() );
    bspliner->SetSpacing( referenceImage->GetSpacing() );
    bspliner->SetSize( this->m_PyramidLevelImageSizes[this->m_CurrentLevel] );
    bspliner->Update();

    typename BSplineFilterType::ArrayType nlevels;
    nlevels.Fill( 2 );
    this->m_TotalDeformationFieldControlPoints = bspliner->RefineControlPointLattice( nlevels );
    }
}

template<class TImage, class TWarpedImage>
void
PerfusionRegistrationFilter<TImage, TWarpedImage>
::IterativeSolve()
{
  const RealType ftol = 1.0e-4;
  const RealType eps = 1.0e-10;

  typename ControlPointLatticeType::Pointer G = ControlPointLatticeType::New();
  typename ControlPointLatticeType::Pointer H = ControlPointLatticeType::New();

  if ( this->m_CurrentLevel > 0 )
    {
    this->m_CurrentDeformationFieldControlPoints = ControlPointLatticeType::New();
    this->m_CurrentDeformationFieldControlPoints->SetRegions
      ( this->m_TotalDeformationFieldControlPoints->GetLargestPossibleRegion() );
    this->m_CurrentDeformationFieldControlPoints->Allocate();
    VectorType V;
    V.Fill( 0 );
    this->m_CurrentDeformationFieldControlPoints->FillBuffer( V );
    }

  this->EvaluateGradientFieldOverImageRegion();
  RealType fp = this->EvaluateMetricOverImageRegion();

  typedef ImageDuplicator<ControlPointLatticeType> DuplicatorType;
  typename DuplicatorType::Pointer duplicatorH = DuplicatorType::New();
  duplicatorH->SetInputImage( this->m_GradientFieldControlPoints );
  duplicatorH->Update();
  H = duplicatorH->GetOutput();

  typename DuplicatorType::Pointer duplicatorG = DuplicatorType::New();
  duplicatorG->SetInputImage( this->m_GradientFieldControlPoints );
  duplicatorG->Update();
  G = duplicatorG->GetOutput();

  ImageRegionIterator<ControlPointLatticeType> ItG( G, G->GetLargestPossibleRegion() );
  ImageRegionIterator<ControlPointLatticeType> ItH( H, H->GetLargestPossibleRegion() );
  ImageRegionIterator<ControlPointLatticeType> ItXi( this->m_GradientFieldControlPoints,
      this->m_GradientFieldControlPoints->GetLargestPossibleRegion() );

  ItG.GoToBegin();
  ItH.GoToBegin();
  ItXi.GoToBegin();
  while ( !ItG.IsAtEnd() )
    {
    ItG.Set( -ItG.Get() );
    ItH.Set( -ItH.Get() );
    ItXi.Set( -ItXi.Get() );
    ++ItG;
    ++ItH;
    ++ItXi;
    }

  for ( unsigned its = 1; its <= this->m_MaximumNumberOfIterations[this->m_CurrentLevel]; its++ )
    {
    std::cout << "  Iteration " << its << ": Current Energy = " << fp << std::endl;
    itkDebugMacro( << "Iteration = " << its << ", Current Energy = " << fp );

    RealType gradientStep = 1.0;
    RealType fret;
    if ( this->m_LineSearchMaximumIterations > 0 )
      {
      this->LineMinimization( &gradientStep, &fret );
      }
    else
      {
      fret = this->EvaluateMetricOverImageRegion( gradientStep );
      }

    if ( 2.0*vnl_math_abs( fret - fp ) <= ftol*( vnl_math_abs( fret ) + vnl_math_abs( fp ) + eps ) )
      {
      std::cout << "  Exit condition (normal):  |fret - fp| = " << vnl_math_abs( fret - fp ) << std::endl;
      return;
      }

    fp = fret;

    typename TimeDependentRealImageType::Pointer lambda = TimeDependentRealImageType::New();
    lambda->SetOrigin( this->m_GradientFieldControlPoints->GetOrigin() );
    lambda->SetSpacing( this->m_GradientFieldControlPoints->GetSpacing() );
    lambda->SetRegions( this->m_GradientFieldControlPoints->GetLargestPossibleRegion().GetSize() );
    lambda->Allocate();
    lambda->FillBuffer( gradientStep );

    typedef MultiplyImageFilter<ControlPointLatticeType, TimeDependentRealImageType,
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

    if ( this->m_CurrentLevel > 0 )
      {
      typename AdderType::Pointer adder = AdderType::New();
      adder->SetInput1( multiplier->GetOutput() );
      adder->SetInput2( this->m_CurrentDeformationFieldControlPoints );
      adder->Update();

      this->m_CurrentDeformationFieldControlPoints = adder->GetOutput();
      }

    ImageRegionIterator<ControlPointLatticeType> ItG( G, G->GetLargestPossibleRegion() );
    ImageRegionIterator<ControlPointLatticeType> ItH( H, H->GetLargestPossibleRegion() );
    ImageRegionIterator<ControlPointLatticeType> ItXi( this->m_GradientFieldControlPoints,
        this->m_GradientFieldControlPoints->GetLargestPossibleRegion() );

    RealType gg = 0.0;
    RealType dgg = 0.0;
    for ( ItG.GoToBegin(), ItXi.GoToBegin(); !ItG.IsAtEnd(); ++ItG, ++ItXi )
      {
      VectorType xi = ItXi.Get();
      VectorType g = ItG.Get();
      gg += g * g;
      // Polak-Ribiere
      dgg += ( xi + g ) * xi;
      }
    if ( gg == 0.0 )
      {
      std::cout << "  Exit condition (gradient = 0): " << std::endl;
      return;
      }
    RealType gamma =
      vnl_math_max( static_cast<RealType>( 0.0 ), static_cast<RealType>( dgg / gg ) );

    ItG.GoToBegin();
    ItH.GoToBegin();
    ItXi.GoToBegin();
    while ( !ItG.IsAtEnd() )
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

template<class TImage, class TWarpedImage>
void
PerfusionRegistrationFilter<TImage, TWarpedImage>
::EvaluateGradientFieldOverImageRegion()
{
  const RealType SmallGradientValue = 1.0e-20;

  typename ImageType::Pointer referenceImage =
    this->EvaluateImageAtPyramidLevel( 0 );

//   std::cout << "Evaluating gradient field." << std::endl;
  itkDebugMacro( "Evaluating gradient field." );

  typename DeformationFieldType::Pointer zeros = DeformationFieldType::New();
  zeros->SetOrigin( referenceImage->GetOrigin() );
  zeros->SetSpacing( referenceImage->GetSpacing() );
  zeros->SetRegions( this->m_PyramidLevelImageSizes[this->m_CurrentLevel] );
  zeros->Allocate();
  VectorType V;
  V.Fill( 0 );
  zeros->FillBuffer( V );

  this->m_PDEDeformableMetric->SetRadius( this->m_MetricRadius );

  typename PointSetType::Pointer fieldPoints = PointSetType::New();
  fieldPoints->Initialize();
  typename BSplineFilterType::WeightsContainerType::Pointer weights
    = BSplineFilterType::WeightsContainerType::New();

  typedef WarpImageFilter<ImageType, ImageType, DeformationFieldType> WarperType;


  unsigned long index = 0;
  for ( unsigned int i = 0; i < this->m_TimePoints.size(); i++ )
    {
//     std::cout << "  Evaluating contribution from image " << i << std::endl;
    itkDebugMacro( "Evaluating contribution from image " << i );

    int ip1 = -1;
    int im1 = -1;

    if( i > 0 && i < this->m_TimePoints.size() - 1 )
      {
      im1 = this->m_TimePoints.size() - 1;
      ip1 = i + 1;
      }
    else if( this->m_WrapTime )
      {
      if( i == 0 )
        {
        im1 = this->m_TimePoints.size() - 1;
        }
      if( i > 0 )
        {
        im1 = ( i - 1 );
        }
      ip1 = ( i + 1 ) % this->m_TimePoints.size();
      }
    else
      {
      if( i == 0 )
        {
        ip1 = 1;
        }
      else if( i == this->m_TimePoints.size() - 1 )
        {
        im1 = i - 1;
        }
      }

    typename ImageType::Pointer fixedImage = NULL;
    typename ImageType::Pointer movingImage = NULL;

    // Create the current estimate for the fixed image
      {
      typename DeformationFieldType::Pointer deformationField = DeformationFieldType::New();
      deformationField = this->EvaluateFieldFromControlPointsAtTimePoint
        ( this->m_TotalDeformationFieldControlPoints, i );

      typename WarperType::Pointer warper = WarperType::New();
      warper->SetInput( this->EvaluateImageAtPyramidLevel( i ) );
      warper->SetDeformationField( deformationField );
      warper->SetInterpolator( this->m_ImageInterpolator );
      warper->SetOutputSpacing( deformationField->GetSpacing() );
      warper->SetOutputOrigin( deformationField->GetOrigin() );
      warper->Update();

      fixedImage = warper->GetOutput();
      this->m_PDEDeformableMetric->SetFixedImage( fixedImage );
      }


    // Create the current estimate for the moving image 1
    if( im1 >= 0 )
      {
      typename DeformationFieldType::Pointer deformationField = DeformationFieldType::New();
      deformationField = this->EvaluateFieldFromControlPointsAtTimePoint
        ( this->m_TotalDeformationFieldControlPoints, im1 );

      typename WarperType::Pointer warper = WarperType::New();
      warper->SetInput( this->EvaluateImageAtPyramidLevel( im1 ) );
      warper->SetDeformationField( deformationField );
      warper->SetInterpolator( this->m_ImageInterpolator );
      warper->SetOutputSpacing( deformationField->GetSpacing() );
      warper->SetOutputOrigin( deformationField->GetOrigin() );
      warper->Update();

      movingImage = warper->GetOutput();

      this->m_PDEDeformableMetric->SetFixedImage( fixedImage );
      this->m_PDEDeformableMetric->SetMovingImage( movingImage );
      this->m_PDEDeformableMetric->SetDeformationField( zeros );
      this->m_PDEDeformableMetric->InitializeIteration();

      typedef typename NeighborhoodAlgorithm
        ::ImageBoundaryFacesCalculator<DeformationFieldType> FaceCalculatorType;
      FaceCalculatorType faceCalculator;

      typename FaceCalculatorType::FaceListType faceList = faceCalculator( deformationField,
        deformationField->GetLargestPossibleRegion(), this->m_MetricRadius );
      typename FaceCalculatorType::FaceListType::iterator fit;

      for ( fit = faceList.begin(); fit != faceList.end(); ++fit )
        {
        NeighborhoodIteratorType It( this->m_MetricRadius, deformationField, *fit );
        for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
          {
          this->m_PDEDeformableMetric->SetEnergy( 0.0 );
          VectorType grad = this->m_PDEDeformableMetric->ComputeUpdate( It, NULL );
          if ( vnl_math_isnan( this->m_PDEDeformableMetric->GetEnergy() ) )
            {
            continue;
            }
          else if ( grad.GetSquaredNorm() >= SmallGradientValue )
            {
            if ( this->m_MaximizeMetric )
              {
              grad *= -1.0;
              }
            typename ImageType::PointType pt;
            fixedImage->TransformIndexToPhysicalPoint( It.GetIndex(), pt );
            typename PointSetType::PointType point;
            typename PointSetType::PixelType data;
            for ( unsigned int d = 0; d < ImageDimension; d++ )
              {
              point[d] = pt[d];
              data[d] = grad[d];
              }
            point[ImageDimension] = im1;
            fieldPoints->SetPoint( index, point );
            fieldPoints->SetPointData( index, data );
            index++;
            }
          }
        }
      }

    // Create the current estimate for the moving image 1
    if( ip1 >= 0 )
      {
      typename DeformationFieldType::Pointer deformationField = DeformationFieldType::New();
      deformationField = this->EvaluateFieldFromControlPointsAtTimePoint
        ( this->m_TotalDeformationFieldControlPoints, ip1 );

      typename WarperType::Pointer warper = WarperType::New();
      warper->SetInput( this->EvaluateImageAtPyramidLevel( ip1 ) );
      warper->SetDeformationField( deformationField );
      warper->SetInterpolator( this->m_ImageInterpolator );
      warper->SetOutputSpacing( deformationField->GetSpacing() );
      warper->SetOutputOrigin( deformationField->GetOrigin() );
      warper->Update();

      movingImage = warper->GetOutput();

      this->m_PDEDeformableMetric->SetFixedImage( fixedImage );
      this->m_PDEDeformableMetric->SetMovingImage( movingImage );
      this->m_PDEDeformableMetric->SetDeformationField( zeros );
      this->m_PDEDeformableMetric->InitializeIteration();

      typedef typename NeighborhoodAlgorithm
        ::ImageBoundaryFacesCalculator<DeformationFieldType> FaceCalculatorType;
      FaceCalculatorType faceCalculator;

      typename FaceCalculatorType::FaceListType faceList = faceCalculator( deformationField,
        deformationField->GetLargestPossibleRegion(), this->m_MetricRadius );
      typename FaceCalculatorType::FaceListType::iterator fit;

      for ( fit = faceList.begin(); fit != faceList.end(); ++fit )
        {
        NeighborhoodIteratorType It( this->m_MetricRadius, deformationField, *fit );
        for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
          {
          this->m_PDEDeformableMetric->SetEnergy( 0.0 );
          VectorType grad = this->m_PDEDeformableMetric->ComputeUpdate( It, NULL );
          if ( vnl_math_isnan( this->m_PDEDeformableMetric->GetEnergy() ) )
            {
            continue;
            }
          else if ( grad.GetSquaredNorm() >= SmallGradientValue )
            {
            if ( this->m_MaximizeMetric )
              {
              grad *= -1.0;
              }
            typename ImageType::PointType pt;
            fixedImage->TransformIndexToPhysicalPoint( It.GetIndex(), pt );
            typename PointSetType::PointType point;
            typename PointSetType::PixelType data;
            for ( unsigned int d = 0; d < ImageDimension; d++ )
              {
              point[d] = pt[d];
              data[d] = grad[d];
              }
            point[ImageDimension] = this->m_TimePoints[ip1];
            fieldPoints->SetPoint( index, point );
            fieldPoints->SetPointData( index, data );
            index++;
            }
          }
        }
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
  std::cout << "Normalizing gradient values." << std::endl;

  RealType sumSquaredNorm = 0.0;
  for ( unsigned int i = 0; i < fieldPoints->GetNumberOfPoints(); i++ )
    {
    typename PointSetType::PixelType data;
    fieldPoints->GetPointData( i, &data );
    sumSquaredNorm += ( data ).GetSquaredNorm();
    }
  RealType sigma;
  for ( unsigned int d = 0; d < ImageDimension; d++ )
    {
    V[d] = referenceImage->GetSpacing()[d];
    }
  if ( ImageDimension == 2 )
    {
    sigma = V.GetNorm();
    }
  else if ( ImageDimension == 3 )
    {
    sigma = V.GetNorm()/sqrt( 2.0 );
    }

  RealType gradientScalingFactor = 0.0;
  if ( sumSquaredNorm > 0.0 )
    {
    gradientScalingFactor = sigma*sqrt
      ( static_cast<RealType>( ImageDimension * fieldPoints->GetNumberOfPoints() )
        / sumSquaredNorm );
    }
  std::cout << "Gradient scaling factor = " << gradientScalingFactor << std::endl;

  for ( unsigned int i = 0; i < fieldPoints->GetNumberOfPoints(); i++ )
    {
    index = i;
    typename PointSetType::PixelType data;
    fieldPoints->GetPointData( i, &data );
    fieldPoints->SetPointData( i, gradientScalingFactor * data );
    weights->InsertElement( i, 1.0 );
    }

  itkDebugMacro( "Fitting B-spline field to gradient values." );

  // Define the rest of the parameters for the B-spline filter

  typename TimeDependentFieldType::PointType origin;
  typename TimeDependentFieldType::SpacingType spacing;
  typename TimeDependentFieldType::SizeType size;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    origin[i] = referenceImage->GetOrigin()[i];
    spacing[i] = referenceImage->GetSpacing()[i];
    size[i] = referenceImage->GetLargestPossibleRegion().GetSize()[i];
    }
  origin[ImageDimension] = this->m_TemporalOrigin;
  spacing[ImageDimension] = ( this->m_TemporalEnd - this->m_TemporalOrigin )
      / static_cast<RealType>( this->m_TimePoints.size() );
  // The size doesn't affect the resolution as does the ncps
  size[ImageDimension] = this->m_TimePoints.size() + 1;

  typename BSplineFilterType::ArrayType close;
  close.Fill( false );
  close[ImageDimension] = this->m_WrapTime;
  typename BSplineFilterType::ArrayType nlevels;
  typename BSplineFilterType::ArrayType ncps;

  nlevels.Fill( 1 );
  for ( unsigned int i = 0; i < ImageDimension + 1; i++ )
    {
    ncps[i] = this->m_TotalDeformationFieldControlPoints->GetLargestPossibleRegion().GetSize()[i];
    }
  if ( this->m_WrapTime )
    {
    ncps[ImageDimension] += this->m_TemporalSplineOrder;
    }

  typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
  bspliner->SetInput( fieldPoints );
  bspliner->SetPointWeights( weights );
  bspliner->SetOrigin( origin );
  bspliner->SetSpacing( spacing );
  bspliner->SetSize( size );
  bspliner->SetNumberOfLevels( nlevels );
  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetNumberOfControlPoints( ncps );
  bspliner->SetCloseDimension( close );
  bspliner->SetGenerateOutputImage( false );
  bspliner->Update();

  this->m_GradientFieldControlPoints = bspliner->GetPhiLattice();
}

template<class TImage, class TWarpedImage>
typename PerfusionRegistrationFilter<TImage, TWarpedImage>::RealType
PerfusionRegistrationFilter<TImage, TWarpedImage>
::EvaluateMetricOverImageRegion( RealType t = 0 )
{
  typename ImageType::Pointer referenceImage =
    this->EvaluateImageAtPyramidLevel( 0 );

  typename TimeDependentRealImageType::Pointer lambda = TimeDependentRealImageType::New();
  lambda->SetOrigin( this->m_GradientFieldControlPoints->GetOrigin() );
  lambda->SetSpacing( this->m_GradientFieldControlPoints->GetSpacing() );
  lambda->SetRegions( this->m_GradientFieldControlPoints->GetLargestPossibleRegion().GetSize() );
  lambda->Allocate();
  lambda->FillBuffer( t );

  typedef MultiplyImageFilter<ControlPointLatticeType, TimeDependentRealImageType,
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

  VectorType V;
  V.Fill( 0 );
  typename DeformationFieldType::Pointer zeros = DeformationFieldType::New();
  zeros->SetOrigin( referenceImage->GetOrigin() );
  zeros->SetSpacing( referenceImage->GetSpacing() );
  zeros->SetRegions( this->m_PyramidLevelImageSizes[this->m_CurrentLevel] );
  zeros->Allocate();
  zeros->FillBuffer( V );

  RealType metricEnergy = 0.0;
  RealType metricCount = 0.0;

  this->m_PDEDeformableMetric->SetRadius( this->m_MetricRadius );

//   std::cout << "Evaluating metric field." << std::endl;;

  this->m_PDEDeformableMetric->SetRadius( this->m_MetricRadius );

  typename PointSetType::Pointer fieldPoints = PointSetType::New();
  fieldPoints->Initialize();
  typename BSplineFilterType::WeightsContainerType::Pointer weights
    = BSplineFilterType::WeightsContainerType::New();

  typedef WarpImageFilter<ImageType, ImageType, DeformationFieldType> WarperType;

  for ( unsigned int i = 0; i < this->m_TimePoints.size(); i++ )
    {
//     std::cout << "   Evaluating contribution from image " << i << std::endl;
    itkDebugMacro( "Evaluating contribution from image " << i );

    int ip1 = -1;
    int im1 = -1;

    if( i > 0 && i < this->m_TimePoints.size() - 1 )
      {
      im1 = i - 1;
      ip1 = i + 1;
      }
    else if( this->m_WrapTime )
      {
      if( i == 0 )
        {
        im1 = this->m_TimePoints.size() - 1;
        }
      if( i > 0 )
        {
        im1 = ( i - 1 );
        }
      ip1 = ( i + 1 ) % this->m_TimePoints.size();
      }
    else
      {
      if( i == 0 )
        {
        ip1 = 1;
        }
      else if( i == this->m_TimePoints.size() - 1 )
        {
        im1 = i - 1;
        }
      }

    typename ImageType::Pointer fixedImage = NULL;
    typename ImageType::Pointer movingImage = NULL;

    // Create the current estimate for the fixed image
      {
      typename DeformationFieldType::Pointer deformationField = DeformationFieldType::New();
      deformationField = this->EvaluateFieldFromControlPointsAtTimePoint
        ( this->m_TotalDeformationFieldControlPoints, i );

      typename WarperType::Pointer warper = WarperType::New();
      warper->SetInput( this->EvaluateImageAtPyramidLevel( i ) );
      warper->SetDeformationField( deformationField );
      warper->SetInterpolator( this->m_ImageInterpolator );
      warper->SetOutputSpacing( deformationField->GetSpacing() );
      warper->SetOutputOrigin( deformationField->GetOrigin() );
      warper->Update();

      fixedImage = warper->GetOutput();
      this->m_PDEDeformableMetric->SetFixedImage( fixedImage );
      }


    // Create the current estimate for the moving image 1
    if( im1 >= 0 )
      {
      typename DeformationFieldType::Pointer deformationField = DeformationFieldType::New();
      deformationField = this->EvaluateFieldFromControlPointsAtTimePoint
        ( this->m_TotalDeformationFieldControlPoints, im1 );

      typename WarperType::Pointer warper = WarperType::New();
      warper->SetInput( this->EvaluateImageAtPyramidLevel( im1 ) );
      warper->SetDeformationField( deformationField );
      warper->SetInterpolator( this->m_ImageInterpolator );
      warper->SetOutputSpacing( deformationField->GetSpacing() );
      warper->SetOutputOrigin( deformationField->GetOrigin() );
      warper->Update();

      movingImage = warper->GetOutput();

      this->m_PDEDeformableMetric->SetFixedImage( fixedImage );
      this->m_PDEDeformableMetric->SetMovingImage( movingImage );
      this->m_PDEDeformableMetric->SetDeformationField( zeros );
      this->m_PDEDeformableMetric->InitializeIteration();

      typedef typename NeighborhoodAlgorithm
        ::ImageBoundaryFacesCalculator<DeformationFieldType> FaceCalculatorType;
      FaceCalculatorType faceCalculator;

      typename FaceCalculatorType::FaceListType faceList = faceCalculator( deformationField,
        deformationField->GetLargestPossibleRegion(), this->m_MetricRadius );
      typename FaceCalculatorType::FaceListType::iterator fit;

      for ( fit = faceList.begin(); fit != faceList.end(); ++fit )
        {
        NeighborhoodIteratorType It( this->m_MetricRadius, deformationField, *fit );
        for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
          {
          this->m_PDEDeformableMetric->SetEnergy( 0.0 );
          VectorType grad = this->m_PDEDeformableMetric->ComputeUpdate( It, NULL );
          if ( vnl_math_isnan( this->m_PDEDeformableMetric->GetEnergy() ) )
            {
            continue;
            }
          else
            {
            metricEnergy += this->m_PDEDeformableMetric->GetEnergy();
            metricCount++;
            }
          }
        }
      }

    // Create the current estimate for the moving image 1
    if( ip1 >= 0 )
      {
      typename DeformationFieldType::Pointer deformationField = DeformationFieldType::New();
      deformationField = this->EvaluateFieldFromControlPointsAtTimePoint
        ( this->m_TotalDeformationFieldControlPoints, ip1 );

      typename WarperType::Pointer warper = WarperType::New();
      warper->SetInput( this->EvaluateImageAtPyramidLevel( ip1 ) );
      warper->SetDeformationField( deformationField );
      warper->SetInterpolator( this->m_ImageInterpolator );
      warper->SetOutputSpacing( deformationField->GetSpacing() );
      warper->SetOutputOrigin( deformationField->GetOrigin() );
      warper->Update();

      movingImage = warper->GetOutput();

      this->m_PDEDeformableMetric->SetFixedImage( fixedImage );
      this->m_PDEDeformableMetric->SetMovingImage( movingImage );
      this->m_PDEDeformableMetric->SetDeformationField( zeros );
      this->m_PDEDeformableMetric->InitializeIteration();

      typedef typename NeighborhoodAlgorithm
        ::ImageBoundaryFacesCalculator<DeformationFieldType> FaceCalculatorType;
      FaceCalculatorType faceCalculator;

      typename FaceCalculatorType::FaceListType faceList = faceCalculator( deformationField,
        deformationField->GetLargestPossibleRegion(), this->m_MetricRadius );
      typename FaceCalculatorType::FaceListType::iterator fit;

      for ( fit = faceList.begin(); fit != faceList.end(); ++fit )
        {
        NeighborhoodIteratorType It( this->m_MetricRadius, deformationField, *fit );
        for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
          {
          this->m_PDEDeformableMetric->SetEnergy( 0.0 );
          VectorType grad = this->m_PDEDeformableMetric->ComputeUpdate( It, NULL );
          if ( vnl_math_isnan( this->m_PDEDeformableMetric->GetEnergy() ) )
            {
            continue;
            }
          else
            {
            metricEnergy += this->m_PDEDeformableMetric->GetEnergy();
            metricCount++;
            }
          }
        }
      }
    }


  RealType energy = 0.0;
  if ( metricCount > 0.0 )
    {
    energy += metricEnergy / metricCount;
    }

  return energy;
}

template<class TImage, class TWarpedImage>
typename PerfusionRegistrationFilter<TImage, TWarpedImage>
::DeformationFieldType::Pointer
PerfusionRegistrationFilter<TImage, TWarpedImage>
::EvaluateFieldFromControlPointsAtTimePoint
  ( typename ControlPointLatticeType::Pointer phiLattice, unsigned int which, bool refineLattice = false )
{
  itkDebugMacro( "Evaluating deformation field at time point " << which );

  typename ImageType::Pointer referenceImage = this->EvaluateImageAtPyramidLevel( which );

  typedef BSplineControlPointImageFilter<ControlPointLatticeType, TimeDependentFieldType>
    BSplineControlPointsFilterType;
  typename BSplineControlPointsFilterType::Pointer bspliner = BSplineControlPointsFilterType::New();
  typename BSplineControlPointsFilterType::ArrayType close;
  typename BSplineControlPointsFilterType::ArrayType nlevels;
  typename BSplineControlPointsFilterType::DirectionType direction;
  typename BSplineFilterType::ArrayType order;

  typename TimeDependentFieldType::PointType origin;
  typename TimeDependentFieldType::SpacingType spacing;
  typename TimeDependentFieldType::SizeType size;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    origin[i] = referenceImage->GetOrigin()[i];
    spacing[i] = referenceImage->GetSpacing()[i];
    size[i] = referenceImage->GetLargestPossibleRegion().GetSize()[i];
    }
  origin[ImageDimension] = this->m_TemporalOrigin;
  spacing[ImageDimension] = ( this->m_TemporalEnd - this->m_TemporalOrigin )
      / static_cast<RealType>( this->m_TimePoints.size() );
  direction.SetIdentity();

  // The size doesn't affect the resolution as does the ncps
  size[ImageDimension] = this->m_TimePoints.size() + 1;
  close.Fill( false );
  close[ImageDimension] = this->m_WrapTime;
  order.Fill( this->m_SplineOrder );
  order[ImageDimension] = this->m_TemporalSplineOrder;

  bspliner->SetSplineOrder( order );
  bspliner->SetCloseDimension( close );
  bspliner->SetDirection( direction );
  bspliner->SetInput( phiLattice );
  bspliner->SetOrigin( origin );
  bspliner->SetSpacing( spacing );
  bspliner->SetSize( size );
  bspliner->Update();

  if ( refineLattice )
    {
    nlevels.Fill( 2 );
    this->m_TotalDeformationFieldControlPoints = bspliner->RefineControlPointLattice( nlevels );
    return NULL;
    }
  else
    {
    typename ControlPointLatticeType::Pointer tPhiLattice =
      ControlPointLatticeType::New();
    typename ControlPointLatticeType::SizeType tPhiSize;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      tPhiSize[d] = phiLattice->GetLargestPossibleRegion().GetSize()[d];
      }
    tPhiSize[ImageDimension] = 1;
    tPhiLattice->SetRegions( tPhiSize );
    tPhiLattice->Allocate();

    RealType numberOfSpans =
      phiLattice->GetLargestPossibleRegion().GetSize()[ImageDimension];
    if( this->m_WrapTime )
      {
      numberOfSpans -= this->m_TemporalSplineOrder;
      }
    RealType t = ( this->m_TimePoints[which] - this->m_TemporalOrigin ) /
      ( this->m_TemporalEnd - this->m_TemporalOrigin ) * numberOfSpans;

    bspliner->CollapsePhiLattice( phiLattice, tPhiLattice, t,
      ImageDimension );

    typedef ExtractImageFilter<TimeDependentFieldType, DeformationFieldType> ExtracterType;
    typename ExtracterType::Pointer extracter = ExtracterType::New();
    extracter->SetInput( tPhiLattice );

    typename ExtracterType::InputImageRegionType region;
    typename ExtracterType::InputImageSizeType regionSize;
    typename ExtracterType::InputImageIndexType regionIndex;

    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      regionSize[d] = tPhiLattice->GetLargestPossibleRegion().GetSize()[d];
      regionIndex[d] = tPhiLattice->GetLargestPossibleRegion().GetIndex()[d];
      }
    regionSize[ImageDimension] = 0;
    regionIndex[ImageDimension] = 0;
    region.SetSize( regionSize );
    region.SetIndex( regionIndex );

    extracter->SetExtractionRegion( region );
    extracter->Update();

    typename DeformationFieldType::PointType origin2;
    typename DeformationFieldType::SpacingType spacing2;
    typename DeformationFieldType::SizeType size2;
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      origin2[i] = referenceImage->GetOrigin()[i];
      spacing2[i] = referenceImage->GetSpacing()[i];
      size2[i] = referenceImage->GetLargestPossibleRegion().GetSize()[i];
      }

    typedef BSplineControlPointImageFilter<DeformationFieldType, DeformationFieldType>
      BSplineControlPointsFilterType2;
    typename BSplineControlPointsFilterType2::Pointer bspliner2 =
      BSplineControlPointsFilterType2::New();
    typename BSplineControlPointsFilterType2::ArrayType close2;
    typename BSplineControlPointsFilterType2::ArrayType order2;
    typename BSplineControlPointsFilterType2::DirectionType direction2;
    close2.Fill( false );
    order2.Fill( this->m_SplineOrder );
    direction2.SetIdentity();

    bspliner2->SetSplineOrder( order2 );
    bspliner2->SetCloseDimension( close2 );
    bspliner2->SetInput( extracter->GetOutput() );
    bspliner2->SetOrigin( origin2 );
    bspliner2->SetSpacing( spacing2 );
    bspliner2->SetDirection( direction2 );
    bspliner2->SetSize( size2 );
    bspliner2->Update();

    return bspliner2->GetOutput();

    }
}

template<class TImage, class TWarpedImage>
typename PerfusionRegistrationFilter<TImage, TWarpedImage>
::ImageType::Pointer
PerfusionRegistrationFilter<TImage, TWarpedImage>
::EvaluateImageAtPyramidLevel( unsigned int which )
{
  typedef ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( this->m_ImageFileNames[which].c_str() );
  reader->Update();

  typename ImagePyramidType::Pointer pyramid = ImagePyramidType::New();
  pyramid->SetInput( reader->GetOutput() );
  pyramid->SetNumberOfLevels( this->m_MaximumNumberOfIterations.size() );
  pyramid->SetStartingShrinkFactors( this->m_ImageShrinkFactors.GetDataPointer() );
  pyramid->Update();

  return pyramid->GetOutput( this->m_CurrentLevel );
}

template<class TImage, class TWarpedImage>
void
PerfusionRegistrationFilter<TImage, TWarpedImage>
::GenerateOutputAtTimePoint( unsigned int which )
{
  RealType t = this->m_TimePoints[which];

  typename DeformationFieldType::Pointer deformationField =
    this->EvaluateFieldFromControlPointsAtTimePoint
    ( this->m_TotalDeformationFieldControlPoints, which );

  typedef WarpImageFilter<ImageType,
                          WarpedImageType,
                          DeformationFieldType> WarperType;
  typename WarperType::Pointer warper = WarperType::New();

  typename ImageType::Pointer image = this->EvaluateImageAtPyramidLevel( which );

  warper->SetInput( image );
  warper->SetDeformationField( deformationField );
  warper->SetInterpolator( this->m_ImageInterpolator );
  warper->SetOutputSpacing( image->GetSpacing() );
  warper->SetOutputOrigin( image->GetOrigin() );
  warper->Update();
  this->GraftNthOutput( 0, warper->GetOutput() );
}

template<class TImage, class TWarpedImage>
void PerfusionRegistrationFilter<TImage, TWarpedImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Maximum iterations = "
     << this->m_MaximumNumberOfIterations << std::endl;
  os << indent << "Number of levels = "
     << this->m_MaximumNumberOfIterations.size() << std::endl;
  os << indent << "Line search maximum iterations = "
     << this->m_LineSearchMaximumIterations << std::endl;
}

template<class TImage, class TWarpedImage>
void
PerfusionRegistrationFilter<TImage, TWarpedImage>
::LineMinimization( RealType *step, RealType *fret )
{
  std::cout << "    Begin line search..." << std::endl;

  // We should now have a, b and c, as well as f(a), f(b), f(c),
  // where b gives the minimum energy position;
  RealType ax, bx, cx;
  this->FindBracketingTriplet( &ax, &bx, &cx );

  this->BrentSearch( ax, bx, cx, step, fret );
//  this->GoldenSectionSearch( ax, bx, cx, step, fret );
}

template<class TImage, class TWarpedImage>
typename PerfusionRegistrationFilter<TImage, TWarpedImage>::RealType
PerfusionRegistrationFilter<TImage, TWarpedImage>
::EvaluateEnergyForLineSearch( RealType lambda )
{
  return this->EvaluateMetricOverImageRegion( lambda );
}

template<class TImage, class TWarpedImage>
void
PerfusionRegistrationFilter<TImage, TWarpedImage>
::FindBracketingTriplet( RealType *ax, RealType *bx, RealType *cx )
{
  const RealType Gold = 1.618034;
  const RealType Glimit = 100.0;
  const RealType Tiny = 1e-20;
  *ax = 0.0;
  *bx = 1.0;
  RealType fa = this->EvaluateEnergyForLineSearch( *ax );
  RealType fb = this->EvaluateEnergyForLineSearch( *bx );

  RealType dum;
  if ( fb > fa )
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
  if ( *cx < *ax )
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
  unsigned int iter = 0;
  while ( fb > fc )
    {
    r = ( *bx-*ax )*( fb-fc );
    q = ( *bx-*cx )*( fb-fa );
    RealType denom = 2.0*vnl_math_max( vnl_math_abs( q-r ), Tiny );
    if ( q-r < 0.0 )
      {
      denom = -denom;
      }
    u = *bx - ( ( *bx-*cx )*q - ( *bx-*ax )*r )/denom;
    ulim = *bx + Glimit*( *cx-*bx );
    if ( ( *bx-u )*( u-*cx ) > 0.0 )
      {
      fu = this->EvaluateEnergyForLineSearch( u );
      if ( fu < fc )
        {
        *ax = *bx;
        *bx = u;
        fa = fb;
        fb = fu;
        if ( *cx < *ax )
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
        return;
        }
      else if ( fu > fb )
        {
        *cx = u;
        fc = fu;
        if ( *cx < *ax )
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
        return;
        }
      u = *cx + Gold*( *cx-*bx );
      fu = this->EvaluateEnergyForLineSearch( u );
      }
    else if ( ( *cx-u )*( u-ulim ) > 0.0 )
      {
      fu = this->EvaluateEnergyForLineSearch( u );
      if ( fu < fc )
        {
        *bx = *cx;
        *cx = u;
        u = *cx + Gold*( *cx-*bx );
        fb = fc;
        fc = fu;
        fu = this->EvaluateEnergyForLineSearch( u );
        }
      }
    else if ( ( u-ulim )*( ulim-*cx ) >=  0.0 )
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
    if ( *cx < *ax )
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
}

template<class TImage, class TWarpedImage>
void
PerfusionRegistrationFilter<TImage, TWarpedImage>
::BrentSearch( RealType ax, RealType bx, RealType cx, RealType *step, RealType *fret )
{
  const unsigned int ITMAX = 100;
  const RealType R = 0.6180339;
  const RealType CGOLD = 1.0 - R;
  const RealType tol = 0.1;
  const RealType ZEPS = 1.0e-10;

  RealType a = 0;
  RealType b = 0;
  RealType d = 0;
  RealType p = 0;
  RealType q = 0;
  RealType r = 0;
  RealType etemp = 0;
  RealType e = 0.0;

  if ( ax < cx )
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
  fx = fw = fv = this->EvaluateEnergyForLineSearch( x );

  for ( unsigned int iter = 1; iter <= 100; iter++ )
    {
    std::cout << "        Iteration (Brent's) " << iter
              << ": f(x = " << x << ") = " << fx << ", x in ["
              << a << ", " << b << "] " << std::endl;

    RealType xm = 0.5 * ( a + b );
    RealType tol1 = tol * vnl_math_abs( x ) + ZEPS;
    RealType tol2 = 2.0 * tol1;
    if ( vnl_math_abs( x - xm ) <= ( tol2 - 0.5 * ( b - a ) ) ||
         ( vnl_math_abs( a ) < 0.01 && vnl_math_abs( b ) < 0.01 ) )
      {
      std::cout <<  "    Results of line search: E(" << x << ") = " << fx << "." << std::endl;
      itkDebugMacro( << "Results of line search: E(" << x << ") = " << fx << "." );
      *step = x;
      *fret = fx;
      return;
      }
    if ( vnl_math_abs( e ) > tol1 )
      {
      r = ( x - w ) * ( fx - fv );
      q = ( x - v ) * ( fx - fw );
      p = ( x - v ) * q - ( x - w ) * r;
      q = 2.0 * ( q - r );
      if ( q > 0.0 )
        {
        p = -p;
        }
      q = vnl_math_abs( q );
      etemp = e;
      e = d;
      if ( vnl_math_abs( p ) >= vnl_math_abs( 0.5 * q * etemp )
           || p <= q * ( a - x ) || p >= q * ( b - x ) )
        {
        if ( x >= xm )
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
        if ( u - a < tol2 || b - u < tol2 )
          {
          d = vnl_math_abs( tol1 );
          if ( xm - x <= 0 )
            {
            d = -d;
            }
          }
        }
      }
    else
      {
      if ( x >= xm )
        {
        e = a - x;
        }
      else
        {
        e = b - x;
        }
      d = CGOLD * e;
      }
    if ( vnl_math_abs( d ) >= tol1 )
      {
      u = x + d;
      }
    else
      {
      u = x + vnl_math_abs( tol1 );
      if ( d <= 0 )
        {
        u = x - vnl_math_abs( tol1 );
        }
      }
    fu = this->EvaluateEnergyForLineSearch( u );
    if ( fu <= fx )
      {
      if ( u >= x )
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
      if ( u < x )
        {
        a = u;
        }
      else
        {
        b = u;
        }
      if ( fu <= fw || w == x )
        {
        v = w;
        w = u;
        fv = fw;
        fw = fu;
        }
      else if ( fu <= fv || v == x || v == w )
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
template<class TImage, class TWarpedImage>
void
PerfusionRegistrationFilter<TImage, TWarpedImage>
::GoldenSectionSearch( RealType ax, RealType bx, RealType cx, RealType *step, RealType *fret )
{
  const RealType R = 0.6180339;
  const RealType C = 1.0 - R;
  const RealType tol = 1.0;

  RealType x0 = ax;
  RealType x1;
  RealType x2;
  RealType x3 = cx;
  if ( vnl_math_abs( cx-bx ) > vnl_math_abs( bx-ax ) )
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
  while ( iters++ < this->m_LineSearchMaximumIterations &&
          vnl_math_abs( x3-x0 ) > tol*( vnl_math_abs( x1 ) + vnl_math_abs( x2 ) ) )
    {
    std::cout << "        Iteration (Golden Section) " << iters << ": f(" << x1 << ") = " << f1
              << "    f(" << x2 << ") = " << f2 <<  std::endl;
    if ( f2 < f1 )
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
  if ( f1 < f2 )
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


template<class TMovingImage, class TFixedImage, class TWarpedImage>
void
FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::BruteForceSearch( RealType ax, RealType bx, RealType cx, RealType *step, RealType *fret )
{
  static unsigned int whichSearch = 0;
  const unsigned int numberOfIntervals = 100;


  ofstream str;
  if ( whichSearch == 0 )
    {
    str.open( "linesearches.txt", std::ios::out );
    }
  else
    {
    str.open( "linesearches.txt", std::ios::app );
    }


  RealType a;
  RealType b;
  if ( ax < cx )
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
  str << minStep << " " << minValue << " " << whichSearch << " " << whichSearch << std::endl;
  for ( RealType x = a; x <= b; x += dx )
    {
    RealType fx = this->EvaluateEnergyForLineSearch( x );
    str << x << " " << fx << " " << whichSearch << " " << whichSearch << std::endl;
    if ( fx < minValue )
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

}
#endif

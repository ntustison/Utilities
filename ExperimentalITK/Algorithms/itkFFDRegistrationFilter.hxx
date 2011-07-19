/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: itkFFDRegistrationFilter.hxx,v $
Language:  C++

Date:      $Date: 2008/11/20 19:04:37 $
Version:   $Revision: 1.6 $

Copyright (c) Insight Software Consortium. All rights reser
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for detail.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef _itkFFDRegistrationFilter_hxx_
#define _itkFFDRegistrationFilter_hxx_

// disable debug warnings in MS compiler
#ifdef _MSC_VER
#pragma warning(disable: 4786)
#endif

#include "itkFFDRegistrationFilter.h"

#include "itkAddImageFilter.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkDecomposeTensorFunction.h"
#include "itkDerivativeImageFilter.h"
#include "itkGridImageSource.h"
#include "itkImageDuplicator.h"
#include "itkImageRandomIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMeanSquareRegistrationFunction.h"
#include "itkMultiplyImageFilter.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkTimeProbe.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkVectorFieldGradientImageFunction.h"
#include "itkVectorImageFileWriter.h"
#include "itkWarpImageFilter.h"

#include "vnl/vnl_math.h"
#include "vnl/algo/vnl_determinant.h"
#include "vnl/algo/vnl_matrix_inverse.h"

#include "fstream.h"

namespace itk {

template<class TMovingImage, class TFixedImage, class TWarpedImage>
FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::FFDRegistrationFilter()
{
  this->m_FixedLandmarkContainer = LandmarkContainerType::New();
  this->m_FixedLandmarkContainer->Initialize();
  this->m_MovingLandmarkContainer = LandmarkContainerType::New();
  this->m_MovingLandmarkContainer->Initialize();

  // Default values
  this->SetNumberOfLevels( 3 );
  this->SetLineSearchMaximumIterations( 10 );
  this->SetMaximumNumberOfIterations( 10 );
  this->SetSplineOrder( 3 );
  this->SetLandmarkWeighting( 1.0 );
  this->SetInitializeWithLandmarks( false );

  this->SetEnforceDiffeomorphism( false );
  this->SetEmploySteepestDescent( false );
  this->SetWhichGradient( 0 );
  this->SetMinimumJacobian( 0.5 );
  this->SetMaximumJacobian( 1.5 );
  this->SetInteriorPenaltyParameter( 1 );

  this->m_MeshResolution.Fill( 4 );
  this->m_DoubleMeshResolutionAtEachLevel = true;
  this->m_Directionality.Fill( 1 );

  this->m_AffineNeighborhoodSamplingDensity = 0.5;
  this->m_NumberOfAffineNeighborhoodSamplesPerIteration = 1;
  this->m_AffineNeighborhoodRadius = 0;

  this->m_WeightImage = NULL;
  this->m_InitialControlPointLattice = NULL;

  this->m_GradientComputationTimes.clear();
  this->m_NumberOfGradientPoints.clear();
  this->m_EnergyValues.clear();
  this->m_LevelNumbers.clear();

  this->m_FilePrefix = "output";
  this->m_Prolificacy = true;


  // Set up the default metric type (mean squares)
  typedef MeanSquareRegistrationFunction
        <FixedImageType, MovingImageType, DeformationFieldType> DefaultMetricType;
  typename DefaultMetricType::Pointer DefaultMetric
          = DefaultMetricType::New();
  this->m_PDEDeformableMetric[0] = DefaultMetric;
  this->m_PDEDeformableMetric[0]->SetNormalizeGradient( false );
  this->SetMaximizeMetric( false, 0 );

  this->m_PDEDeformableMetric[1] = NULL;

  // Set up the default image interpolator (linear)
  typedef LinearInterpolateImageFunction
          <MovingImageType, double> DefaultImageInterpolatorType;
  typename DefaultImageInterpolatorType::Pointer interpolator
          = DefaultImageInterpolatorType::New();
  this->m_ImageInterpolator = interpolator;

  this->m_Randomizer = RandomizerType::New();
  this->m_Randomizer->SetSeed(static_cast <ITK_UINT32> (0));
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::~FFDRegistrationFilter()
{
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void
FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::SetNumberOfLevels( unsigned int n )
{
  if ( this->m_NumberOfLevels != n )
    {
    this->m_NumberOfLevels = n;

    unsigned int tmp_MI;











    ArrayType tmp_MR;

    tmp_MI = 1;
    if ( this->m_MaximumNumberOfIterations.size() > 0 )
      {
      tmp_MI = this->m_MaximumNumberOfIterations[0];
      }

    this->m_MaximumNumberOfIterations.set_size( n );
    if ( this->m_NumberOfLevels != 0 )
      {
      this->m_MaximumNumberOfIterations.fill( tmp_MI );
      }
    this->Modified();
  }
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void
FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::GenerateData()
{
  if ( this->GetNumberOfInputs() < 2 )
    {
    itkExceptionMacro( << "Images are not specified." );
    }

  for ( this->m_CurrentLevel = 0; this->m_CurrentLevel < this->m_NumberOfLevels; this->m_CurrentLevel++ )
    {

    this->InitializeImages();

    std::cout << "Current Level = " << this->m_CurrentLevel+1 << " (" << this->m_NumberOfLevels << " total levels).  "
              << "Image Size = " << this->m_PyramidLevelImageSizes[this->m_CurrentLevel] << ".  "
              << "Grid Size = " << this->m_TotalDeformationFieldControlPoints->GetLargestPossibleRegion().GetSize() << ". " << std::endl;

    itkDebugMacro( "Current Level = " << this->m_CurrentLevel+1 << " (" << this->m_NumberOfLevels << " total levels).  "
              << "Image Size = " << this->m_PyramidLevelImageSizes[this->m_CurrentLevel] << ".  "
              << "Grid Size = " << this->m_TotalDeformationFieldControlPoints->GetLargestPossibleRegion().GetSize() << ". " );

    if ( this->m_MaximumNumberOfIterations[this->m_CurrentLevel] > 0 )
      {
      itkDebugMacro( "Begin iterative solve. " );
      this->IterativeSolve();
      itkDebugMacro( "Finish iterative solve. " );
      }

    itk::OStringStream buf;
    buf << "_" << this->m_CurrentLevel << "_" << this->GetSplineOrder();

    std::string filename = std::string( this->m_FilePrefix )
      + std::string( "ControlPointLattice" ) + buf.str()
      + std::string( ".nii.gz" );

    typedef VectorImageFileWriter<ControlPointLatticeType, RealImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( filename.c_str() );
    writer->SetInput( this->m_TotalDeformationFieldControlPoints );
    writer->Update();
    }

  typedef BSplineControlPointImageFilter
    <ControlPointLatticeType, DeformationFieldType> BSplineControlPointsFilterType;
  typename BSplineControlPointsFilterType::Pointer bspliner = BSplineControlPointsFilterType::New();
  typename BSplineControlPointsFilterType::ArrayType close;

  close.Fill( false );
  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetCloseDimension( close );
  bspliner->SetInput( this->m_TotalDeformationFieldControlPoints );
  bspliner->SetOrigin( this->m_CurrentFixedImage[0]->GetOrigin() );
  bspliner->SetSize( this->m_CurrentFixedImage[0]->GetLargestPossibleRegion().GetSize() );
  bspliner->SetSpacing( this->m_CurrentFixedImage[0]->GetSpacing() );
  bspliner->Update();

  typedef WarpImageFilter<MovingImageType,
                          WarpedImageType,
                          DeformationFieldType> WarperType;
  typename WarperType::Pointer warper = WarperType::New();

  warper->SetInput( this->GetInput( 1 ) );
  warper->SetDeformationField( bspliner->GetOutput() );
  warper->SetInterpolator( this->m_ImageInterpolator );
  warper->SetOutputSpacing( this->m_CurrentFixedImage[0]->GetSpacing() );
  warper->SetOutputOrigin( this->m_CurrentFixedImage[0]->GetOrigin() );
  warper->SetEdgePaddingValue( 0.0 );
  warper->Update();
  this->SetNthOutput( 0, warper->GetOutput() );
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void
FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::InitializeImages()
{
  if ( this->m_CurrentLevel == 0 )
    {
    for ( unsigned int m = 0; m < 2; m++ )
      {
      if ( !this->m_PDEDeformableMetric[m] )
        {
        continue;
        }
      typename FixedImagePyramidType::Pointer fixedPyramid = FixedImagePyramidType::New();
      fixedPyramid->SetInput( this->GetInput( 2*m ) );
      fixedPyramid->SetNumberOfLevels( this->m_NumberOfLevels );
      fixedPyramid->SetStartingShrinkFactors( this->m_FixedImageShrinkFactors.GetDataPointer() );

      fixedPyramid->Update();
      this->m_CurrentFixedImage[m] = fixedPyramid->GetOutput( 0 );

      typename MovingImagePyramidType::Pointer movingPyramid = MovingImagePyramidType::New();
      movingPyramid->SetInput( this->GetInput( 2*m+1 ) );
      movingPyramid->SetNumberOfLevels( this->m_NumberOfLevels );
      movingPyramid->SetStartingShrinkFactors( this->m_MovingImageShrinkFactors.GetDataPointer() );

      movingPyramid->Update();
      this->m_CurrentMovingImage[m] = movingPyramid->GetOutput( 0 );

      if ( this->m_Prolificacy )
        {
        OStringStream buf;
        buf << "_" << m << "_" << this->m_CurrentLevel;

        std::string filename = std::string( this->m_FilePrefix )
          + std::string( "FixedPyramid" ) + buf.str()
          + std::string( ".nii.gz" );

        typedef ImageFileWriter<FixedImageType> FixedWriterType;
        typename FixedWriterType::Pointer fixedWriter = FixedWriterType::New();
        fixedWriter->SetInput( fixedPyramid->GetOutput( 0 ) );
        fixedWriter->SetFileName( filename.c_str() );
        fixedWriter->Update();

        filename = std::string( this->m_FilePrefix )
          + std::string( "MovingPyramid" ) + buf.str()
          + std::string( ".nii.gz" );

        typedef ImageFileWriter<MovingImageType> MovingWriterType;
        typename MovingWriterType::Pointer movingWriter = MovingWriterType::New();
        movingWriter->SetInput( movingPyramid->GetOutput( 0 ) );
        movingWriter->SetFileName( filename.c_str() );
        movingWriter->Update();
        }

      if ( m == 0 )
        {
        this->m_FixedPyramidSchedule = fixedPyramid->GetSchedule();
        this->m_MovingPyramidSchedule = movingPyramid->GetSchedule();
        if ( this->m_NumberOfLevels == 0 )
          {
          this->m_PyramidLevelImageSizes.push_back(
             fixedPyramid->GetOutput( 0 )->GetLargestPossibleRegion().GetSize() );
          }
        else
          {
          for ( unsigned int i = 0; i < this->m_NumberOfLevels; i++ )
            {
            this->m_PyramidLevelImageSizes.push_back(
               fixedPyramid->GetOutput( i )->GetLargestPossibleRegion().GetSize() );
            }
          }
        }
      }

    if ( this->m_WeightImage )
      {
      typename WeightImagePyramidType::Pointer weightPyramid = WeightImagePyramidType::New();
      weightPyramid->SetInput( this->m_WeightImage );
      weightPyramid->SetNumberOfLevels( 1 );
      weightPyramid->SetStartingShrinkFactors( this->m_FixedImageShrinkFactors.GetDataPointer() );
      weightPyramid->Update();
      this->m_CurrentWeightImage = weightPyramid->GetOutput( 0 );

      if ( this->m_Prolificacy )
        {
        OStringStream buf;
        buf << "_" << this->m_CurrentLevel;

        std::string filename = std::string( this->m_FilePrefix )
          + std::string( "WeightPyramid" ) + buf.str()
          + std::string( ".nii.gz" );

        typedef ImageFileWriter<WeightImageType> WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetInput( this->m_CurrentWeightImage );
        writer->SetFileName( filename.c_str() );
        writer->Update();

        }
      }

    if ( this->m_InitializeWithLandmarks )
      {
      std::cout << "Initialize deformation field with the landmarks. " << std::endl;
      itkDebugMacro( "Initializing deformation field with the landmarks." );
      this->InitializeDeformationFieldWithLandmarks();
      }


    if ( !this->m_InitialControlPointLattice )
      {
      typename ControlPointLatticeType::RegionType::SizeType ncps;
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        unsigned int totalNumberOfBSplineLevels
          = static_cast<unsigned int>(
              floor( vcl_log( static_cast<RealType>( this->m_PyramidLevelImageSizes[this->m_NumberOfLevels-1][i]
                       /this->m_MeshResolution[i] ) )/vnl_math::ln2 ) );
        unsigned int nlevels = totalNumberOfBSplineLevels - this->m_NumberOfLevels + 2;
        ncps[i] = static_cast<unsigned int>( pow( 2, nlevels - 1 ) ) + this->m_SplineOrder;
        if ( ncps[i] < this->m_SplineOrder+1 )
          {
          ncps[i] = this->m_SplineOrder+1;
          }
        }

      VectorType V;
      V.Fill( 0 );

      this->m_TotalDeformationFieldControlPoints = ControlPointLatticeType::New();
      this->m_TotalDeformationFieldControlPoints->SetRegions( ncps );
      this->m_TotalDeformationFieldControlPoints->Allocate();
      this->m_TotalDeformationFieldControlPoints->FillBuffer( V );
      }
    else
      {
      this->m_TotalDeformationFieldControlPoints = this->m_InitialControlPointLattice;
      typename ControlPointLatticeType::SizeType size
        = this->m_TotalDeformationFieldControlPoints->GetLargestPossibleRegion().GetSize();
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        this->m_MeshResolution[i] = static_cast<unsigned int>(
          static_cast<RealType>( this->m_PyramidLevelImageSizes[0][i] ) /
          static_cast<RealType>( size[i] - this->m_SplineOrder ) + 0.5 );
        }
      }
    }
  else
    {
    ArrayType fixedImageShrinkFactors;
    ArrayType movingImageShrinkFactors;
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      fixedImageShrinkFactors[i] = this->m_FixedPyramidSchedule[this->m_CurrentLevel][i];
      movingImageShrinkFactors[i] = this->m_MovingPyramidSchedule[this->m_CurrentLevel][i];
      }

    for ( unsigned int m = 0; m < 2; m++ )
      {
      if ( !this->m_PDEDeformableMetric[m] )
        {
        continue;
        }
      typename FixedImagePyramidType::Pointer fixedPyramid = FixedImagePyramidType::New();
      fixedPyramid->SetInput( this->GetInput( 2*m ) );
      fixedPyramid->SetNumberOfLevels( 1 );
      fixedPyramid->SetStartingShrinkFactors( fixedImageShrinkFactors.GetDataPointer() );
      fixedPyramid->Update();
      this->m_CurrentFixedImage[m] = fixedPyramid->GetOutput( 0 );

      typename MovingImagePyramidType::Pointer movingPyramid = MovingImagePyramidType::New();
      movingPyramid->SetInput( this->GetInput( 2*m+1 ) );
      movingPyramid->SetNumberOfLevels( 1 );
      movingPyramid->SetStartingShrinkFactors( movingImageShrinkFactors.GetDataPointer() );
      movingPyramid->Update();
      this->m_CurrentMovingImage[m] = movingPyramid->GetOutput( 0 );

      if ( this->m_Prolificacy )
        {
        OStringStream buf;
        buf << "_" << m << "_" << this->m_CurrentLevel;

        std::string filename = std::string( this->m_FilePrefix )
          + std::string( "FixedPyramid" ) + buf.str()
          + std::string( ".nii.gz" );

        typedef ImageFileWriter<FixedImageType> FixedWriterType;
        typename FixedWriterType::Pointer fixedWriter = FixedWriterType::New();
        fixedWriter->SetInput( fixedPyramid->GetOutput( 0 ) );
        fixedWriter->SetFileName( filename.c_str() );
        fixedWriter->Update();

        filename = std::string( this->m_FilePrefix )
          + std::string( "MovingPyramid" ) + buf.str()
          + std::string( ".nii.gz" );

        typedef ImageFileWriter<MovingImageType> MovingWriterType;
        typename MovingWriterType::Pointer movingWriter = MovingWriterType::New();
        movingWriter->SetInput( movingPyramid->GetOutput( 0 ) );
        movingWriter->SetFileName( filename.c_str() );
        movingWriter->Update();
        }

      }
    if ( this->m_WeightImage )
      {
      typename WeightImagePyramidType::Pointer weightPyramid = WeightImagePyramidType::New();
      weightPyramid->SetInput( this->m_WeightImage );
      weightPyramid->SetNumberOfLevels( 1 );
      weightPyramid->SetStartingShrinkFactors( fixedImageShrinkFactors.GetDataPointer() );
      weightPyramid->Update();
      this->m_CurrentWeightImage = weightPyramid->GetOutput( 0 );
      }

    // Refine the control lattice

    if ( this->m_DoubleMeshResolutionAtEachLevel )
      {
      typedef BSplineControlPointImageFilter
        <ControlPointLatticeType, DeformationFieldType> BSplineControlPointsFilterType;
      typename BSplineControlPointsFilterType::Pointer bspliner = BSplineControlPointsFilterType::New();
      typename BSplineControlPointsFilterType::ArrayType close;

      close.Fill( false );
      bspliner->SetSplineOrder( this->m_SplineOrder );
      bspliner->SetCloseDimension( close );
      bspliner->SetInput( this->m_TotalDeformationFieldControlPoints );
      bspliner->SetOrigin( this->m_CurrentFixedImage[0]->GetOrigin() );
      bspliner->SetSpacing( this->m_CurrentFixedImage[0]->GetSpacing() );
      bspliner->SetSize( this->m_PyramidLevelImageSizes[this->m_CurrentLevel] );

      typename BSplineFilterType::ArrayType nlevels;
      nlevels.Fill( 2 );
      this->m_TotalDeformationFieldControlPoints = bspliner->RefineControlLattice( nlevels );
      }
    }
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void
FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
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

  this->m_GradientComputationTimes.push_back( 0.0 );
  this->m_NumberOfGradientPoints.push_back( 0 );

  RealType fp;
  if ( this->m_AffineNeighborhoodRadius > 0 )
    {
    this->EvaluateGradientFieldOverImageRegionForAffineTransform();
    fp = this->EvaluateMetricOverImageRegion( 0.0 );
    }
  else
    {
    fp = this->EvaluateGradientFieldOverImageRegion();
    }

  this->m_EnergyValues.push_back( fp );
  this->m_LevelNumbers.push_back( this->m_CurrentLevel );

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
//    std::cout << "  Current landmark error = " << this->CalculateLandmarkError() << std::endl;
    std::cout << "  Iteration " << its << ": Current Energy = " << fp << std::endl;
    itkDebugMacro( << "Iteration = " << its << ", Current Energy = " << fp );

    RealType gradientStep = 1.0;
    RealType fret;
    if ( this->m_LineSearchMaximumIterations > 0 )
      {
      this->LineMinimization( &gradientStep, &fret );
      }

    this->m_EnergyValues.push_back( fret );
    this->m_LevelNumbers.push_back( this->m_CurrentLevel );

    if ( 2.0*vnl_math_abs( fret - fp ) <= ftol*( vnl_math_abs( fret ) + vnl_math_abs( fp ) + eps ) &&
         this->m_LineSearchMaximumIterations > 0 )
      {
      std::cout << "  Exit condition (normal):  |fret - fp| = " << vnl_math_abs( fret - fp ) << std::endl;
      return;
      }

    fp = fret;

    typename RealImageType::Pointer lambda = RealImageType::New();
    lambda->SetOrigin( this->m_GradientFieldControlPoints->GetOrigin() );
    lambda->SetSpacing( this->m_GradientFieldControlPoints->GetSpacing() );
    lambda->SetRegions( this->m_GradientFieldControlPoints->GetLargestPossibleRegion().GetSize() );
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

    if ( this->m_CurrentLevel > 0 )
      {
      typename AdderType::Pointer adder = AdderType::New();
      adder->SetInput1( multiplier->GetOutput() );
      adder->SetInput2( this->m_CurrentDeformationFieldControlPoints );
      adder->Update();

      this->m_CurrentDeformationFieldControlPoints = adder->GetOutput();
      }

    if ( this->m_AffineNeighborhoodRadius > 0 )
      {
      this->EvaluateGradientFieldOverImageRegionForAffineTransform();
      }
    else
      {
      this->EvaluateGradientFieldOverImageRegion();
      }

    if ( this->m_Prolificacy )
      {
      OStringStream buf;
      buf << "_" << this->m_CurrentLevel << "_" << its-1;

      std::string filename = std::string( this->m_FilePrefix )
        + std::string( "DeformedGrid" ) + buf.str()
        + std::string( ".nii.gz" );

      typedef BSplineControlPointImageFilter
        <ControlPointLatticeType, DeformationFieldType> BSplineControlPointsFilterType;
      typename BSplineControlPointsFilterType::Pointer bspliner = BSplineControlPointsFilterType::New();
      typename BSplineControlPointsFilterType::ArrayType close;

      close.Fill( false );
      bspliner->SetSplineOrder( this->m_SplineOrder );
      bspliner->SetCloseDimension( close );
      bspliner->SetInput( this->m_TotalDeformationFieldControlPoints );
      bspliner->SetOrigin( this->GetInput( 0 )->GetOrigin() );
      bspliner->SetSize( this->GetInput( 0 )->GetLargestPossibleRegion().GetSize() );
      bspliner->SetSpacing( this->GetInput( 0 )->GetSpacing() );
      bspliner->Update();

      typedef GridImageSource<RealImageType> GridSourceType;
      typename GridSourceType::Pointer gridder = GridSourceType::New();
      gridder->SetSpacing( bspliner->GetOutput()->GetSpacing() );
      gridder->SetOrigin( bspliner->GetOutput()->GetOrigin() );
      gridder->SetSize( bspliner->GetOutput()->GetLargestPossibleRegion().GetSize() );

      typename GridSourceType::ArrayType gridSpacing;
      typename GridSourceType::ArrayType sigma;
      typename GridSourceType::BoolArrayType which;
      which.Fill( true );

      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        gridSpacing[i] = bspliner->GetOutput()->GetSpacing()[i] *
          ( bspliner->GetOutput()->GetLargestPossibleRegion().GetSize()[i] - 1 )
          / static_cast<RealType>(
           this->m_TotalDeformationFieldControlPoints->GetLargestPossibleRegion().GetSize()[i]
             - this->m_SplineOrder );
        sigma[i] = bspliner->GetOutput()->GetSpacing()[i] * 2.0;
        }

      gridder->SetGridSpacing( gridSpacing );
      gridder->SetSigma( sigma );
      gridder->SetWhichDimensions( which );
      gridder->Update();

      typedef WarpImageFilter<RealImageType, RealImageType, DeformationFieldType> WarperType;
      typename WarperType::Pointer warper = WarperType::New();
      warper->SetDeformationField( bspliner->GetOutput() );
      warper->SetInput( gridder->GetOutput() );
      warper->SetOutputOrigin( gridder->GetOutput()->GetOrigin() );
      warper->SetOutputSpacing( gridder->GetOutput()->GetSpacing() );
      warper->Update();

      typedef ImageFileWriter<RealImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( warper->GetOutput() );
      writer->SetFileName( filename.c_str() );
      writer->Update();

      filename = std::string( this->m_FilePrefix )
        + std::string( "DeformedMovingImage" ) + buf.str()
        + std::string( ".nii.gz" );

      warper->SetInput( this->GetInput( 1 ) );
      warper->SetOutputOrigin( this->GetInput( 1 )->GetOrigin() );
      warper->SetOutputSpacing( this->GetInput( 1 )->GetSpacing() );
      warper->Update();

      typedef ImageFileWriter<MovingImageType> MovingWriterType;
      typename MovingWriterType::Pointer movingWriter = MovingWriterType::New();
      movingWriter->SetInput( warper->GetOutput() );
      movingWriter->SetFileName( filename.c_str() );
      movingWriter->Update();
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
    RealType gamma = 0.0;
    if ( !this->m_EmploySteepestDescent )
      {
      gamma = vnl_math_max( static_cast<RealType>( 0.0 ), static_cast<RealType>( dgg / gg ) );
      }

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

template<class TMovingImage, class TFixedImage, class TWarpedImage>
typename FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>::RealType
FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::EvaluateGradientFieldOverImageRegion()
{
  const RealType SmallGradientValue = 1.0e-20;

  itkDebugMacro( "Evaluating gradient fields." );

  typedef BSplineControlPointImageFilter
    <ControlPointLatticeType, DeformationFieldType> BSplineControlPointsFilterType;
  typename BSplineControlPointsFilterType::Pointer bspliner = BSplineControlPointsFilterType::New();
  typename BSplineControlPointsFilterType::ArrayType close;

  close.Fill( false );
  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetCloseDimension( close );
  bspliner->SetInput( this->m_TotalDeformationFieldControlPoints );
  bspliner->SetOrigin( this->m_CurrentFixedImage[0]->GetOrigin() );
  bspliner->SetSize( this->m_CurrentFixedImage[0]->GetLargestPossibleRegion().GetSize() );
  bspliner->SetSpacing( this->m_CurrentFixedImage[0]->GetSpacing() );
  bspliner->Update();

  typename DeformationFieldType::Pointer deformationField = bspliner->GetOutput();  

  typename PointSetType::Pointer fieldPoints = PointSetType::New();
  fieldPoints->Initialize();

  RealType metricEnergy[2];
  RealType metricCount[2];

  for ( unsigned int m = 0; m < 2; m++ )
    {
    metricEnergy[m] = 0.0;
    metricCount[m] = 0.0;
    if ( !this->m_PDEDeformableMetric[m] )
      {
      continue;
      }

    typedef WarpImageFilter<MovingImageType,
                            MovingImageType,
                            DeformationFieldType> WarperType;
    typename WarperType::Pointer warper = WarperType::New();
    warper->SetInput( this->m_CurrentMovingImage[m] );
    warper->SetDeformationField( deformationField );
    warper->SetInterpolator( this->m_ImageInterpolator );
    warper->SetOutputSpacing( this->m_CurrentFixedImage[m]->GetSpacing() );
    warper->SetOutputOrigin( this->m_CurrentFixedImage[m]->GetOrigin() );
    warper->SetEdgePaddingValue( 0.0 );
    warper->Update();

    this->m_PDEDeformableMetric[m]->SetRadius( this->m_MetricRadius[m] );
    this->m_PDEDeformableMetric[m]->SetMovingImage( warper->GetOutput() );
    this->m_PDEDeformableMetric[m]->SetFixedImage( this->m_CurrentFixedImage[m] );
    this->m_PDEDeformableMetric[m]->SetDeformationField( NULL );
//    this->m_PDEDeformableMetric[m]->SetMaskImage( this->m_CurrentWeightImage );
    this->m_PDEDeformableMetric[m]->InitializeIteration();

    typedef typename NeighborhoodAlgorithm
      ::ImageBoundaryFacesCalculator<DeformationFieldType> FaceCalculatorType;
    FaceCalculatorType faceCalculator;

    typename FaceCalculatorType::FaceListType faceList = faceCalculator( deformationField,
      deformationField->GetLargestPossibleRegion(), this->m_MetricRadius[m] );
    typename FaceCalculatorType::FaceListType::iterator fit;

    for ( fit = faceList.begin(); fit != faceList.end(); ++fit )
      {
      NeighborhoodIteratorType It( this->m_MetricRadius[m], deformationField, *fit );
      for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
        {
        if ( this->m_WeightImage &&
             ( this->m_CurrentWeightImage->GetPixel( It.GetIndex() ) < 0 ||
               ( !this->m_PDEDeformableMetric[1] &&
                 this->m_CurrentWeightImage->GetPixel( It.GetIndex() ) <= 0 ) ) )
          {
          continue;
          }

        this->m_PDEDeformableMetric[m]->SetEnergy( 0.0 );
        VectorType grad = this->m_PDEDeformableMetric[m]->ComputeUpdate( It, NULL );
        RealType metric = this->m_PDEDeformableMetric[m]->GetEnergy();
        if ( !vnl_math_isnan( metric ) )
          {
//          if ( this->m_MaximizeMetric[m] )
//            {
//            metric = -metric;
//            }
          if ( this->m_WeightImage && this->m_PDEDeformableMetric[1] )
            {
            if ( m == 0 )
              {
              metric *= this->m_CurrentWeightImage->GetPixel( It.GetIndex() );
              }
            else
              {
              metric *= ( 1.0 - this->m_CurrentWeightImage->GetPixel( It.GetIndex() ) );
              }
            }
          metricEnergy[m] += metric;
          metricCount[m] += 1.0;
          }
          
        for ( unsigned int i = 0; i < ImageDimension; i++ )
          {
          if ( !this->m_Directionality[i] )
            {
            grad[i] = 0.0;
            }
          }
          
        if ( grad.GetSquaredNorm() >= SmallGradientValue &&
             grad[0] < NumericTraits<RealType>::max()-1 )
          {
          typename PointSetType::PointType point;
          this->m_CurrentFixedImage[m]->TransformIndexToPhysicalPoint( It.GetIndex(), point );

          unsigned int count = fieldPoints->GetNumberOfPoints();
          fieldPoints->SetPoint( count, point );
//          if ( this->m_MaximizeMetric[m] )
//            {
            fieldPoints->SetPointData( count, -grad );
//            }
//          else
//            {
//            fieldPoints->SetPointData( count, grad );
//            }
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

  typename PointSetType::PointDataContainer::Iterator 
    ItP = fieldPoints->GetPointData()->Begin();

  RealType sumSquaredNorm = 0.0;
  RealType N = 0.0;
  while ( ItP != fieldPoints->GetPointData()->End() )
    {
    sumSquaredNorm += ( ItP.Value() ).GetSquaredNorm();
    N += 1.0;
    ++ItP;
    }

  RealType sigma;
  VectorType V;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    V[i] = this->m_CurrentFixedImage[0]->GetSpacing()[i];
    }
  if ( ImageDimension == 2 )
    {
    sigma = V.GetNorm();
    }
  else if ( ImageDimension == 3 )
    {
    sigma = V.GetNorm()/sqrt( 2.0 );
    }
  RealType gradientScalingFactor = sigma*vcl_sqrt
    ( static_cast<RealType>( ImageDimension )*N / sumSquaredNorm );

  ItP = fieldPoints->GetPointData()->Begin();
  while ( ItP != fieldPoints->GetPointData()->End() )
    {
    fieldPoints->GetPointData()->InsertElement( ItP.Index(), gradientScalingFactor*ItP.Value() );
    ++ItP;
    }

  /**
   * Calculate contribution from Jacobian constraint
   */
  RealType jacobianEnergy = 0.0;
  RealType jacobianCount = 0.0;

  if ( this->m_EnforceDiffeomorphism )
    {
    typedef VectorFieldGradientImageFunction<DeformationFieldType> FunctionType;
    typename FunctionType::Pointer function = FunctionType::New();

    if ( this->m_CurrentLevel == 0 )
      {
      function->SetInputImage( deformationField );
      }
    else
      {
      typedef BSplineControlPointImageFilter
        <ControlPointLatticeType, DeformationFieldType> BSplineControlPointsFilterType;
      typename BSplineControlPointsFilterType::Pointer bspliner = BSplineControlPointsFilterType::New();
      typename BSplineControlPointsFilterType::ArrayType close;

      close.Fill( false );
      bspliner->SetSplineOrder( this->m_SplineOrder );
      bspliner->SetCloseDimension( close );
      bspliner->SetInput( this->m_CurrentDeformationFieldControlPoints );
      bspliner->SetOrigin( this->m_CurrentFixedImage[0]->GetOrigin() );
      bspliner->SetSize( this->m_CurrentFixedImage[0]->GetLargestPossibleRegion().GetSize() );
      bspliner->SetSpacing( this->m_CurrentFixedImage[0]->GetSpacing() );
      bspliner->Update();

      function->SetInputImage( bspliner->GetOutput() );
      }

    typename RealImageType::Pointer jacobianImage = RealImageType::New();
    jacobianImage->SetRegions( this->m_CurrentFixedImage[0]->GetLargestPossibleRegion() );
    jacobianImage->SetOrigin( this->m_CurrentFixedImage[0]->GetOrigin() );
    jacobianImage->SetSpacing( this->m_CurrentFixedImage[0]->GetSpacing() );
    jacobianImage->Allocate();
    jacobianImage->FillBuffer( 1.0 );

    ImageRegionIterator<RealImageType> ItJ
      ( jacobianImage, jacobianImage->GetLargestPossibleRegion() );

    for ( ItJ.GoToBegin(); !ItJ.IsAtEnd(); ++ItJ )
      {
      ItJ.Set( function->EvaluateJacobianDeterminantAtIndex( ItJ.GetIndex() ) );
      }

    typedef DerivativeImageFilter<RealImageType, RealImageType> DerivativeFilterType;
    typename DerivativeFilterType::Pointer derivativeFilter[ImageDimension];

    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      derivativeFilter[i] = DerivativeFilterType::New();
      derivativeFilter[i]->SetInput( jacobianImage );
      derivativeFilter[i]->SetUseImageSpacingOn();
      derivativeFilter[i]->SetOrder( 1 );
      derivativeFilter[i]->SetDirection( i );
      derivativeFilter[i]->Update();
      }

    typedef LinearInterpolateImageFunction<RealImageType> InterpolatorType;

    typename InterpolatorType::Pointer jacobianInterpolator = InterpolatorType::New();
    jacobianInterpolator->SetInputImage( jacobianImage );

    typename PointSetType::PointDataContainer::Iterator ItP;
    typename PointSetType::PointsContainer::ConstIterator ItQ;
    ItP = fieldPoints->GetPointData()->Begin();
    ItQ = fieldPoints->GetPoints()->Begin();

    while ( ItP != fieldPoints->GetPointData()->End() )
      {
      typename InterpolatorType::PointType point;
      for ( unsigned int d = 0; d < ImageDimension; d++ )
        {
        point[d] = ItQ.Value()[d];
        }
      RealType detJ = jacobianInterpolator->Evaluate( point );
      if ( !vnl_math_isnan( detJ ) )
        {
        if ( detJ < this->m_MinimumJacobian || detJ > this->m_MaximumJacobian )
          {
          jacobianEnergy = NumericTraits<RealType>::max();
          jacobianCount = 1.0;
          break;
          }
        else
          {
/*
          RealType g = vcl_exp( 100*( detJ - this->m_MinimumJacobian ) );
          RealType h = vcl_exp( 100*( this->m_MaximumJacobian - detJ ) );
          jacobianEnergy += 1e10 * this->m_InteriorPenaltyParameter *
            ( 1.0 / ( 1.0 + g ) + 1.0 / ( 1.0 + h ) );
*/
          jacobianEnergy -= this->m_InteriorPenaltyParameter *
            ( vcl_log( detJ - this->m_MinimumJacobian )
            + vcl_log( this->m_MaximumJacobian - detJ ) );
          jacobianCount += 1.0;
          }
        }
      if ( detJ > this->m_MinimumJacobian && detJ < this->m_MaximumJacobian )
        {
        VectorType V;
        RealType g = this->m_MinimumJacobian - detJ;
        RealType h = detJ - this->m_MaximumJacobian;

        for ( unsigned int d = 0; d < ImageDimension; d++ )
          {
          typename InterpolatorType::Pointer derivativeInterpolator = InterpolatorType::New();
          derivativeInterpolator->SetInputImage( derivativeFilter[d]->GetOutput() );
          RealType derJ = derivativeInterpolator->Evaluate( point );
          V[d] = derJ / g + derJ / h;
          }
        fieldPoints->GetPointData()->InsertElement( ItP.Index(), ItP.Value() - this->m_InteriorPenaltyParameter * V );
        }
      ++ItP;
      ++ItQ;
      }
    }

  itkDebugMacro( "Fitting B-spline field to gradient values." );


  RealType landmarkCount = 0.0;
  RealType landmarkEnergy = 0.0;

/*
  if ( this->m_FixedLandmarkContainer->Size() > 0 )
    {
    itkDebugMacro( "Extracting landmark information." );

    typedef BSplineControlPointImageFilter
      <ControlPointLatticeType, DeformationFieldType> BSplineControlPointsFilterType;
    typename BSplineControlPointsFilterType::Pointer bspliner = BSplineControlPointsFilterType::New();
    typename BSplineControlPointsFilterType::ArrayType close;

    close.Fill( false );
    bspliner->SetSplineOrder( this->m_SplineOrder );
    bspliner->SetCloseDimension( close );
    bspliner->SetInput( this->m_TotalDeformationFieldControlPoints );
    bspliner->SetOrigin( this->m_CurrentMovingImage[0]->GetOrigin() );
    bspliner->SetSize(
      this->m_CurrentMovingImage[0]->GetLargestPossibleRegion().GetSize() );
    bspliner->SetSpacing( this->m_CurrentMovingImage[0]->GetSpacing() );


    typename LandmarkContainerType::Pointer deformedMovingLandmarkContainer
      = LandmarkContainerType::New();
    deformedMovingLandmarkContainer->Initialize();

    for ( unsigned int i = 0; i < this->m_MovingLandmarkContainer->Size(); i++ )
      {
      LandmarkType moving = this->m_MovingLandmarkContainer->GetElement( i );
      typename BSplineControlPointsFilterType::PointDataType point;
      typename BSplineControlPointsFilterType::PointType params;

      for ( unsigned int d = 0; d < ImageDimension; d++ )
        {
        point[d] = moving[d];
        params[d] = ( point[d] - bspliner->GetOrigin()[d] ) /
          ( bspliner->GetSpacing()[d] * ( bspliner->GetSize()[d] - 1 ) );
        }

      bspliner->EvaluateParametersAtPoint( point, params );
      LandmarkType movingWarped;
      for ( unsigned int d = 0; d < ImageDimension; d++ )
        {
        movingWarped[d] = params[d] *
          ( bspliner->GetSpacing()[d] * ( bspliner->GetSize()[d] - 1 ) )
          + bspliner->GetOrigin()[d];
        }
      movingWarped[ImageDimension] = moving[ImageDimension];
      deformedMovingLandmarkContainer->InsertElement( i, movingWarped );
      }

    for ( unsigned int i = 0; i < this->m_FixedLandmarkContainer->Size(); i++ )
      {
      LandmarkType fixed = this->m_FixedLandmarkContainer->GetElement( i );

      RealType min_distance = NumericTraits<RealType>::max();
      int min_idx = -1;

      for ( unsigned int j = 0; j < deformedMovingLandmarkContainer->Size(); j++ )
        {
        LandmarkType moving = deformedMovingLandmarkContainer->GetElement( j );
        if ( static_cast<unsigned int>( fixed[ImageDimension] ) ==
             static_cast<unsigned int>( moving[ImageDimension] ) )
          {
          if ( ( moving - fixed ).GetSquaredNorm() < min_distance )
            {
            min_idx = j;
            min_distance = ( moving - fixed ).GetSquaredNorm();
            }
          }
        }
      if ( min_idx != -1 )
        {
        landmarkEnergy += ( 0.5 * this->m_LandmarkWeighting * min_distance );
        landmarkCount += 1.0;

        LandmarkType moving = this->m_MovingLandmarkContainer->GetElement( min_idx );
        typename PointSetType::PointType point;
        typename DeformationFieldType::PixelType data;
        for ( unsigned int j = 0; j < ImageDimension; j++ )
          {
          point[j] = fixed[j];
          data[j] = moving[j] - fixed[j];
          }
        fieldPoints->SetPoint( index, point );
        fieldPoints->SetPointData( index, data );
        weights->InsertElement( index, this->m_LandmarkWeighting );
        index++;
        }
      }
    }
*/
  // Define the rest of the parameters for the B-spline filter

  typename DeformationFieldType::PointType origin;
  typename DeformationFieldType::SpacingType spacing;
  typename DeformationFieldType::SizeType size;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    origin[i] = this->m_CurrentFixedImage[0]->GetOrigin()[i];
    spacing[i] = this->m_CurrentFixedImage[0]->GetSpacing()[i];
    size[i] = this->m_CurrentFixedImage[0]->GetLargestPossibleRegion().GetSize()[i];
    }

  typename BSplineFilterType::ArrayType nlevels;
  typename BSplineFilterType::ArrayType ncps;

  nlevels.Fill( 1 );
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    ncps[i] = this->m_TotalDeformationFieldControlPoints->GetLargestPossibleRegion().GetSize()[i];
    }

  TimeProbe timer;
  timer.Start();

  switch ( this->m_WhichGradient )
    {
    case 0: default:
      {

      typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
      bspliner->SetInput( fieldPoints );
      bspliner->SetOrigin( origin );
      bspliner->SetSpacing( spacing );
      bspliner->SetSize( size );
      bspliner->SetNumberOfLevels( nlevels );
      bspliner->SetSplineOrder( this->m_SplineOrder );
      bspliner->SetNumberOfControlPoints( ncps );
      bspliner->SetGenerateOutputImage( false );
      bspliner->Update();

      this->m_GradientFieldControlPoints = bspliner->GetPhiLattice();
      break;
      }
    case 1:
      {
      typedef BSplineControlPointImageFilter
        <ControlPointLatticeType, DeformationFieldType> BSplineControlPointsFilterType;
      typename BSplineControlPointsFilterType::Pointer bspliner = BSplineControlPointsFilterType::New();
      typename BSplineControlPointsFilterType::ArrayType close;

      close.Fill( false );
      bspliner->SetSplineOrder( this->m_SplineOrder );
      bspliner->SetInput( this->m_TotalDeformationFieldControlPoints );
      bspliner->SetOrigin( origin );
      bspliner->SetSize( size );
      bspliner->SetSpacing( spacing );

      this->m_GradientFieldControlPoints
        = bspliner->CalculateLatticeGradientFromPoints( fieldPoints );

      ImageRegionIterator<ControlPointLatticeType> ItCP( this->m_GradientFieldControlPoints,
        this->m_GradientFieldControlPoints->GetLargestPossibleRegion() );
  
      RealType avg = 0.0;
      RealType N = 0.0;
      for ( ItCP.GoToBegin(); !ItCP.IsAtEnd(); ++ItCP )
        {
        avg += ItCP.Get().GetNorm();
        N += 1.0;
        }
      avg /= N;
  
      RealType scaleFactor = 0.0;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        scaleFactor += vnl_math_sqr( spacing[d] );
        }
      scaleFactor = vcl_sqrt( scaleFactor );
  
      for ( ItCP.GoToBegin(); !ItCP.IsAtEnd(); ++ItCP )
        {
        ItCP.Set( scaleFactor * ItCP.Get() / avg  );
        }
        
      break;
      }
    }

  timer.Stop();
  this->m_GradientComputationTimes.push_back( timer.GetMeanTime() );
  this->m_NumberOfGradientPoints.push_back( fieldPoints->GetNumberOfPoints() );


  RealType energy = 0.0;
  if ( metricCount[0] > 0.0 )
    {
    energy += metricEnergy[0] / metricCount[0];
//    std::cout << " metric = " << metricEnergy / metricCount << std::flush;
    }
  if ( metricCount[1] > 0.0 )
    {
    energy += metricEnergy[1] / metricCount[1];
//    std::cout << " metric = " << metricEnergy / metricCount << std::flush;
    }
  if ( landmarkCount > 0.0 )
    {
    energy += landmarkEnergy / landmarkCount;
//    std::cout << " landmark = " << landmarkEnergy / landmarkCount << std::flush;
    }
  if ( jacobianCount > 0.0 )
    {
    energy += jacobianEnergy / jacobianCount;
//    std::cout << " jacobian = " << jacobianEnergy / jacobianCount << std::flush;
    }
//  std::cout << std::endl;


  // Normalize the gradient

  return energy;
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void
FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::EvaluateGradientFieldOverImageRegionForAffineTransform()
{
  const RealType SmallGradientValue = 1.0e-20;

  itkDebugMacro( "Evaluating gradient field." );

  typedef BSplineControlPointImageFilter
    <ControlPointLatticeType, DeformationFieldType> BSplineControlPointsFilterType;
  typename BSplineControlPointsFilterType::Pointer bspliner = BSplineControlPointsFilterType::New();
  typename BSplineControlPointsFilterType::ArrayType close;

  close.Fill( false );
  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetCloseDimension( close );
  bspliner->SetInput( this->m_TotalDeformationFieldControlPoints );
  bspliner->SetOrigin( this->m_CurrentFixedImage[0]->GetOrigin() );
  bspliner->SetSize( this->m_CurrentFixedImage[0]->GetLargestPossibleRegion().GetSize() );
  bspliner->SetSpacing( this->m_CurrentFixedImage[0]->GetSpacing() );
  bspliner->Update();

  typename DeformationFieldType::Pointer deformationField = bspliner->GetOutput();

  typename PointSetType::Pointer fieldPoints = PointSetType::New();
  fieldPoints->Initialize();

  for ( unsigned int m = 0; m < 2; m++ )
    {
    if ( !this->m_PDEDeformableMetric[m] )
      {
      continue;
      }

    typedef WarpImageFilter<MovingImageType,
                            MovingImageType,
                            DeformationFieldType> WarperType;
    typename WarperType::Pointer warper = WarperType::New();
    warper->SetInput( this->m_CurrentMovingImage[m] );
    warper->SetDeformationField( deformationField );
    warper->SetInterpolator( this->m_ImageInterpolator );
    warper->SetOutputSpacing( this->m_CurrentFixedImage[m]->GetSpacing() );
    warper->SetOutputOrigin( this->m_CurrentFixedImage[m]->GetOrigin() );
    warper->SetEdgePaddingValue( 0 );
    warper->Update();

    this->m_PDEDeformableMetric[m]->SetRadius( this->m_MetricRadius[m] );
    this->m_PDEDeformableMetric[m]->SetMovingImage( warper->GetOutput() );
    this->m_PDEDeformableMetric[m]->SetFixedImage( this->m_CurrentFixedImage[m] );
    this->m_PDEDeformableMetric[m]->SetDeformationField( NULL );
//    this->m_PDEDeformableMetric[m]->SetMaskImage( this->m_CurrentWeightImage );
    this->m_PDEDeformableMetric[m]->InitializeIteration();

    NeighborhoodIteratorType It( this->m_MetricRadius[m], deformationField,
      deformationField->GetLargestPossibleRegion() );
    It.GoToBegin();

    typename ShapedNeighborhoodIteratorType::RadiusType affineRadius;
    affineRadius.Fill( this->m_AffineNeighborhoodRadius );

    ShapedNeighborhoodIteratorType ItS( affineRadius, deformationField,
      deformationField->GetLargestPossibleRegion() );
    ItS.GoToBegin();
    ItS.ClearActiveList();
    unsigned int numberOfNeighbors = 1;
    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      numberOfNeighbors *= ( 2*this->m_AffineNeighborhoodRadius + 1 );
      }

    Array<unsigned char> active;
    active.SetSize( numberOfNeighbors );
    active.Fill( false );

    unsigned int n = 0;
    while ( n < this->m_AffineNeighborhoodSamplingDensity * numberOfNeighbors )
      {
      unsigned int whichNeighbor
        = this->m_Randomizer->GetIntegerVariate( numberOfNeighbors-1 );
      if ( active[whichNeighbor] == false )
        {
        ItS.ActivateOffset( ItS.GetOffset( whichNeighbor ) );
        active[whichNeighbor] = true;
        n++;
        }
      }

    typedef ImageRandomIteratorWithIndex<DeformationFieldType> RandomIteratorType;
    RandomIteratorType ItR( deformationField, deformationField->GetLargestPossibleRegion() );
    ItR.SetNumberOfSamples( this->m_NumberOfAffineNeighborhoodSamplesPerIteration );
    ItR.ReinitializeSeed();
    ItR.GoToBegin();

    unsigned int count = 0;
    while ( !ItR.IsAtEnd() )
      {
      if ( this->m_WeightImage &&
           ( this->m_CurrentWeightImage->GetPixel( ItR.GetIndex() ) < 0 ||
             ( !this->m_PDEDeformableMetric[1] &&
               this->m_CurrentWeightImage->GetPixel( ItR.GetIndex() ) <= 0 ) ) )
        {
        ++ItR;
        ItR.SetNumberOfSamples( ItR.GetNumberOfSamples() + 1 );
        continue;
        }
      ItS.SetLocation( ItR.GetIndex() );

      typename ShapedNeighborhoodIteratorType::IndexListType indexList
        = ItS.GetActiveIndexList();
      typename ShapedNeighborhoodIteratorType::IndexListType::const_iterator iter;
      for ( iter = indexList.begin(); iter != indexList.end(); iter++ )
        {
        typename DeformationFieldType::IndexType index = ItS.GetIndex( *iter );
        if ( !deformationField->GetLargestPossibleRegion().IsInside( index ) )
          {
          continue;
          }
        It.SetLocation( index );
        this->m_PDEDeformableMetric[m]->SetEnergy( 0.0 );
        VectorType grad = this->m_PDEDeformableMetric[m]->ComputeUpdate( It, NULL );

        if ( vnl_math_isnan( this->m_PDEDeformableMetric[m]->GetEnergy() ) )
          {
          continue;
          }
        for ( unsigned int i = 0; i < ImageDimension; i++ )
          {
          if ( !this->m_Directionality[i] )
            {
            grad[i] = 0.0;
            }
          }
        if ( grad.GetSquaredNorm() >= SmallGradientValue &&
             grad[0] < NumericTraits<RealType>::max()-1 )
          {
          typename PointSetType::PointType point;
          this->m_CurrentFixedImage[m]->TransformIndexToPhysicalPoint( index, point );

          fieldPoints->SetPoint( count, point );
//          if ( this->m_MaximizeMetric[m] )
//            {
            fieldPoints->SetPointData( count, -grad );
//            }
//          else
//            {
//            fieldPoints->SetPointData( count, grad );
//            }
          count++;
          }
        }
      if ( count == 0 )
        {
        ItR.SetNumberOfSamples( ItR.GetNumberOfSamples() + 1 );
        }
      ++ItR;
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
  RealType N = 0.0;
  typename PointSetType::PointDataContainer::Iterator ItP;
  ItP = fieldPoints->GetPointData()->Begin();
  while ( ItP != fieldPoints->GetPointData()->End() )
    {
    sumSquaredNorm += ( ItP.Value() ).GetSquaredNorm();
    N += 1.0;
    ++ItP;
    }

  RealType sigma;
  VectorType V;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    V[i] = this->m_CurrentFixedImage[0]->GetSpacing()[i];
    }
  if ( ImageDimension == 2 )
    {
    sigma = V.GetNorm();
    }
  else if ( ImageDimension == 3 )
    {
    sigma = V.GetNorm()/sqrt( 2.0 );
    }
  RealType gradientScalingFactor = sigma*sqrt
    ( static_cast<RealType>( ImageDimension )*N / sumSquaredNorm );


  ItP = fieldPoints->GetPointData()->Begin();
  while ( ItP != fieldPoints->GetPointData()->End() )
    {
    fieldPoints->GetPointData()->InsertElement( ItP.Index(), gradientScalingFactor*ItP.Value() );
    ++ItP;
    }

  itkDebugMacro( "Fitting B-spline field to gradient values." );

  // Define the rest of the parameters for the B-spline filter

  typename DeformationFieldType::PointType origin;
  typename DeformationFieldType::SpacingType spacing;
  typename DeformationFieldType::SizeType size;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    origin[i] = this->m_CurrentFixedImage[0]->GetOrigin()[i];
    spacing[i] = this->m_CurrentFixedImage[0]->GetSpacing()[i];
    size[i] = this->m_CurrentFixedImage[0]->GetLargestPossibleRegion().GetSize()[i];
    }

  typename BSplineFilterType::ArrayType nlevels;
  typename BSplineFilterType::ArrayType ncps;

  nlevels.Fill( 1 );
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    ncps[i] = this->m_TotalDeformationFieldControlPoints->GetLargestPossibleRegion().GetSize()[i];
    }

  TimeProbe timer;
  timer.Start();

  switch ( this->m_WhichGradient )
    {
    case 0: default:
      {

      typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
      bspliner->SetInput( fieldPoints );
      bspliner->SetOrigin( origin );
      bspliner->SetSpacing( spacing );
      bspliner->SetSize( size );
      bspliner->SetNumberOfLevels( nlevels );
      bspliner->SetSplineOrder( this->m_SplineOrder );
      bspliner->SetNumberOfControlPoints( ncps );
      bspliner->SetGenerateOutputImage( false );
      bspliner->Update();

      this->m_GradientFieldControlPoints = bspliner->GetPhiLattice();
      break;
      }
    case 1:
      {
      typedef BSplineControlPointImageFilter
        <ControlPointLatticeType, DeformationFieldType> BSplineControlPointsFilterType;
      typename BSplineControlPointsFilterType::Pointer bspliner = BSplineControlPointsFilterType::New();
      typename BSplineControlPointsFilterType::ArrayType close;

      close.Fill( false );
      bspliner->SetSplineOrder( this->m_SplineOrder );
      bspliner->SetInput( this->m_TotalDeformationFieldControlPoints );
      bspliner->SetOrigin( origin );
      bspliner->SetSize( size );
      bspliner->SetSpacing( spacing );

      this->m_GradientFieldControlPoints
        = bspliner->CalculateLatticeGradientFromPoints( fieldPoints );

      break;
      }
    }
  timer.Stop();

  this->m_GradientComputationTimes.push_back( timer.GetMeanTime() );
  this->m_NumberOfGradientPoints.push_back( fieldPoints->GetNumberOfPoints() );

}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
typename FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>::RealType
FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::EvaluateMetricOverImageRegion( RealType t = 0 )
{

  typename RealImageType::Pointer lambda = RealImageType::New();
  lambda->SetOrigin( this->m_GradientFieldControlPoints->GetOrigin() );
  lambda->SetSpacing( this->m_GradientFieldControlPoints->GetSpacing() );
  lambda->SetRegions( this->m_GradientFieldControlPoints->GetLargestPossibleRegion().GetSize() );
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

  typedef BSplineControlPointImageFilter
    <ControlPointLatticeType, DeformationFieldType> BSplineControlPointsFilterType;
  typename BSplineControlPointsFilterType::Pointer bspliner = BSplineControlPointsFilterType::New();
  typename BSplineControlPointsFilterType::ArrayType close;

  close.Fill( false );
  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetCloseDimension( close );
  bspliner->SetInput( adder->GetOutput() );
  bspliner->SetOrigin( this->m_CurrentFixedImage[0]->GetOrigin() );
  bspliner->SetSize(
    this->m_CurrentFixedImage[0]->GetLargestPossibleRegion().GetSize() );
  bspliner->SetSpacing( this->m_CurrentFixedImage[0]->GetSpacing() );
  bspliner->Update();

  typename DeformationFieldType::Pointer deformationField = bspliner->GetOutput();

/*
  ImageRegionIterator<DeformationFieldType> It( deformationField,
    deformationField->GetLargestPossibleRegion() );
  this->m_MeanDeformationFieldVectorLength = 0.0;
  RealType N = 0.0;
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    this->m_MeanDeformationFieldVectorLength += It.Get().GetSquaredNorm();
    N += 1.0;
    }
  this->m_MeanDeformationFieldVectorLength /= N;
  this->m_MeanDeformationFieldVectorLength = vcl_sqrt( this->m_MeanDeformationFieldVectorLength );
*/
  VectorType V;
  V.Fill( 0 );

  RealType metricEnergy[2];
  RealType metricCount[2];

  for ( unsigned int m = 0; m < 2; m++ )
    {
    metricEnergy[m] = 0.0;
    metricCount[m] = 0.0;
    if ( !this->m_PDEDeformableMetric[m] )
      {
      continue;
      }

    typedef WarpImageFilter<MovingImageType,
                            MovingImageType,
                            DeformationFieldType> WarperType;
    typename WarperType::Pointer warper = WarperType::New();
    warper->SetInput( this->m_CurrentMovingImage[m] );
    warper->SetDeformationField( deformationField );
    warper->SetInterpolator( this->m_ImageInterpolator );
    warper->SetOutputSpacing( this->m_CurrentMovingImage[0]->GetSpacing() );
    warper->SetOutputOrigin( this->m_CurrentMovingImage[0]->GetOrigin() );
    warper->SetEdgePaddingValue( 0.0 );
    warper->Update();

    this->m_PDEDeformableMetric[m]->SetRadius( this->m_MetricRadius[m] );
    this->m_PDEDeformableMetric[m]->SetMovingImage( warper->GetOutput() );
    this->m_PDEDeformableMetric[m]->SetFixedImage( this->m_CurrentFixedImage[m] );
    this->m_PDEDeformableMetric[m]->SetDeformationField( NULL );
    this->m_PDEDeformableMetric[m]->InitializeIteration();

    typedef typename NeighborhoodAlgorithm
      ::ImageBoundaryFacesCalculator<DeformationFieldType> FaceCalculatorType;
    FaceCalculatorType faceCalculator;

    typename FaceCalculatorType::FaceListType faceList = faceCalculator( deformationField,
      deformationField->GetLargestPossibleRegion(), this->m_MetricRadius[m] );
    typename FaceCalculatorType::FaceListType::iterator fit;

    for ( fit = faceList.begin(); fit != faceList.end(); ++fit )
      {
      NeighborhoodIteratorType It( this->m_MetricRadius[m], deformationField, *fit );
      for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
        {
        if ( this->m_WeightImage &&
             ( this->m_CurrentWeightImage->GetPixel( It.GetIndex() ) < 0 ||
               ( !this->m_PDEDeformableMetric[1] &&
                 this->m_CurrentWeightImage->GetPixel( It.GetIndex() ) <= 0 ) ) )
          {
          continue;
          }
        this->m_PDEDeformableMetric[m]->SetEnergy( 0.0 );
        this->m_PDEDeformableMetric[m]->ComputeUpdate( It, NULL );
        RealType metric = this->m_PDEDeformableMetric[m]->GetEnergy();
        if ( !vnl_math_isnan( metric ) )
          {
//          if ( this->m_MaximizeMetric[m] )
//            {
//            metric = -metric;
//            }
          if ( this->m_WeightImage && this->m_PDEDeformableMetric[1] )
            {
            if ( m == 0 )
              {
              metric *= this->m_CurrentWeightImage->GetPixel( It.GetIndex() );
              }
            else
              {
              metric *= ( 1.0 - this->m_CurrentWeightImage->GetPixel( It.GetIndex() ) );
              }
            }
          metricEnergy[m] += metric;
          metricCount[m] += 1.0;
          }
        }
      }
    }


  RealType landmarkCount = 0.0;
  RealType landmarkEnergy = 0.0;

/*
  if ( this->m_FixedLandmarkContainer->Size() > 0 )
    {
    itkDebugMacro( "Extracting landmark information." );

    typename LandmarkContainerType::Pointer deformedMovingLandmarkContainer
      = LandmarkContainerType::New();
    deformedMovingLandmarkContainer->Initialize();

    for ( unsigned int i = 0; i < this->m_MovingLandmarkContainer->Size(); i++ )
      {
      LandmarkType moving = this->m_MovingLandmarkContainer->GetElement( i );
      typename BSplineControlPointsFilterType::PointDataType point;
      typename BSplineControlPointsFilterType::PointType params;

      for ( unsigned int d = 0; d < ImageDimension; d++ )
        {
        point[d] = moving[d];
        params[d] = ( point[d] - bspliner->GetOrigin()[d] ) /
          ( bspliner->GetSpacing()[d] * ( bspliner->GetSize()[d] - 1 ) );
        }

      bspliner->EvaluateParametersAtPoint( point, params );
      LandmarkType movingWarped;
      for ( unsigned int d = 0; d < ImageDimension; d++ )
        {
        movingWarped[d] = params[d] *
          ( bspliner->GetSpacing()[d] * ( bspliner->GetSize()[d] - 1 ) )
          + bspliner->GetOrigin()[d];
        }
      movingWarped[ImageDimension] = moving[ImageDimension];
      deformedMovingLandmarkContainer->InsertElement( i, movingWarped );
      }

    for ( unsigned int i = 0; i < this->m_FixedLandmarkContainer->Size(); i++ )
      {
      LandmarkType fixed = this->m_FixedLandmarkContainer->GetElement( i );

      RealType min_distance = NumericTraits<RealType>::max();
      int min_idx = -1;

      for ( unsigned int j = 0; j < deformedMovingLandmarkContainer->Size(); j++ )
        {
        LandmarkType moving = deformedMovingLandmarkContainer->GetElement( j );
        if ( static_cast<unsigned int>( fixed[ImageDimension] ) ==
             static_cast<unsigned int>( moving[ImageDimension] ) )
          {
          if ( ( moving - fixed ).GetSquaredNorm() < min_distance )
            {
            min_idx = j;
            min_distance = ( moving - fixed ).GetSquaredNorm();
            }
          }
        }
      if ( min_idx != -1 )
        {
        landmarkEnergy += ( 0.5 * this->m_LandmarkWeighting * min_distance );
        landmarkCount += 1.0;
        }
      }
    }
*/
  /**
   * Calculate contribution from Jacobian constraint
   */
  RealType jacobianEnergy = 0.0;
  RealType jacobianCount = 0.0;

  if ( this->m_EnforceDiffeomorphism )
    {
    typedef VectorFieldGradientImageFunction<DeformationFieldType> FunctionType;
    typename FunctionType::Pointer function = FunctionType::New();

    if ( this->m_CurrentLevel == 0 )
      {
      function->SetInputImage( deformationField );
      }
    else
      {
      typedef BSplineControlPointImageFilter
        <ControlPointLatticeType, DeformationFieldType> BSplineControlPointsFilterType;
      typename BSplineControlPointsFilterType::Pointer bspliner = BSplineControlPointsFilterType::New();
      typename BSplineControlPointsFilterType::ArrayType close;

      close.Fill( false );
      bspliner->SetSplineOrder( this->m_SplineOrder );
      bspliner->SetCloseDimension( close );
      bspliner->SetInput( this->m_CurrentDeformationFieldControlPoints );
      bspliner->SetOrigin( this->m_CurrentFixedImage[0]->GetOrigin() );
      bspliner->SetSize(
        this->m_CurrentFixedImage[0]->GetLargestPossibleRegion().GetSize() );
      bspliner->SetSpacing( this->m_CurrentFixedImage[0]->GetSpacing() );
      bspliner->Update();

      function->SetInputImage( bspliner->GetOutput() );
      }

    typedef DecomposeTensorFunction<typename BSplineFilterType::GradientType> DecomposerType;
    DecomposerType decomposer;

    ImageRegionIteratorWithIndex<FixedImageType> ItP(
      this->m_CurrentFixedImage[0],
      this->m_CurrentFixedImage[0]->GetLargestPossibleRegion() );

    for ( ItP.GoToBegin(); !ItP.IsAtEnd(); ++ItP )
      {
      if ( this->m_WeightImage &&
           ( this->m_CurrentWeightImage->GetPixel( ItP.GetIndex() ) < 0 ||
             ( !this->m_PDEDeformableMetric[1] &&
               this->m_CurrentWeightImage->GetPixel( ItP.GetIndex() ) <= 0 ) ) )
        {
        continue;
        }
/*
      typename BSplineFilterType::GradientType jac;
      bspliner->EvaluateSpatialJacobianAtIndex( ItP.GetIndex(), jac );
      RealType detJ = decomposer.EvaluateDeterminant( jac );
*/
      RealType detJ = function->EvaluateJacobianDeterminantAtIndex( ItP.GetIndex() );
      if ( !vnl_math_isnan( detJ ) )
        {
        if ( detJ < this->m_MinimumJacobian || detJ > this->m_MaximumJacobian )
          {
          jacobianEnergy = NumericTraits<RealType>::max();
          jacobianCount = 1.0;
          break;
          }
        else
          {
/*
          RealType g = vcl_exp( 100*( detJ - this->m_MinimumJacobian ) );
          RealType h = vcl_exp( 100*( this->m_MaximumJacobian - detJ ) );
          jacobianEnergy += 1e10 * this->m_InteriorPenaltyParameter *
            ( 1.0 / ( 1.0 + g ) + 1.0 / ( 1.0 + h ) );
*/
          jacobianEnergy -= this->m_InteriorPenaltyParameter *
            ( vcl_log( detJ - this->m_MinimumJacobian )
            + vcl_log( this->m_MaximumJacobian - detJ ) );
          jacobianCount += 1.0;
          }
        }
      }
    }

//  std::cout << "          Energy (" << metricCount << ", " << landmarkCount << ", " << jacobianCount << "): " << std::flush;
  RealType energy = 0.0;
  if ( metricCount[0] > 0.0 )
    {
    energy += metricEnergy[0] / metricCount[0];
//    std::cout << " metric = " << metricEnergy / metricCount << std::flush;
    }
  if ( metricCount[1] > 0.0 )
    {
    energy += metricEnergy[1] / metricCount[1];
//    std::cout << " metric = " << metricEnergy / metricCount << std::flush;
    }
  if ( landmarkCount > 0.0 )
    {
    energy += landmarkEnergy / landmarkCount;
//    std::cout << " landmark = " << landmarkEnergy / landmarkCount << std::flush;
    }
  if ( jacobianCount > 0.0 )
    {
    energy += jacobianEnergy / jacobianCount;
//    std::cout << " jacobian = " << jacobianEnergy / jacobianCount << std::flush;
    }
//  std::cout << std::endl;


  return energy;
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void
FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::InitializeDeformationFieldWithLandmarks()
{
  typename PointSetType::Pointer fieldPoints = PointSetType::New();
  fieldPoints->Initialize();

  unsigned int index = 0;
  for ( unsigned int i = 0; i < this->m_FixedLandmarkContainer->Size(); i++ )
    {
    LandmarkType fixed = this->m_FixedLandmarkContainer->GetElement( i );

    RealType min_distance = NumericTraits<RealType>::max();
    int min_idx = -1;
    for ( unsigned int j = 0; j < this->m_MovingLandmarkContainer->Size(); j++ )
      {
      LandmarkType moving = this->m_MovingLandmarkContainer->GetElement( j );
      if ( static_cast<unsigned int>( fixed[ImageDimension] ) ==
           static_cast<unsigned int>( moving[ImageDimension] ) )
        {
        if ( ( moving - fixed ).GetSquaredNorm() < min_distance )
          {
          min_idx = j;
          min_distance = ( moving - fixed ).GetSquaredNorm();
          }
        }
      }
    if ( min_idx != -1 )
      {
      LandmarkType moving = this->m_MovingLandmarkContainer->GetElement( min_idx );
      typename PointSetType::PointType point;
      typename DeformationFieldType::PixelType data;
      for ( unsigned int j = 0; j < ImageDimension; j++ )
        {
        point[j] = fixed[j];
        data[j] = moving[j] - fixed[j];
        }
      fieldPoints->SetPoint( index, point );
      fieldPoints->SetPointData( index, data );
      index++;
      }
    }

  typename DeformationFieldType::PointType origin;
  typename DeformationFieldType::SpacingType spacing;
  typename DeformationFieldType::SizeType size;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    origin[i] = this->m_CurrentFixedImage[0]->GetOrigin()[i];
    spacing[i] = this->m_CurrentFixedImage[0]->GetSpacing()[i];
    size[i] = this->m_CurrentFixedImage[0]->GetLargestPossibleRegion().GetSize()[i];
    }

  typename BSplineFilterType::ArrayType close;
  close.Fill( false );
  typename BSplineFilterType::ArrayType nlevels;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    unsigned int totalNumberOfBSplineLevels
      = static_cast<unsigned int>(
          floor( log( this->m_PyramidLevelImageSizes[this->m_NumberOfLevels-1][i]
                   /this->m_MeshResolution[i] )/vnl_math::ln2 ) );
    nlevels[i] = totalNumberOfBSplineLevels - this->m_NumberOfLevels + 2;
    }
  typename BSplineFilterType::ArrayType ncps;
  ncps.Fill( this->m_SplineOrder + 1 );

  typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();

//  bspliner->DebugOn();
  bspliner->SetOrigin( origin );
  bspliner->SetSpacing( spacing );
  bspliner->SetSize( size );
  bspliner->SetNumberOfLevels( nlevels );
  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetNumberOfControlPoints( ncps );
  bspliner->SetCloseDimension( close );
  bspliner->SetInput( fieldPoints );
  bspliner->SetGenerateOutputImage( false );
  bspliner->Update();

  this->m_TotalDeformationFieldControlPoints = bspliner->GetPhiLattice();
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
typename FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>::RealType
FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::CalculateLandmarkError()
{
  RealType error = 0.0;

  if ( this->m_FixedLandmarkContainer->Size() > 0 )
    {

    typedef BSplineControlPointImageFilter
      <ControlPointLatticeType, DeformationFieldType> BSplineControlPointsFilterType;
    typename BSplineControlPointsFilterType::Pointer bspliner = BSplineControlPointsFilterType::New();
    typename BSplineControlPointsFilterType::ArrayType close;

    close.Fill( false );
    bspliner->SetSplineOrder( this->m_SplineOrder );
    bspliner->SetCloseDimension( close );
    bspliner->SetInput( this->m_TotalDeformationFieldControlPoints );
    bspliner->SetOrigin( this->m_CurrentFixedImage[0]->GetOrigin() );
    bspliner->SetSize(
      this->m_CurrentFixedImage[0]->GetLargestPossibleRegion().GetSize() );
    bspliner->SetSpacing( this->m_CurrentFixedImage[0]->GetSpacing() );
    /* Deform all the moving landmarks according the the deformation field  */

    typename LandmarkContainerType::Pointer deformedMovingLandmarkContainer
      = LandmarkContainerType::New();
    deformedMovingLandmarkContainer->Initialize();

    for ( unsigned int i = 0; i < this->m_MovingLandmarkContainer->Size(); i++ )
      {
      LandmarkType moving = this->m_MovingLandmarkContainer->GetElement( i );
      typename BSplineControlPointsFilterType::PointDataType point;
      typename BSplineControlPointsFilterType::PointType params;

      for ( unsigned int d = 0; d < ImageDimension; d++ )
        {
        point[d] = moving[d];
        params[d] = ( point[d] - bspliner->GetOrigin()[d] ) /
          ( bspliner->GetSpacing()[d] * ( bspliner->GetSize()[d] - 1 ) );
        }

      bspliner->EvaluateParametersAtPoint( point, params );
      LandmarkType movingWarped;
      for ( unsigned int d = 0; d < ImageDimension; d++ )
        {
        movingWarped[d] = params[d] *
          ( bspliner->GetSpacing()[d] * ( bspliner->GetSize()[d] - 1 ) )
          + bspliner->GetOrigin()[d];
        }
      movingWarped[ImageDimension] = moving[ImageDimension];
      deformedMovingLandmarkContainer->InsertElement( i, movingWarped );
      }

    for ( unsigned int i = 0; i < this->m_FixedLandmarkContainer->Size(); i++ )
      {
      LandmarkType fixed = this->m_FixedLandmarkContainer->GetElement( i );

      RealType min_distance = NumericTraits<RealType>::max();
      int min_idx = -1;

      for ( unsigned int j = 0; j < deformedMovingLandmarkContainer->Size(); j++ )
        {
        LandmarkType moving = deformedMovingLandmarkContainer->GetElement( j );
        if ( static_cast<unsigned int>( fixed[ImageDimension] ) ==
             static_cast<unsigned int>( moving[ImageDimension] ) )
          {
          if ( ( moving - fixed ).GetSquaredNorm() < min_distance )
            {
            min_idx = j;
            min_distance = ( moving - fixed ).GetSquaredNorm();
            }
          }
        }
      if ( min_idx != -1 )
        {
        LandmarkType moving = deformedMovingLandmarkContainer->GetElement( min_idx );
        typename PointSetType::PointType point;
        typename DeformationFieldType::PixelType data;
        for ( unsigned int j = 0; j < ImageDimension; j++ )
          {
          point[j] = fixed[j];
          data[j] = moving[j] - fixed[j];
          }
        error += data.GetNorm();
        }
      }
    }
  return error;
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
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

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void
FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::LineMinimization( RealType *step, RealType *fret )
{
  std::cout << "    Begin line search..." << std::endl;

  // We should now have a, b and c, as well as f(a), f(b), f(c),
  // where b gives the minimum energy position;
  RealType ax, bx, cx;
  this->FindBracketingTriplet( &ax, &bx, &cx );

  this->BrentSearch( ax, bx, cx, step, fret );
//  this->GoldenSectionSearch( ax, bx, cx, step, fret );
 //this->BruteForceSearch( -0.5, 0.0, 1.5, step, fret );
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
typename FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>::RealType
FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::EvaluateEnergyForLineSearch( RealType lambda )
{
  return this->EvaluateMetricOverImageRegion( lambda );
}

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void
FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
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

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void
FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::BrentSearch( RealType ax, RealType bx, RealType cx, RealType *step, RealType *fret )
{
  const RealType R = 0.6180339;
  const RealType CGOLD = 1.0 - R;
  const RealType tol = 0.1;
  const RealType ZEPS = 1.0e-10;

  RealType a;
  RealType b;
  RealType d = 0.0;
  RealType p;
  RealType q;
  RealType r;
  RealType etemp;
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

  for ( unsigned int iter = 1; iter <= this->m_LineSearchMaximumIterations; iter++ )
    {
    std::cout << "        Iteration (Brent's) " << iter
              << ": f(x = " << x << ") = " << fx << ", x in ["
              << a << ", " << b << "] " << std::endl;

    RealType xm = 0.5 * ( a + b );
    RealType tol1 = tol * vnl_math_abs( x ) + ZEPS;
    RealType tol2 = 2.0 * tol1;
    if ( vnl_math_abs( x - xm ) <= ( tol2 - 0.5 * ( b - a ) ) )
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
template<class TMovingImage, class TFixedImage, class TWarpedImage>
void
FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
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
*/

template<class TMovingImage, class TFixedImage, class TWarpedImage>
void
FFDRegistrationFilter<TMovingImage, TFixedImage, TWarpedImage>
::BruteForceSearch( RealType ax, RealType bx, RealType cx, RealType *step, RealType *fret )
{
  static unsigned int whichSearch = 0;
  const unsigned int numberOfIntervals = 10;


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
  str << minStep << " " << minValue << " " << whichSearch << " " << this->m_CurrentLevel << std::endl;
  for ( RealType x = a; x <= b; x += dx )
    {
    RealType fx = this->EvaluateEnergyForLineSearch( x );
    str << x << " " << fx << " " << whichSearch << " " << this->m_CurrentLevel << std::endl;
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

}
#endif

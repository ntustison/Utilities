/*=========================================================================

  Program:   Advanced Normalization Tools
  Module:    $RCSfile: itkLabelOverlapMeasuresImageFilter.txx,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) ConsortiumOfANTS. All rights reserved.
  See accompanying COPYING.txt or
 http://sourceforge.net/projects/advants/files/ANTS/ANTSCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDiReCTImageFilter_txx
#define __itkDiReCTImageFilter_txx

#include "itkDiReCTImageFilter.h"

#include "itkBinaryContourImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkMaximumImageFilter.h"
#include "itkMultiplyByConstantImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkWarpImageFilter.h"

namespace itk
{

template<class TInputImage, class TOutputImage>
DiReCTImageFilter<TInputImage, TOutputImage>
::DiReCTImageFilter() :
  m_MaximumNumberOfIterations( 50 ),
  m_ThicknessPriorEstimate( 6.0 ),
  m_SmoothingSigma( 1.0 ),
  m_GradientStep( 0.5 ),
  m_NumberOfIntegrationPoints( 10 ),
  m_GrayMatterLabel( 2 ),
  m_WhiteMatterLabel( 3 )
{
  this->SetNumberOfRequiredInputs( 3 );
}

template<class TInputImage, class TOutputImage>
DiReCTImageFilter<TInputImage, TOutputImage>
::~DiReCTImageFilter()
{
}

template<class TInputImage, class TOutputImage>
void
DiReCTImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  // Estimation of cortical thickness involves the following steps:
  //   1. Diffuse the white matter region
  //   2. Calculate gradient of output of step 1.
  //   3.

  // Diffuse the white matter region and calculate the gradient

  RealImagePointer diffusedWhiteMatter =
    this->DiffuseWhiteMatterRegion( this->GetSegmentationImage(), 100 );

  // Dilate the gray matter

  InputImagePointer dilatedGrayMatter = this->DilateRegion(
    this->GetSegmentationImage(), this->m_GrayMatterLabel, 1 );

  // Extract the gray and white matter contours

  InputImagePointer grayMatterContours = this->ExtractRegionalContours(
    this->GetSegmentationImage(), this->m_GrayMatterLabel );
  InputImagePointer whiteMatterContours = this->ExtractRegionalContours(
    this->GetSegmentationImage(), this->m_WhiteMatterLabel );

  // Iterate

    // Initialize fields

  VectorType zeroVector( 0.0 );

  VectorImagePointer forwardField = VectorImageType::New();
  forwardField->CopyInformation( this->GetInput() );
  forwardField->SetRegions( this->GetInput()->GetRequestedRegion() );
  forwardField->Allocate();

  VectorImagePointer forwardIncrementalField = VectorImageType::New();
  forwardIncrementalField->CopyInformation( this->GetInput() );
  forwardIncrementalField->SetRegions( this->GetInput()->GetRequestedRegion() );
  forwardIncrementalField->Allocate();

  VectorImagePointer inverseField = VectorImageType::New();
  inverseField->CopyInformation( this->GetInput() );
  inverseField->SetRegions( this->GetInput()->GetRequestedRegion() );
  inverseField->Allocate();

  VectorImagePointer inverseIncrementalField = VectorImageType::New();
  inverseIncrementalField->CopyInformation( this->GetInput() );
  inverseIncrementalField->SetRegions( this->GetInput()->GetRequestedRegion() );
  inverseIncrementalField->Allocate();


  unsigned int iterations = 0;
  while( iterations++ < this->m_MaximumNumberOfIterations ) // && badct < 4 )
    {
    forwardField->FillBuffer( zeroVector );
    inverseField->FillBuffer( zeroVector );
    inverseIncrementalField->FillBuffer( zeroVector );

    unsigned int integrationPoints = 0;
    while( integrationPoints++ < this->m_NumberOfIntegrationPoints )
      {
      inverseField = this->ComposeDiffeomorphisms( inverseField,
        inverseIncrementalField );

   	  RealImagePointer warpedWhiteMatter = this->WarpImage(
   	    this->GetWhiteMatterProbabilityImage(), inverseField );
   	  RealImagePointer warpedThicknessImage = this->WarpImage(
   	    thicknessImage, inverseField );
   	  RealImagePointer warpedWhiteMatterContours = this->WarpImage(
   	    whiteMatterContours, inverseField );

      typedef GradientRecursiveGaussianImageFilter<RealImageType, VectorImageType>
        GradientImageFilterType;
      typename GradientImageFilterType::Pointer gradient =
        GradientImageFilterType::New();
      gradient->SetInput( warpedWhiteMatterContours );
      gradient->SetSigma( this->m_SmoothingSigma );
      gradient->Update();

      // Generate speed image

      RealImagePointer speedImage = RealImageType::New();
      speedImage->CopyInformation( this->GetInput() );
      speedImage->SetRegions( this->GetInput()->GetRequestedRegion() );
      speedImage->Allocate();
      speedImage->FillBuffer( 0.0 );

      ImageRegionIteratorWithIndex<RealImageType> ItS( speedImage,
        speedImage->GetRequestedRegion() );
      ImageRegionConstIterator<InputImageType> It( this->GetSegmentationImage(),
        this->GetSegmentationImage()->GetRequestedRegion() );
      ImageRegionIterator<VectorImageType> ItG( gradient->GetOutput(),
        gradient->GetOutput()->GetRequestedRegion() );

      It.GoToBegin();
      ItS.GoToBegin();
      ItG.GoToBegin();

      while( !It.IsAtEnd() )
        {
        if( It.Get() == this->m_GrayMatterLabel )
          {
          RealType norm = ( ItG.Get() ).GetNorm();
          if( norm > 0.0 && !vnl_math_isnan( norm ) && !vnl_math_isinf( norm ) )
            {
            ItG.Set( ItG.Get() / norm );
            }
          else
            {
            ItG.Set( zeroVector );
            }
          }
        RealType difference = 1.0 - warpedWhiteMatterContours;
 	      totalError += vnl_math_abs( difference );

        ItS.Set( difference );

        ++It;
        ++ItS;
        ++ItG;
        }

      // Calculate objective function value

      ImageRegionIterator<VectorImageType> ItI( incrementalField,
        incrementalField->GetRequestedRegion() );

      It.GoToBegin();
      ItS.GoToBegin();
      ItG.GoToBegin();
      ItI.GoToBegin();

      while( !It.IsAtEnd() )
        {
        ItI.Set( ItI.Get() + ItG.Get() * ItS.Get() );
        }
	     }
    }

  typedef CastImageFilter<InputImageType, OutputImageType> CasterType;
  typename CasterType::Pointer caster = CasterType::New();
  caster->SetInput( grayMatterContours );
  caster->Update();

  this->SetNthOutput( 0, caster->GetOutput() );
//  this->SetNthOutput( 0, diffusedWhiteMatter );
}

template<class TInputImage, class TOutputImage>
typename DiReCTImageFilter<TInputImage, TOutputImage>::RealImagePointer
DiReCTImageFilter<TInputImage, TOutputImage>
::DiffuseWhiteMatterRegion( const InputImageType *segmentationImage,
  unsigned int numberOfIterations )
{
  typedef BinaryThresholdImageFilter<InputImageType, RealImageType>
    ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( segmentationImage );
  thresholder->SetLowerThreshold( this->m_WhiteMatterLabel );
  thresholder->SetUpperThreshold( this->m_WhiteMatterLabel );
  thresholder->SetInsideValue( 1 );
  thresholder->SetOutsideValue( 0 );
  thresholder->Update();

  typedef ImageDuplicator<RealImageType> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage( thresholder->GetOutput() );
  duplicator->Update();
  RealImagePointer diffusedWhiteMatter = duplicator->GetOutput();

  typedef DiscreteGaussianImageFilter<RealImageType, RealImageType> SmootherType;
  typedef MaximumImageFilter<RealImageType, RealImageType, RealImageType>
    MaxFilterType;

  for( unsigned int n = 0; n < numberOfIterations; n++ )
    {
    typename SmootherType::Pointer smoother = SmootherType::New();
    smoother->SetVariance( vnl_math_sqr( this->m_SmoothingSigma ) );
    smoother->SetUseImageSpacingOn();
    smoother->SetMaximumError( 0.01 );
    smoother->SetInput( diffusedWhiteMatter );
    smoother->Update();

    typename MaxFilterType::Pointer maxFilter = MaxFilterType::New();
    maxFilter->SetInput1( smoother->GetOutput() );
    maxFilter->SetInput2( thresholder->GetOutput() );

    diffusedWhiteMatter = maxFilter->GetOutput();
    diffusedWhiteMatter->Update();
    diffusedWhiteMatter->DisconnectPipeline();
    }

  return diffusedWhiteMatter;
}

template<class TInputImage, class TOutputImage>
typename DiReCTImageFilter<TInputImage, TOutputImage>::InputImagePointer
DiReCTImageFilter<TInputImage, TOutputImage>
::DilateRegion( const InputImageType *segmentationImage,
  unsigned int whichRegion, unsigned int radius )
{
  typedef BinaryThresholdImageFilter<InputImageType, InputImageType>
    ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( segmentationImage );
  thresholder->SetLowerThreshold( whichRegion );
  thresholder->SetUpperThreshold( whichRegion );
  thresholder->SetInsideValue( 1 );
  thresholder->SetOutsideValue( 0 );
  thresholder->Update();

  typedef BinaryBallStructuringElement<InputPixelType, ImageDimension>
    StructuringElementType;
  typedef BinaryDilateImageFilter<InputImageType, InputImageType,
    StructuringElementType> DilatorType;

  StructuringElementType  structuringElement;
  structuringElement.SetRadius( radius );
  structuringElement.CreateStructuringElement();

  typename DilatorType::Pointer dilator = DilatorType::New();
  dilator->SetInput( thresholder->GetOutput() );
  dilator->SetKernel( structuringElement );
  dilator->SetDilateValue( 1 );

  InputImagePointer dilatedRegion = dilator->GetOutput();
  dilatedRegion->Update();
  dilatedRegion->DisconnectPipeline();
  dilatedRegion->SetRegions( segmentationImage->GetRequestedRegion() );

  return dilatedRegion;
}

template<class TInputImage, class TOutputImage>
typename DiReCTImageFilter<TInputImage, TOutputImage>::InputImagePointer
DiReCTImageFilter<TInputImage, TOutputImage>
::ExtractRegionalContours( const InputImageType *segmentationImage,
  unsigned int whichRegion )
{
  typedef BinaryThresholdImageFilter<InputImageType, InputImageType>
    ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( segmentationImage );
  thresholder->SetLowerThreshold( whichRegion );
  thresholder->SetUpperThreshold( whichRegion );
  thresholder->SetInsideValue( 1 );
  thresholder->SetOutsideValue( 0 );
  thresholder->Update();

  typedef BinaryContourImageFilter<InputImageType, InputImageType>
    ContourFilterType;
  typename ContourFilterType::Pointer contourFilter = ContourFilterType::New();
  contourFilter->SetInput( thresholder->GetOutput() );
  contourFilter->SetFullyConnected( true );
  contourFilter->SetBackgroundValue( 0 );
  contourFilter->SetForegroundValue( 1 );

  InputImagePointer contours = contourFilter->GetOutput();
  contours->Update();
  contours->DisconnectPipeline();
  contours->SetRegions( segmentationImage->GetRequestedRegion() );

  return contours;
}

template<class TInputImage, class TOutputImage>
typename DiReCTImageFilter<TInputImage, TOutputImage>::VectorImagePointer
DiReCTImageFilter<TInputImage, TOutputImage>
::ComposeDiffeomorphisms( const VectorImageType *inputField,
  const VectorImageType *warp )
{
  VectorType zeroVector( 0.0 );

  VectorImagePointer outputField = VectorImageType::New();
  outputField->CopyInformation( warp );
  outputField->SetRegions( warp->GetRequestedRegion() );
  outputField->Allocate();
  outputField->FillBuffer( zeroVector );

  typedef VectorLinearInterpolateImageFunction<VectorImageType, RealType>
    InterpolatorType;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage( inputField );

  ImageRegionConstIteratorWithIndex<VectorImageType> ItW( warp,
    warp->GetRequestedRegion() );
  ImageRegionIterator<VectorImageType> ItF( outputField,
    outputField->GetRequestedRegion() );
  for( ItW.GoToBegin(), ItF.GoToBegin(); !ItW.IsAtEnd(); ++ItW, ++ItF )
    {
    PointType point1;
    warp->TransformIndexToPhysicalPoint( ItW.GetIndex(), point1 );

    PointType point2 = point1 + ItW.Get();

    typename InterpolatorType::OutputType displacement;
    if( interpolator->IsInsideBuffer( point2 ) )
      {
      displacement = interpolator->Evaluate( point2 );
      }
    else
      {
      displacement.Fill( 0.0 );
      }

    ItF.Set( ( point2 + displacement ) - point1 );
    }
}

template<class TInputImage, class TOutputImage>
typename DiReCTImageFilter<TInputImage, TOutputImage>::RealImagePointer
DiReCTImageFilter<TInputImage, TOutputImage>
::WarpImage( const RealImageType *inputImage,
  const VectorImageType *deformationField )
{
  typedef WarpImageFilter<RealImageType, RealImageType, VectorImageType>
    WarperType;
  typename WarperType::Pointer warper = WarperType::New();
  warper->SetInput( inputImage );
  warper->SetDeformationField( deformationField );
  warper->SetEdgePaddingValue( 0 );
  warper->SetOutputSpacing( inputImage->GetSpacing() );
  warper->SetOutputOrigin( inputImage->GetOrigin() );
  warper->SetOutputDirection( inputImage->GetDirection() );

  RealImagePointer warpedImage = warper->GetOutput();
  warpedImage->Update();
  warpedImage->DisconnectPipeline();

  return warpedImage;
}

template<class TInputImage, class TOutputImage>
void
DiReCTImageFilter<TInputImage, TOutputImage>
::InvertDeformationField( const VectorImageType *deformationField,
  VectorImageType *inverseField )
{
  typename VectorImageType::SpacingType spacing =
    deformationField->GetSpacing();

  RealType maxNorm = 1.0;
  RealType meanNorm = 1.0;
  unsigned int iterations = 0;
  while( iterations++ < 20 && maxNorm > 0.1 && meanNorm > 0.001 )
    {
    meanNorm = 0.0;
    maxNorm = 0.0;

    VectorImagePointer eulerianField =
      this->ComposeDiffeomorphisms( inverseField, deformationField );

    ImageRegionIterator<VectorImageType> ItE( eulerianField,
      eulerianField->GetRequestedRegion() );
    for( ItE.GoToBegin(); !ItE.IsAtEnd(); ++ItE )
      {
     	VectorType update = ItE.Get();
     	for( unsigned int d = 0; d < ImageDimension; d++ )
	       {
   	    update[d] *= -1.0 / spacing[d];
        }
      RealType norm = update.GetNorm();
     	meanNorm += norm;
     	if( norm > maxNorm )
     	  {
     	  maxNorm = norm;
     	  }
     	ItE.Set( update );
      }
    meanNorm /= static_cast<RealType>(
      eulerianField->GetRequestedRegion().GetNumberOfPixels() );

    ImageRegionIterator<VectorImageType> ItI( inverseField,
      inverseField->GetRequestedRegion() );
    for( ItI.GoToBegin(), ItE.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItE )
      {
      VectorType update = ItE.Get();
      RealType updateNorm = update.GetNorm();

      if( updateNorm > 0.5 * maxNorm )
        {
        update *= ( 0.5 * maxNorm / updateNorm );
        }
      ItI.Set( ItI.Get() + update * 0.5 );
      }
    }
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TOutputImage>
void
DiReCTImageFilter<TInputImage, TOutputImage>
::PrintSelf( std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  std::cout << indent << "Maximum number of iterations = "
    << this->m_MaximumNumberOfIterations << std::endl;
  std::cout << indent << "Thickness prior estimate = "
    << this->m_ThicknessPriorEstimate << std::endl;
  std::cout << indent << "Smoothing sigma = "
    << this->m_SmoothingSigma << std::endl;
  std::cout << indent << "Gradient step = "
    << this->m_GradientStep << std::endl;
}

} // end namespace itk

#endif
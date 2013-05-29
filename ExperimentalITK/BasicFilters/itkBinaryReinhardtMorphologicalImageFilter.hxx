/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBinaryReinhardtMorphologicalImageFilter.hxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:50 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBinaryReinhardtMorphologicalImageFilter_hxx
#define __itkBinaryReinhardtMorphologicalImageFilter_hxx

#include "itkBinaryReinhardtMorphologicalImageFilter.h"

#include "itkBinaryBoundedSpaceDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryBoxStructuringElement.h"
#include "itkBinaryDiamondStructuringElement.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkStatisticsImageFilter.h"

namespace itk
{

template <class TInputImage, class TOutputImage, class TKernel>
BinaryReinhardtMorphologicalImageFilter<TInputImage, TOutputImage, TKernel>
::BinaryReinhardtMorphologicalImageFilter()
{
  this->m_EmploySaltAndPepperRepair = true;
  this->m_SaltAndPepperMinimumSizeInPixels = 5;

  this->m_EmployMinimumDiameterFilter = true;
  this->m_MinimumDiameterStructuringElementRadius = 1;

  this->m_EmployUnwantedCavityDeletion = true;

  this->m_EmployMinimumSizeFilter = true;
  this->m_MinimumSizeStructuringElementRadius = 1;

  this->m_EmployMaximumDiameterFilter = true;
  this->m_MaximumDiameterStructuringElementRadius = 1;

  this->m_EmployConnectivityFilter = true;
  this->m_NumberOfConnectedComponents = 1;

  this->m_EmployBoundarySmoother = true;
  this->m_BoundarySmootherStructuringElementRadius = 1;

   this->m_EmployUnclassifiedPixelProcessing = true;
}

template< class TInputImage, class TOutputImage, class TKernel>
void
BinaryReinhardtMorphologicalImageFilter< TInputImage, TOutputImage, TKernel>
::GenerateData()
{
  typedef CastImageFilter<InputImageType, OutputImageType> CasterType;
  typename CasterType::Pointer caster = CasterType::New();
  caster->SetInput( this->GetInput() );
  caster->Update();

  typename OutputImageType::Pointer image = caster->GetOutput();

  if ( this->m_EmploySaltAndPepperRepair )
    {
    this->SaltAndPepperRepair( image );
    }
  if ( this->m_EmployMinimumDiameterFilter )
    {
    this->MinimumDiameterFilter( image );
    }
  if ( this->m_EmployUnwantedCavityDeletion )
    {
    this->UnwantedCavityDeletion( image );
    }
  if ( this->m_EmployMinimumSizeFilter )
    {
    this->MinimumSizeFilter( image );
    }
  if ( this->m_EmployMaximumDiameterFilter )
    {
    this->MaximumDiameterFilter( image );
    }
  if ( this->m_EmployConnectivityFilter )
    {
    this->ConnectivityFilter( image );
    }
  if ( this->m_EmployBoundarySmoother )
    {
    this->BoundarySmoother( image );
    }
  if ( this->m_EmployUnclassifiedPixelProcessing )
    {
    this->UnclassifiedPixelProcessing( image );
    }

  this->GraftOutput( image );
}

template <class TInputImage, class TOutput, class TKernel>
void
BinaryReinhardtMorphologicalImageFilter<TInputImage, TOutput, TKernel>
::SaltAndPepperRepair( typename OutputImageType::Pointer image )
{
  typedef BinaryThresholdImageFilter<OutputImageType, OutputImageType> ThresholderType;
  typename ThresholderType::Pointer thresholder1 = ThresholderType::New();
  thresholder1->SetInput( image );
  thresholder1->SetLowerThreshold(
    static_cast<OutputPixelType>( this->GetForegroundValue() ) );
  thresholder1->SetUpperThreshold(
    static_cast<OutputPixelType>( this->GetForegroundValue() ) );
  thresholder1->SetInsideValue( NumericTraits<OutputPixelType>::One );
  thresholder1->SetOutsideValue( NumericTraits<OutputPixelType>::Zero );
  thresholder1->Update();

  typedef ConnectedComponentImageFilter<OutputImageType, OutputImageType> ConnectedComponentType;
  typename ConnectedComponentType::Pointer connecter1 = ConnectedComponentType::New();
  connecter1->SetInput( thresholder1->GetOutput() );
  connecter1->FullyConnectedOff();
  connecter1->Update();

  typedef RelabelComponentImageFilter<OutputImageType, OutputImageType> LabelerType;
  typename LabelerType::Pointer labeler1 = LabelerType::New();
  labeler1->SetInput( connecter1->GetOutput() );
  labeler1->SetMinimumObjectSize( this->m_SaltAndPepperMinimumSizeInPixels );
  labeler1->Update();

  typename ThresholderType::Pointer thresholder2 = ThresholderType::New();
  thresholder2->SetInput( labeler1->GetOutput() );
  thresholder2->SetLowerThreshold( NumericTraits<OutputPixelType>::One );
  if ( labeler1->GetNumberOfObjects() > 1 )
    {
    thresholder2->SetUpperThreshold(
      static_cast<OutputPixelType>( labeler1->GetNumberOfObjects() ) );
    }
  else
    {
    thresholder2->SetUpperThreshold( NumericTraits<OutputPixelType>::One );
    }
  thresholder2->SetInsideValue( NumericTraits<OutputPixelType>::Zero );
  thresholder2->SetOutsideValue( NumericTraits<OutputPixelType>::One );
  thresholder2->Update();

  typename ConnectedComponentType::Pointer connecter2 = ConnectedComponentType::New();
  connecter2->SetInput( thresholder2->GetOutput() );
  connecter2->FullyConnectedOff();
  connecter2->Update();

  typename LabelerType::Pointer labeler2 = LabelerType::New();
  labeler2->SetInput( connecter2->GetOutput() );
  labeler2->SetMinimumObjectSize( this->m_SaltAndPepperMinimumSizeInPixels );
  labeler2->Update();

  ImageRegionIterator<OutputImageType> ItI( image,
    image->GetRequestedRegion() );
  ImageRegionIterator<OutputImageType> ItL( labeler2->GetOutput(),
    labeler2->GetOutput()->GetRequestedRegion() );
  for ( ItI.GoToBegin(), ItL.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItL )
    {
    if ( ItL.Get() == 0 )
      {
      ItI.Set( static_cast<OutputPixelType>( this->GetForegroundValue() ) );
      }
    }
}

template <class TInputImage, class TOutput, class TKernel>
void
BinaryReinhardtMorphologicalImageFilter<TInputImage, TOutput, TKernel>
::MinimumDiameterFilter( typename OutputImageType::Pointer image )
{
  typedef BinaryThresholdImageFilter<OutputImageType, OutputImageType> ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( image );
  thresholder->SetLowerThreshold(
    static_cast<OutputPixelType>( this->GetForegroundValue() ) );
  thresholder->SetUpperThreshold(
    static_cast<OutputPixelType>( this->GetForegroundValue() ) );
  thresholder->SetInsideValue( NumericTraits<OutputPixelType>::One );
  thresholder->SetOutsideValue( NumericTraits<OutputPixelType>::Zero );
  thresholder->Update();

  typedef BinaryBallStructuringElement<OutputPixelType,
    OutputImageDimension> BallStructuringElementType;
  BallStructuringElementType ballStructuringElement;
  ballStructuringElement.SetRadius( this->m_MinimumDiameterStructuringElementRadius );
  ballStructuringElement.CreateStructuringElement();

  typedef BinaryErodeImageFilter<OutputImageType, OutputImageType,
    BallStructuringElementType> EroderType;
  typename EroderType::Pointer eroder = EroderType::New();
  eroder->SetInput( thresholder->GetOutput() );
  eroder->SetKernel( ballStructuringElement );
  eroder->SetForegroundValue( NumericTraits<OutputPixelType>::One );
  eroder->SetBackgroundValue( NumericTraits<OutputPixelType>::Zero );
  eroder->Update();

  typedef ConnectedComponentImageFilter<OutputImageType, OutputImageType> ConnectedComponentType;
  typename ConnectedComponentType::Pointer connecter = ConnectedComponentType::New();
  connecter->SetInput( eroder->GetOutput() );
  connecter->FullyConnectedOff();
  connecter->Update();

  typedef RelabelComponentImageFilter<OutputImageType, OutputImageType> LabelerType;
  typename LabelerType::Pointer labeler = LabelerType::New();
  labeler->SetInput( connecter->GetOutput() );
  labeler->Update();

  ImageRegionIterator<OutputImageType> It( labeler->GetOutput(),
    labeler->GetOutput()->GetRequestedRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( It.Get() > static_cast<OutputPixelType>(
         this->m_NumberOfConnectedComponents ) || It.Get() == 0 )
      {
      It.Set( NumericTraits<OutputPixelType>::Zero );
      }
    else
      {
      It.Set( NumericTraits<OutputPixelType>::One );
      }
    }

  typedef BinaryDilateImageFilter<OutputImageType, OutputImageType,
    BallStructuringElementType> DilaterType;
  typename DilaterType::Pointer dilater = DilaterType::New();
  dilater->SetInput( labeler->GetOutput() );
  dilater->SetKernel( ballStructuringElement );
  dilater->SetForegroundValue( NumericTraits<OutputPixelType>::One );
  dilater->SetBackgroundValue( NumericTraits<OutputPixelType>::Zero );
  dilater->Update();

  typedef BinaryBoxStructuringElement<OutputPixelType,
    OutputImageDimension> BoxStructuringElementType;
  BoxStructuringElementType boxStructuringElement;
  boxStructuringElement.SetRadius( 1 );
  boxStructuringElement.CreateStructuringElement();

  typedef BinaryDilateImageFilter<OutputImageType, OutputImageType,
    BoxStructuringElementType> BoxDilaterType;
  typename BoxDilaterType::Pointer boxDilater = BoxDilaterType::New();
  boxDilater->SetInput( dilater->GetOutput() );
  boxDilater->SetKernel( boxStructuringElement );
  boxDilater->SetForegroundValue( NumericTraits<OutputPixelType>::One );
  boxDilater->SetBackgroundValue( NumericTraits<OutputPixelType>::Zero );
  boxDilater->Update();

  ImageRegionIterator<OutputImageType> ItB( boxDilater->GetOutput(),
    boxDilater->GetOutput()->GetRequestedRegion() );
  ImageRegionIterator<OutputImageType> ItI( image,
    image->GetRequestedRegion() );
  for ( ItB.GoToBegin(), ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItB, ++ItI )
    {
    if ( ItB.Get() != NumericTraits<OutputPixelType>::One &&
         ItI.Get() == static_cast<OutputPixelType>( this->GetForegroundValue() ) )
      {
      ItI.Set( this->GetBackgroundValue() );
      ItB.Set( NumericTraits<OutputPixelType>::One );
      }
    else
      {
      ItB.Set( NumericTraits<OutputPixelType>::Zero );
      }
    }

  typedef StatisticsImageFilter<OutputImageType> StatsType;
  typename StatsType::Pointer stats = StatsType::New();
  stats->SetInput( image );
  stats->Update();

  /**
   * Add attachments to other regions
   */
  for ( unsigned int i = 1; i <= static_cast<unsigned int>(
          stats->GetMaximum() ); i++ )
    {
    typedef BinaryThresholdImageFilter<OutputImageType, OutputImageType> ThresholderType;
    typename ThresholderType::Pointer thresholder = ThresholderType::New();
    thresholder->SetInput( image );
    thresholder->SetOutsideValue( NumericTraits<OutputPixelType>::Zero );
    thresholder->SetInsideValue( NumericTraits<OutputPixelType>::One );
    thresholder->SetLowerThreshold( static_cast<OutputPixelType>( i ) );
    thresholder->SetUpperThreshold( static_cast<OutputPixelType>( i ) );
    thresholder->Update();

    typedef BinaryMorphologicalClosingImageFilter<OutputImageType, OutputImageType,
      BallStructuringElementType> CloserType;
    typename CloserType::Pointer closer = CloserType::New();
    closer->SetInput( thresholder->GetOutput() );
    closer->SetKernel( ballStructuringElement );
    closer->SetForegroundValue( NumericTraits<OutputPixelType>::One );
    closer->Update();

    ImageRegionIterator<OutputImageType> ItC( closer->GetOutput(),
      closer->GetOutput()->GetRequestedRegion() );
    ImageRegionIterator<OutputImageType> ItI( image,
      image->GetRequestedRegion() );
    ItC.GoToBegin();
    ItI.GoToBegin();
    ItB.GoToBegin();
    while ( !ItC.IsAtEnd() )
      {
      if ( ( ItB.Get() == NumericTraits<OutputPixelType>::One &&
           ItC.Get() == NumericTraits<OutputPixelType>::One ) ||
           ( ItI.Get() == static_cast<OutputPixelType>( i ) ) )
        {
        ItI.Set( static_cast<OutputPixelType>( i ) );
        }

      ++ItC;
      ++ItI;
      ++ItB;
      }
    }
}

template <class TInputImage, class TOutput, class TKernel>
void
BinaryReinhardtMorphologicalImageFilter<TInputImage, TOutput, TKernel>
::UnwantedCavityDeletion( typename OutputImageType::Pointer image )
{
  typedef BinaryThresholdImageFilter<OutputImageType, OutputImageType> ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( image );
  thresholder->SetLowerThreshold(
    static_cast<OutputPixelType>( this->GetForegroundValue() ) );
  thresholder->SetUpperThreshold(
    static_cast<OutputPixelType>( this->GetForegroundValue() ) );
  thresholder->SetInsideValue( NumericTraits<OutputPixelType>::Zero );
  thresholder->SetOutsideValue( NumericTraits<OutputPixelType>::One );
  thresholder->Update();

  typedef ConnectedComponentImageFilter<OutputImageType, OutputImageType> ConnectedComponentType;
  typename ConnectedComponentType::Pointer connecter = ConnectedComponentType::New();
  connecter->SetInput( thresholder->GetOutput() );
  connecter->FullyConnectedOff();
  connecter->Update();

  typedef RelabelComponentImageFilter<OutputImageType, OutputImageType> LabelerType;
  typename LabelerType::Pointer labeler = LabelerType::New();
  labeler->SetInput( connecter->GetOutput() );
  labeler->Update();

  ImageRegionIterator<OutputImageType> ItL( labeler->GetOutput(),
    labeler->GetOutput()->GetRequestedRegion() );
  ImageRegionIterator<OutputImageType> ItI( image,
    image->GetRequestedRegion() );
  for ( ItL.GoToBegin(), ItI.GoToBegin(); !ItL.IsAtEnd(); ++ItL, ++ItI )
    {
    if ( ItL.Get() > 1 || ItL.Get() == 0 )
      {
      ItI.Set( static_cast<OutputPixelType>( this->GetForegroundValue() ) );
      }
    }
}

template <class TInputImage, class TOutput, class TKernel>
void
BinaryReinhardtMorphologicalImageFilter<TInputImage, TOutput, TKernel>
::MinimumSizeFilter( typename OutputImageType::Pointer image )
{
  typedef BinaryThresholdImageFilter<OutputImageType, OutputImageType> ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( image );
  thresholder->SetLowerThreshold(
    static_cast<OutputPixelType>( this->GetForegroundValue() ) );
  thresholder->SetUpperThreshold(
    static_cast<OutputPixelType>( this->GetForegroundValue() ) );
  thresholder->SetInsideValue( NumericTraits<OutputPixelType>::One );
  thresholder->SetOutsideValue( NumericTraits<OutputPixelType>::Zero );
  thresholder->Update();

  typedef BinaryBallStructuringElement<OutputPixelType,
    OutputImageDimension> BallStructuringElementType;
  BallStructuringElementType ballStructuringElement;
  ballStructuringElement.SetRadius( this->m_MinimumSizeStructuringElementRadius );
  ballStructuringElement.CreateStructuringElement();

  typedef BinaryErodeImageFilter<OutputImageType, OutputImageType,
    BallStructuringElementType> EroderType;
  typename EroderType::Pointer eroder = EroderType::New();
  eroder->SetInput( thresholder->GetOutput() );
  eroder->SetKernel( ballStructuringElement );
  eroder->SetForegroundValue( NumericTraits<OutputPixelType>::One );
  eroder->SetBackgroundValue( NumericTraits<OutputPixelType>::Zero );
  eroder->Update();

  typedef BinaryBoxStructuringElement<OutputPixelType,
    OutputImageDimension> BoxStructuringElementType;
  BoxStructuringElementType boxStructuringElement;
  boxStructuringElement.SetRadius( 1 );
  boxStructuringElement.CreateStructuringElement();

  typedef itk::BinaryBoundedSpaceDilateImageFilter<
    OutputImageType, OutputImageType, BoxStructuringElementType> ConditionalDilaterType;

  typename ConditionalDilaterType::Pointer dilater = ConditionalDilaterType::New();
  dilater->SetInput( eroder->GetOutput() );
  dilater->SetScaling( 50 );
  dilater->SetKernel( boxStructuringElement );
  dilater->SetBoundedSpaceImage( thresholder->GetOutput() );
  dilater->SetBoundedSpaceValue( NumericTraits<OutputPixelType>::One );
  dilater->SetForegroundValue( NumericTraits<OutputPixelType>::One );
  dilater->Update();

  ImageRegionIterator<OutputImageType> ItD( dilater->GetOutput(),
    dilater->GetOutput()->GetRequestedRegion() );
  ImageRegionIterator<OutputImageType> ItI( image,
    image->GetRequestedRegion() );
  for ( ItD.GoToBegin(), ItI.GoToBegin(); !ItD.IsAtEnd(); ++ItD, ++ItI )
    {
    if ( ItD.Get() == NumericTraits<OutputPixelType>::One )
      {
      ItI.Set( static_cast<OutputPixelType>( this->GetForegroundValue() ) );
      }
    }
}

template <class TInputImage, class TOutput, class TKernel>
void
BinaryReinhardtMorphologicalImageFilter<TInputImage, TOutput, TKernel>
::MaximumDiameterFilter( typename OutputImageType::Pointer image )
{
  typedef BinaryThresholdImageFilter<OutputImageType, OutputImageType> ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( image );
  thresholder->SetLowerThreshold(
    static_cast<OutputPixelType>( this->GetForegroundValue() ) );
  thresholder->SetUpperThreshold(
    static_cast<OutputPixelType>( this->GetForegroundValue() ) );
  thresholder->SetInsideValue( NumericTraits<OutputPixelType>::One );
  thresholder->SetOutsideValue( NumericTraits<OutputPixelType>::Zero );
  thresholder->Update();

  typedef BinaryBallStructuringElement<OutputPixelType,
    OutputImageDimension> BallStructuringElementType;
  BallStructuringElementType ballStructuringElement;
  ballStructuringElement.SetRadius( this->m_MaximumDiameterStructuringElementRadius );
  ballStructuringElement.CreateStructuringElement();

  typedef BinaryMorphologicalOpeningImageFilter<OutputImageType, OutputImageType,
    BallStructuringElementType> OpenerType;
  typename OpenerType::Pointer opener = OpenerType::New();
  opener->SetInput( thresholder->GetOutput() );
  opener->SetKernel( ballStructuringElement );
  opener->SetForegroundValue( NumericTraits<OutputPixelType>::One );
  opener->SetBackgroundValue( NumericTraits<OutputPixelType>::Zero );
  opener->Update();

  typedef BinaryBoxStructuringElement<OutputPixelType,
    OutputImageDimension> BoxStructuringElementType;
  BoxStructuringElementType boxStructuringElement;
  boxStructuringElement.SetRadius( 2 );
  boxStructuringElement.CreateStructuringElement();

  typedef BinaryDilateImageFilter<OutputImageType, OutputImageType,
    BoxStructuringElementType> DilaterType;
  typename DilaterType::Pointer dilater = DilaterType::New();
  dilater->SetInput( opener->GetOutput() );
  dilater->SetKernel( boxStructuringElement );
  dilater->SetForegroundValue( NumericTraits<OutputPixelType>::One );
  dilater->SetBackgroundValue( NumericTraits<OutputPixelType>::Zero );
  dilater->Update();

  ImageRegionIterator<OutputImageType> ItB( dilater->GetOutput(),
    dilater->GetOutput()->GetRequestedRegion() );
  ImageRegionIterator<OutputImageType> ItI( image,
    image->GetRequestedRegion() );
  for ( ItB.GoToBegin(), ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItB, ++ItI )
    {
    if ( ItB.Get() == NumericTraits<OutputPixelType>::One &&
         ItI.Get() == static_cast<OutputPixelType>( this->GetForegroundValue() ) )
      {
      ItB.Set( NumericTraits<OutputPixelType>::One );
      }
    else
      {
      ItB.Set( NumericTraits<OutputPixelType>::Zero );
      }
    }

  for ( ItB.GoToBegin(), ItI.GoToBegin(); !ItB.IsAtEnd(); ++ItB, ++ItI )
    {
    OutputPixelType currentLabel = ItI.Get();
    if ( currentLabel == static_cast<OutputPixelType>( this->GetForegroundValue() ) &&
         ItB.Get() == NumericTraits<OutputPixelType>::Zero )
      {
      ItI.Set( static_cast<OutputPixelType>( this->GetForegroundValue() ) );
      }
    }
}

template <class TInputImage, class TOutput, class TKernel>
void
BinaryReinhardtMorphologicalImageFilter<TInputImage, TOutput, TKernel>
::ConnectivityFilter( typename OutputImageType::Pointer image )
{
  typedef BinaryThresholdImageFilter<OutputImageType, OutputImageType> ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( image );
  thresholder->SetLowerThreshold(
    static_cast<OutputPixelType>( this->GetForegroundValue() ) );
  thresholder->SetUpperThreshold(
    static_cast<OutputPixelType>( this->GetForegroundValue() ) );
  thresholder->SetInsideValue( NumericTraits<OutputPixelType>::One );
  thresholder->SetOutsideValue( NumericTraits<OutputPixelType>::Zero );
  thresholder->Update();

  typedef ConnectedComponentImageFilter<OutputImageType, OutputImageType> ConnectedComponentType;
  typename ConnectedComponentType::Pointer connecter = ConnectedComponentType::New();
  connecter->SetInput( thresholder->GetOutput() );
  connecter->FullyConnectedOff();
  connecter->Update();

  typedef RelabelComponentImageFilter<OutputImageType, OutputImageType> LabelerType;
  typename LabelerType::Pointer labeler = LabelerType::New();
  labeler->SetInput( connecter->GetOutput() );
  labeler->Update();

  ImageRegionIterator<OutputImageType> ItL( labeler->GetOutput(),
    labeler->GetOutput()->GetRequestedRegion() );
  ImageRegionIterator<OutputImageType> ItI( image,
    image->GetRequestedRegion() );

  for ( ItI.GoToBegin(), ItL.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItL )
    {
    OutputPixelType currentLabel = ItI.Get();
    if ( ItL.Get() > 0 && ItL.Get() <= static_cast<OutputPixelType>(
          this->m_NumberOfConnectedComponents ) )
      {
      ItI.Set( static_cast<OutputPixelType>( this->GetForegroundValue() ) );
      }
    else if ( currentLabel == static_cast<OutputPixelType>( this->GetForegroundValue() ) )
      {
      ItI.Set( NumericTraits<OutputPixelType>::Zero );
      }
    }
}

template <class TInputImage, class TOutput, class TKernel>
void
BinaryReinhardtMorphologicalImageFilter<TInputImage, TOutput, TKernel>
::BoundarySmoother( typename OutputImageType::Pointer image )
{
  typedef BinaryThresholdImageFilter<OutputImageType, OutputImageType> ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( image );
  thresholder->SetLowerThreshold(
    static_cast<OutputPixelType>( this->GetForegroundValue() ) );
  thresholder->SetUpperThreshold(
    static_cast<OutputPixelType>( this->GetForegroundValue() ) );
  thresholder->SetInsideValue( NumericTraits<OutputPixelType>::Zero );
  thresholder->SetOutsideValue( NumericTraits<OutputPixelType>::One );
  thresholder->Update();

  typedef ConnectedComponentImageFilter<OutputImageType, OutputImageType> ConnectedComponentType;
  typename ConnectedComponentType::Pointer connecter = ConnectedComponentType::New();
  connecter->SetInput( thresholder->GetOutput() );
  connecter->FullyConnectedOff();
  connecter->Update();

  typedef RelabelComponentImageFilter<OutputImageType, OutputImageType> LabelerType;
  typename LabelerType::Pointer labeler = LabelerType::New();
  labeler->SetInput( connecter->GetOutput() );
  labeler->Update();

  ImageRegionIterator<OutputImageType> ItL( labeler->GetOutput(),
    labeler->GetOutput()->GetRequestedRegion() );
  ImageRegionIterator<OutputImageType> ItI( image,
    image->GetRequestedRegion() );
  ImageRegionIterator<OutputImageType> ItT( thresholder->GetOutput(),
    thresholder->GetOutput()->GetRequestedRegion() );

  ItI.GoToBegin();
  ItL.GoToBegin();
  ItT.GoToBegin();
  while ( !ItI.IsAtEnd() )
    {
    if ( ItL.Get() != 1 && ItI.Get() !=
          static_cast<OutputPixelType>( this->GetForegroundValue() ) )
      {
      ItL.Set( NumericTraits<OutputPixelType>::One );
      }
    else
      {
      ItL.Set( NumericTraits<OutputPixelType>::Zero );
      }
    if ( ItT.Get() == NumericTraits<OutputPixelType>::One )
      {
      ItT.Set( NumericTraits<OutputPixelType>::Zero );
      }
    else
      {
      ItT.Set( NumericTraits<OutputPixelType>::One );
      }
    ++ItI;
    ++ItL;
    ++ItT;
    }

  typedef BinaryBallStructuringElement<OutputPixelType,
    OutputImageDimension> BallStructuringElementType;
  BallStructuringElementType ballStructuringElement;
  ballStructuringElement.SetRadius( this->m_BoundarySmootherStructuringElementRadius );
  ballStructuringElement.CreateStructuringElement();

  typedef BinaryMorphologicalOpeningImageFilter<OutputImageType, OutputImageType,
    BallStructuringElementType> OpenerType;
  typename OpenerType::Pointer opener = OpenerType::New();
  opener->SetInput( thresholder->GetOutput() );
  opener->SetKernel( ballStructuringElement );
  opener->SetForegroundValue( NumericTraits<OutputPixelType>::One );
  opener->SetBackgroundValue( NumericTraits<OutputPixelType>::Zero );
  opener->Update();

  typedef BinaryMorphologicalClosingImageFilter<OutputImageType, OutputImageType,
    BallStructuringElementType> CloserType;
  typename OpenerType::Pointer closer = OpenerType::New();
  closer->SetInput( opener->GetOutput() );
  closer->SetKernel( ballStructuringElement );
  closer->SetForegroundValue( NumericTraits<OutputPixelType>::One );
  closer->SetBackgroundValue( NumericTraits<OutputPixelType>::Zero );
  closer->Update();

  ImageRegionIterator<OutputImageType> ItC( closer->GetOutput(),
    closer->GetOutput()->GetRequestedRegion() );
  ItI.GoToBegin();
  ItC.GoToBegin();
  ItL.GoToBegin();
  while ( !ItC.IsAtEnd() )
    {
    OutputPixelType currentLabel = ItI.Get();
    if ( ItC.Get() == NumericTraits<OutputPixelType>::One &&
         ItL.Get() != NumericTraits<OutputPixelType>::One )
      {
      ItI.Set( static_cast<OutputPixelType>( this->GetForegroundValue() ) );
      }
    else if ( currentLabel == static_cast<OutputPixelType>( this->GetForegroundValue() ) )
      {
      ItI.Set( NumericTraits<OutputPixelType>::Zero );
      }
    ++ItI;
    ++ItC;
    ++ItL;
    }
}

template <class TInputImage, class TOutput, class TKernel>
void
BinaryReinhardtMorphologicalImageFilter<TInputImage, TOutput, TKernel>
::UnclassifiedPixelProcessing( typename OutputImageType::Pointer image )
{
  ImageRegionIterator<OutputImageType> ItI( image,
    image->GetRequestedRegion() );
  for ( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
    {
    if ( ItI.Get() == NumericTraits<OutputPixelType>::Zero )
      {
      ItI.Set( static_cast<OutputPixelType>( this->GetForegroundValue() ) );
      }
    }
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TOutput, class TKernel>
void
BinaryReinhardtMorphologicalImageFilter<TInputImage, TOutput, TKernel>
::PrintSelf( std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  unsigned int step = 1;

  if ( this->m_EmploySaltAndPepperRepair )
    {
    os << indent << "Step " << step++
       << ": EmploySaltAndPepperRepair" << std::endl;
    os << indent << "  SaltAndPepperMinimumSizeInPixels: "
       << this->m_SaltAndPepperMinimumSizeInPixels << std::endl;
    }

  if ( this->m_EmployMinimumDiameterFilter )
    {
    os << indent << "Step " << step++
       << ":  " << "EmployMinimumDiameterFilter" << std::endl;
    os << indent << "  MinimumDiameterStructuringElementRadius: "
       << this->m_MinimumDiameterStructuringElementRadius << std::endl;
    os << indent << "  NumberOfConnectedComponents: "
       << this->m_NumberOfConnectedComponents << std::endl;
    }

  if ( this->m_EmployUnwantedCavityDeletion )
    {
    os << indent << "Step " << step++
       << ":  " << "EmployUnwantedCavityDeletion" << std::endl;
    }

  if ( this->m_EmployMinimumSizeFilter )
    {
    os << indent << "Step " << step++
       << ":  " << "EmployMinimumSizeFilter" << std::endl;
    os << indent << "  MinimumSizeStructuringElementRadius: "
       << this->m_MinimumSizeStructuringElementRadius << std::endl;
    }

  if ( this->m_EmployMaximumDiameterFilter )
    {
    os << indent << "Step " << step++
       << ":  " << "EmployMaximumDiameterFilter" << std::endl;
    os << indent << "  MaximumDiameterStructuringElementRadius: "
       << this->m_MaximumDiameterStructuringElementRadius << std::endl;
    }

  if ( this->m_EmployConnectivityFilter )
    {
    os << indent << "Step " << step++
       << ":  " << "EmployConnectivityFilter" << std::endl;
    os << indent << "  NumberOfConnectedComponents: "
       << this->m_NumberOfConnectedComponents << std::endl;
    }

  if ( this->m_EmployBoundarySmoother )
    {
    os << indent << "Step " << step++
       << ":  " << "EmployBoundarySmoother" << std::endl;
    os << indent << "  BoundarySmootherStructuringElementRadius: "
       << this->m_BoundarySmootherStructuringElementRadius << std::endl;
    }

  if ( this->m_EmployUnclassifiedPixelProcessing )
    {
    os << indent << "Step " << step++
       << ":  " << "EmployUnclassifiedPixelProcessing" << std::endl;
    }

}

} // end namespace itk

#endif

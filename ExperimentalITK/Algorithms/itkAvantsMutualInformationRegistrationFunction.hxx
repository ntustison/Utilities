/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAvantsMutualInformationRegistrationFunction.hxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:13:39 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkAvantsMutualInformationRegistrationFunction_hxx
#define _itkAvantsMutualInformationRegistrationFunction_hxx

#include "itkAvantsMutualInformationRegistrationFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkCovariantVector.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageIterator.h"
#include "vnl/vnl_math.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkBSplineDeformableTransform.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{


/**
 * Constructor
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::AvantsMutualInformationRegistrationFunction()
{
  this->m_NumberOfSpatialSamples = 500;
  this->m_NumberOfHistogramBins = 50;
  this->m_NormalizeGradient = true;
  this->m_InterpolatorIsBSpline = false;
  this->m_TransformIsBSpline = false;
  this->m_OpticalFlow = false;

  // Initialize PDFs to NULL
  this->m_JointPDF = NULL;
  this->m_JointPDFDerivatives = NULL;

  typename TransformType::Pointer transformer 
          = TransformType::New();
  this->SetTransform(transformer);

  typename BSplineInterpolatorType::Pointer interpolator 
          = BSplineInterpolatorType::New();
  this->SetInterpolator( interpolator );

  this->m_FixedImageMask = NULL;
  this->m_MovingImageMask = NULL;

  // Initialize memory
  this->m_MovingImageNormalizedMin = 0.0;
  this->m_FixedImageNormalizedMin = 0.0;
  this->m_MovingImageTrueMin = 0.0;
  this->m_MovingImageTrueMax = 0.0;
  this->m_FixedImageBinSize = 0.0;
  this->m_MovingImageBinSize = 0.0;
  this->m_CubicBSplineDerivativeKernel = NULL;
  this->m_BSplineInterpolator = NULL;
  this->m_DerivativeCalculator = NULL;
  this->m_NumParametersPerDim = 0;
  this->m_NumBSplineWeights = 0;
  this->m_BSplineTransform = NULL;
  this->m_NumberOfParameters = ImageDimension;

  this->m_FixedImageGradientCalculator = GradientCalculatorType::New();
  this->m_MovingImageGradientCalculator = GradientCalculatorType::New();

  typename DefaultInterpolatorType::Pointer interp =  DefaultInterpolatorType::New();
  typename DefaultInterpolatorType::Pointer interp2 = DefaultInterpolatorType::New();

  this->m_MovingImageInterpolator = static_cast<InterpolatorType*>(
    interp.GetPointer() );
  this->m_FixedImageInterpolator = static_cast<InterpolatorType*>(
    interp2.GetPointer() );
  this->m_Interpolator = static_cast<InterpolatorType*>(
    interp.GetPointer() );
}

template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "NumberOfSpatialSamples: ";
  os << this->m_NumberOfSpatialSamples << std::endl;
  os << indent << "NumberOfHistogramBins: ";
  os << this->m_NumberOfHistogramBins << std::endl;

  // Debugging information
  os << indent << "NumberOfParameters: ";
  os << this->m_NumberOfParameters << std::endl;
  os << indent << "FixedImageNormalizedMin: ";
  os << this->m_FixedImageNormalizedMin << std::endl;
  os << indent << "MovingImageNormalizedMin: ";
  os << this->m_MovingImageNormalizedMin << std::endl;
  os << indent << "MovingImageTrueMin: ";
  os << this->m_MovingImageTrueMin << std::endl;
  os << indent << "MovingImageTrueMax: ";
  os << this->m_MovingImageTrueMax << std::endl;
  os << indent << "FixedImageBinSize: "; 
  os << this->m_FixedImageBinSize << std::endl;
  os << indent << "MovingImageBinSize: ";
  os << this->m_MovingImageBinSize << std::endl;
  os << indent << "InterpolatorIsBSpline: ";
  os << this->m_InterpolatorIsBSpline << std::endl;
  os << indent << "TransformIsBSpline: ";
  os << this->m_TransformIsBSpline << std::endl;
}

template <class TFixedImage, class TMovingImage, class TDeformationField> 
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::InitializeIteration()
{
  this->m_CubicBSplineKernel = CubicBSplineFunctionType::New();
  this->m_CubicBSplineDerivativeKernel = CubicBSplineDerivativeFunctionType::New();

  pdfinterpolator = pdfintType::New();
  dpdfinterpolator = dpdfintType::New();
  pdfinterpolator2 = pdfintType2::New();
  pdfinterpolator3 = pdfintType2::New();
  this->m_DerivativeCalculator = DerivativeFunctionType::New();

  this->m_FixedImageGradientCalculator->SetInputImage( this->m_FixedImage );
  this->m_MovingImageGradientCalculator->SetInputImage( this->m_MovingImage );
  this->m_FixedImageInterpolator->SetInputImage( this->m_FixedImage );
  this->m_Interpolator->SetInputImage( this->m_MovingImage );

  this->m_FixedImageSpacing = this->m_FixedImage->GetSpacing();
  this->m_FixedImageOrigin = this->m_FixedImage->GetOrigin();
  this->m_Normalizer = 0.0;
  this->m_NumberOfSpatialSamples = 1;
  for( unsigned int k = 0; k < ImageDimension; k++ )
    {
      this->m_Normalizer += this->m_FixedImageSpacing[k] * this->m_FixedImageSpacing[k];
      this->m_NumberOfSpatialSamples *= this->m_FixedImage->
               GetLargestPossibleRegion().GetSize()[k];
    }
  this->m_Normalizer /= static_cast<double>( ImageDimension );

  this->m_Energy = 0.0;
  
  /**
   * Compute binsize for the histograms.
   *
   * The binsize for the image intensities needs to be adjusted so that 
   * we can avoid dealing with boundary conditions using the cubic 
   * spline as the Parzen window.  We do this by increasing the size
   * of the bins so that the joint histogram becomes "padded" at the 
   * borders. Because we are changing the binsize, 
   * we also need to shift the minimum by the padded amount in order to 
   * avoid minimum values filling in our padded region.
   *
   * Note that there can still be non-zero bin values in the padded region,
   * it's just that these bins will never be a central bin for the Parzen
   * window.
   *
   */

  double movingImageMin = NumericTraits<double>::max();
  double movingImageMax = NumericTraits<double>::NonpositiveMin();
  double fixedImageMin = NumericTraits<double>::max();
  double fixedImageMax = NumericTraits<double>::NonpositiveMin();
 
  typedef ImageRegionConstIterator<MovingImageType> MovingIteratorType;
  MovingIteratorType movingImageIterator(
     this->m_MovingImage, this->m_MovingImage->GetBufferedRegion() );

  for ( movingImageIterator.GoToBegin(); 
        !movingImageIterator.IsAtEnd(); ++movingImageIterator)
    {
    double sample = static_cast<double>( movingImageIterator.Get() );
    double fsample = static_cast<double>( 
        this->m_FixedImage->GetPixel( movingImageIterator.GetIndex() ) );
    if ( sample < movingImageMin )
      {
      movingImageMin = sample;
      }
    if ( sample > movingImageMax )
      {
      movingImageMax = sample;
      }
    if ( fsample < fixedImageMin )
      {
      fixedImageMin = fsample;
      }
    if ( fsample > fixedImageMax )
      {
      fixedImageMax = fsample;
      }  
    }
  
  fixedImageMax = ( fixedImageMax - fixedImageMin ) + fixedImageMin;
  movingImageMax = ( movingImageMax - movingImageMin ) + movingImageMin;
  const int padding = 2;  // this will pad by 2 bins

  this->m_FixedImageBinSize = ( fixedImageMax - fixedImageMin ) /
    static_cast<double>( this->m_NumberOfHistogramBins - 2 * padding );
  this->m_FixedImageNormalizedMin = fixedImageMin / this->m_FixedImageBinSize - 
    static_cast<double>( padding );

  this->m_MovingImageBinSize = ( movingImageMax - movingImageMin ) /
    static_cast<double>( this->m_NumberOfHistogramBins - 2 * padding );
  this->m_MovingImageNormalizedMin = movingImageMin / this->m_MovingImageBinSize -
    static_cast<double>( padding );

  this->m_FixedImageMarginalPDF = MarginalPDFType::New();
  this->m_MovingImageMarginalPDF = MarginalPDFType::New();
  typename MarginalPDFType::RegionType             mPDFRegion;
  typename MarginalPDFType::SizeType               mPDFSize;
  typename MarginalPDFType::IndexType              mPDFIndex;
  mPDFIndex.Fill( 0 ); 
  mPDFSize.Fill( this->m_NumberOfHistogramBins ); 
  mPDFRegion.SetIndex( mPDFIndex );
  mPDFRegion.SetSize( mPDFSize );
  this->m_FixedImageMarginalPDF->SetRegions( mPDFRegion );
  this->m_FixedImageMarginalPDF->Allocate();
  this->m_MovingImageMarginalPDF->SetRegions( mPDFRegion );
  this->m_MovingImageMarginalPDF->Allocate();

  /**
   * Allocate memory for the joint PDF and joint PDF derivatives.
   * The joint PDF and joint PDF derivatives are store as itk::Image.
   */
  this->m_JointPDF = JointPDFType::New();
  this->m_JointPDFDerivatives = JointPDFDerivativesType::New();

  // Instantiate a region, index, size
  JointPDFRegionType            jointPDFRegion;
  JointPDFIndexType              jointPDFIndex;
  JointPDFSizeType              jointPDFSize;

  JointPDFDerivativesRegionType  jointPDFDerivativesRegion;
  JointPDFDerivativesIndexType  jointPDFDerivativesIndex;
  JointPDFDerivativesSizeType    jointPDFDerivativesSize;

  // For the joint PDF define a region starting from {0,0} 
  // with size {this->m_NumberOfHistogramBins, this->m_NumberOfHistogramBins}.
  // The dimension represents fixed image parzen window index
  // and moving image parzen window index, respectively.
  jointPDFIndex.Fill( 0 ); 
  jointPDFSize.Fill( this->m_NumberOfHistogramBins ); 

  jointPDFRegion.SetIndex( jointPDFIndex );
  jointPDFRegion.SetSize( jointPDFSize );

  // Set the regions and allocate
  this->m_JointPDF->SetRegions( jointPDFRegion );
  this->m_JointPDF->Allocate();

  this->m_JointHist = JointPDFType::New();
  this->m_JointHist->SetRegions(jointPDFRegion );
  this->m_JointHist->Allocate();

  // For the derivatives of the joint PDF define a region starting from {0,0,0} 
  // with size {this->m_NumberOfParameters,this->m_NumberOfHistogramBins, 
  // this->m_NumberOfHistogramBins}. The dimension represents transform parameters,
  // fixed image parzen window index and moving image parzen window index,
  // respectively. 
  jointPDFDerivativesIndex.Fill( 0 ); 
  jointPDFDerivativesSize[0] = this->m_NumberOfParameters;
  jointPDFDerivativesSize[1] = this->m_NumberOfHistogramBins;
  jointPDFDerivativesSize[2] = this->m_NumberOfHistogramBins;

  jointPDFDerivativesRegion.SetIndex( jointPDFDerivativesIndex );
  jointPDFDerivativesRegion.SetSize( jointPDFDerivativesSize );

  // Set the regions and allocate
  this->m_JointPDFDerivatives->SetRegions( jointPDFDerivativesRegion );

  this->m_InterpolatorIsBSpline = true;

  BSplineInterpolatorType * testPtr = dynamic_cast<BSplineInterpolatorType *>(
    this->m_Interpolator.GetPointer() );
  if ( !testPtr )
    {
    this->m_InterpolatorIsBSpline = false;

    this->m_DerivativeCalculator = DerivativeFunctionType::New();
    this->m_DerivativeCalculator->SetInputImage( this->m_MovingImage );
    this->m_BSplineInterpolator = NULL;
    } 
  else
    {
    this->m_BSplineInterpolator = testPtr;
    this->m_DerivativeCalculator = NULL;
    }

  this->m_TransformIsBSpline = false;
  this->m_BSplineTransform = NULL;
  
  this->m_NormalizeMetric = 1.0;
  for ( int i = 0; i < ImageDimension; i++ )
    {
    this->m_NormalizeMetric *= this->m_FixedImage->
         GetLargestPossibleRegion().GetSize()[i];
    }
  this->GetProbabilities();
  
  pdfinterpolator->SetInputImage( this->m_JointPDF );
  pdfinterpolator2->SetInputImage( this->m_FixedImageMarginalPDF );
  pdfinterpolator3->SetInputImage( this->m_MovingImageMarginalPDF );
  pdfinterpolator->SetSplineOrder( 3 );
  pdfinterpolator2->SetSplineOrder( 3 );
  pdfinterpolator3->SetSplineOrder( 3 );
  
  this->ComputeMutualInformation();  
}

template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::GetProbabilities() 
{
  for ( unsigned int j = 0; j < this->m_NumberOfHistogramBins; j++ )
    {
      MarginalPDFIndexType mind;
      mind[0] = j;
      this->m_FixedImageMarginalPDF->SetPixel(mind,0);
      this->m_MovingImageMarginalPDF->SetPixel(mind,0);
    }

  // Reset the joint pdfs to zero
  this->m_JointPDF->FillBuffer( 0.0 );
  this->m_JointHist->FillBuffer( 0.0 );

  unsigned long nSamples = 0;
  unsigned long nFixedImageSamples = 0;
  SizeType imagesize = this->m_FixedImage->GetLargestPossibleRegion().GetSize();

  typedef ImageRegionConstIteratorWithIndex<FixedImageType> RandomIterator;
  RandomIterator iter(this->m_FixedImage, 
     this->m_FixedImage->GetLargestPossibleRegion() );  
  for( iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
    {
    FixedImageIndexType index = iter.GetIndex();
    double fixedImageValue = static_cast<double>( iter.Get() );
    CovariantVectorType fixedGradient;
    fixedGradient.Fill( 0 );

    bool inimage = true;
    for ( unsigned int dd = 0; dd < ImageDimension; dd++ )
      {
      if ( index[dd] < 1 || index[dd] >= imagesize[dd]-1 ) 
        {
        inimage = false;
        }
      }    
    if ( inimage ) 
      {
      fixedGradient = this->m_FixedImageGradientCalculator->EvaluateAtIndex( index );
      }
      
    MovingImagePointType mappedPoint;
    double movingImageValue = static_cast<double>( 
               this->m_MovingImage->GetPixel( index ) );
    double movingImageParzenWindowTerm = movingImageValue 
            / static_cast<double>( this->m_MovingImageBinSize )
            - static_cast<double>( this->m_MovingImageNormalizedMin );
    unsigned int movingImageParzenWindowIndex = 
      static_cast<unsigned int>( floor( movingImageParzenWindowTerm ) );

    // Make sure the extreme values are in valid bins
    if ( movingImageParzenWindowIndex < 2 )
      {
      movingImageParzenWindowIndex = 2;
      }
    else if ( movingImageParzenWindowIndex > this->m_NumberOfHistogramBins - 3 )
      {
      movingImageParzenWindowIndex = this->m_NumberOfHistogramBins - 3;
      }

    // Determine parzen window arguments (see eqn 6 of Mattes paper [2]).    
    double fixedImageParzenWindowTerm = fixedImageValue 
         / static_cast<double>( this->m_FixedImageBinSize )
         - static_cast<double>( this->m_FixedImageNormalizedMin );
    unsigned int fixedImageParzenWindowIndex = 
      static_cast<unsigned int>( floor( fixedImageParzenWindowTerm ) );

    // Make sure the extreme values are in valid bins
    if ( fixedImageParzenWindowIndex < 2 )
      {
      fixedImageParzenWindowIndex = 2;
      }
    else if ( fixedImageParzenWindowIndex > this->m_NumberOfHistogramBins - 3 )
      {
      fixedImageParzenWindowIndex = this->m_NumberOfHistogramBins - 3;
      }

    JointPDFValueType *pdfPtr = this->m_JointPDF->GetBufferPointer() +
                  ( fixedImageParzenWindowIndex*this->m_NumberOfHistogramBins );
    int pdfMovingIndex = static_cast<int>( movingImageParzenWindowIndex );
    pdfPtr += pdfMovingIndex;
    *(pdfPtr) += static_cast<PDFValueType>( 1 );
    ++nSamples;
    }

  /**
   * Normalize the PDFs, compute moving image marginal PDF
   *
   */
  typedef ImageRegionIterator<JointPDFType> JointPDFIteratorType;
  JointPDFIteratorType jointPDFIterator ( this->m_JointPDF, 
          this->m_JointPDF->GetBufferedRegion() );

  // Compute joint PDF normalization factor (to ensure joint PDF sum adds to 1.0)
  double jointPDFSum = 0.0;

  jointPDFIterator.GoToBegin();
  while( !jointPDFIterator.IsAtEnd() )
    {
    float temp = jointPDFIterator.Get();
    jointPDFIterator.Set( temp );
    jointPDFSum += temp;
    ++jointPDFIterator;
    }

  if ( jointPDFSum == 0.0 )
    {
    itkExceptionMacro( "Joint PDF summed to zero" );
    }

  // Normalize the PDF bins
  jointPDFIterator.GoToEnd();
  while( !jointPDFIterator.IsAtBegin() )
    {
    --jointPDFIterator;
    jointPDFIterator.Value() /= static_cast<PDFValueType>( jointPDFSum );
    }

  typedef DiscreteGaussianImageFilter<JointPDFType, JointPDFType> dgtype;
  typename dgtype::Pointer dg = dgtype::New();
  dg->SetInput( this->m_JointPDF );
  dg->SetVariance( 1.0 );
  dg->SetUseImageSpacingOff();
  dg->SetMaximumError( .01f );
  dg->Update();
  this->m_JointPDF = dg->GetOutput();

  // Compute moving image marginal PDF by summing over fixed image bins.
  typedef ImageLinearIteratorWithIndex<JointPDFType> JointPDFLinearIterator;
  JointPDFLinearIterator linearIter( 
    this->m_JointPDF, this->m_JointPDF->GetBufferedRegion() );

  linearIter.SetDirection( 0 );
  linearIter.GoToBegin();
  unsigned int fixedIndex = 0;
  
  while( !linearIter.IsAtEnd() )
    {      
    double sum = 0.0;
    while( !linearIter.IsAtEndOfLine() )
      {
      sum += linearIter.Get();
      ++linearIter;
      }
    MarginalPDFIndexType mind;
    mind[0] = fixedIndex;
    this->m_FixedImageMarginalPDF->SetPixel( mind, static_cast<PDFValueType>( sum ) );
    linearIter.NextLine();
    ++fixedIndex;
    }

  linearIter.SetDirection( 1 );
  linearIter.GoToBegin();
  unsigned int movingIndex = 0;

  while( !linearIter.IsAtEnd() )
    {
    double sum = 0.0;
    while( !linearIter.IsAtEndOfLine() )
      {
      sum += linearIter.Get();
      ++linearIter;
      }
    MarginalPDFIndexType mind;
    mind[0] = movingIndex;
    this->m_MovingImageMarginalPDF->SetPixel( mind, static_cast<PDFValueType>( sum ) );

    linearIter.NextLine();
    ++movingIndex;
    }
}

template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::NormalizeJointHist( )
{ 
  typedef ImageRegionIterator<JointPDFType> JointPDFIteratorType;
  JointPDFIteratorType jointPDFIterator ( this->m_JointPDF, this->m_JointPDF->GetBufferedRegion() );


  jointPDFIterator.GoToBegin();
  while( !jointPDFIterator.IsAtEnd() )
    {
    float temp = this->m_JointHist->GetPixel( jointPDFIterator.GetIndex() );
    this->m_JointPDF->SetPixel( jointPDFIterator.GetIndex(), temp );
    ++jointPDFIterator;
    }

  double jointPDFSum = 0.0;
  jointPDFIterator.GoToBegin();
  while( !jointPDFIterator.IsAtEnd() )
    {
      float temp = jointPDFIterator.Get();
      //      if (temp > max) temp=max;
      jointPDFIterator.Set(temp);
      jointPDFSum += temp;
      ++jointPDFIterator;
    }

  jointPDFIterator.GoToEnd();
  while( !jointPDFIterator.IsAtBegin() )
    {
    --jointPDFIterator;
    jointPDFIterator.Value() /= static_cast<PDFValueType>( jointPDFSum );
    }

  typedef DiscreteGaussianImageFilter<JointPDFType, JointPDFType> dgtype;
  typename dgtype::Pointer dg = dgtype::New();
  dg->SetInput( this->m_JointPDF );
  dg->SetVariance( 1.0 );
  dg->SetUseImageSpacingOn();
  dg->SetMaximumError( .01f );
  dg->Update();
  this->m_JointPDF = dg->GetOutput();
 
  // Compute moving image marginal PDF by summing over fixed image bins.
  typedef ImageLinearIteratorWithIndex<JointPDFType> JointPDFLinearIterator;
  JointPDFLinearIterator linearIter( 
    this->m_JointPDF, this->m_JointPDF->GetBufferedRegion() );


  unsigned int fixedIndex = 0;
  linearIter.SetDirection( 0 );
  linearIter.GoToBegin();
  while( !linearIter.IsAtEnd() )
    {      
    double sum = 0.0;
    while ( !linearIter.IsAtEndOfLine() )
     	{
	     sum += linearIter.Get();
	     ++linearIter;
	     } 
    typename MarginalPDFType::IndexType mind;  
    mind[0] = fixedIndex;
    this->m_FixedImageMarginalPDF->GetPixel( mind ) = static_cast<PDFValueType>( sum );
    linearIter.NextLine();
    ++fixedIndex;
    }

  linearIter.SetDirection( 1 );
  linearIter.GoToBegin();
  unsigned int movingIndex = 0;

  while( !linearIter.IsAtEnd() )
    {
    double sum = 0.0;
    while( !linearIter.IsAtEndOfLine() )
     	{
	     sum += linearIter.Get();
	     ++linearIter;
     	}
    typename MarginalPDFType::IndexType mind;  
    mind[0] = movingIndex;
    this->m_MovingImageMarginalPDF->GetPixel( mind ) = static_cast<PDFValueType>( sum );
    
    linearIter.NextLine();
    ++movingIndex;
    }
}

template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::SampleFixedImageDomain( FixedImageSpatialSampleContainer& samples )
{
  // Set up a random interator within the user specified fixed image region.
  typedef ImageRandomConstIteratorWithIndex<FixedImageType> RandomIterator2;
  typedef ImageRegionConstIteratorWithIndex<FixedImageType> RandomIterator;
  RandomIterator randIter( this->m_FixedImage, 
         this->m_FixedImage->GetLargestPossibleRegion() );
  randIter.GoToBegin();

  typename FixedImageSpatialSampleContainer::iterator iter;
  typename FixedImageSpatialSampleContainer::const_iterator end = samples.end();

  for( iter = samples.begin(); iter != end; ++iter )
    {
    // Get sampled index
    FixedImageIndexType index = randIter.GetIndex();
    // Get sampled fixed image value
    (*iter).FixedImageValue = randIter.Get();
    // Translate index to point
    (*iter).FixedImageIndex = index;
    this->m_FixedImage->TransformIndexToPhysicalPoint( index,
                  (*iter).FixedImagePointValue );
    // Jump to random position
    ++randIter;
    }
}

template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::SampleFixedImageDomainLocal( FixedImageSpatialSampleContainer& samples, 
			                            FixedImageIndexType oindex )
{
  SizeType imagesize = this->m_FixedImage->GetLargestPossibleRegion().GetSize();
  SizeType hradius = this->GetRadius();
  hradius.Fill( 2 );

  NeighborhoodIterator<TDeformationField> asamIt( hradius,
          	    this->GetDeformationField(),
	              this->GetDeformationField()->GetLargestPossibleRegion() );
  asamIt.SetLocation( oindex );
  
  int indct = 0;
  typename FixedImageSpatialSampleContainer::iterator iter;
  typename FixedImageSpatialSampleContainer::const_iterator end = samples.end();
  iter = samples.begin();
  while ( indct < m_NumberOfSpatialSamples && indct < asamIt.Size() )
    {
    IndexType index = asamIt.GetIndex(indct);
    bool inimage = true;
    for ( unsigned int dd = 0; dd < ImageDimension; dd++ )
     	{
	     if ( index[dd] < 0 || index[dd] > imagesize[dd]-1 ) 
        {
        inimage = false;
        }
      }
    if ( inimage )
      {
      (*iter).FixedImageValue = this->m_FixedImage->GetPixel( index );
       (*iter).FixedImageIndex = index;
       this->m_FixedImage->TransformIndexToPhysicalPoint( index,
                  (*iter).FixedImagePointValue );	  
       ++iter;
      }
    indct++;
    }
}

template < class TFixedImage, class TMovingImage , class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::ComputeFixedImageParzenWindowIndices( FixedImageSpatialSampleContainer& samples )
{
  typename FixedImageSpatialSampleContainer::iterator iter;
  typename FixedImageSpatialSampleContainer::const_iterator end = samples.end();

  for( iter = samples.begin(); iter != end; ++iter )
    {
    // Determine parzen window arguments (see eqn 6 of Avants paper [2]).  
    double windowTerm = static_cast<double>( (*iter).FixedImageValue ) 
                      / static_cast<double>( this->m_FixedImageBinSize ) 
                      - static_cast<double>( this->m_FixedImageNormalizedMin );
    unsigned int pindex = static_cast<unsigned int>( floor( windowTerm ) );

    // Make sure the extreme values are in valid bins
    if ( pindex < 2 )
      {
      pindex = 2;
      }
    else if ( pindex > this->m_NumberOfHistogramBins - 3 )
      {
      pindex = this->m_NumberOfHistogramBins - 3;
      }
    (*iter).FixedImageParzenWindowIndex = pindex;
    }
}

template <class TFixedImage, class TMovingImage, class TDeformationField>
double
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::GetValueAndDerivative( IndexType oindex,
			                      MeasureType& value,
	                      		DerivativeType& derivative1,
                         DerivativeType& derivative2 ) 
{
  DerivativeType zero( ImageDimension );
  zero.Fill( 0 );

  double fixedImageValue = static_cast<double>( this->m_FixedImage->GetPixel( oindex ) );
  double movingImageValue = static_cast<double>( this->GetMovingImageValue( oindex, zero ) );  
  unsigned int fixedIndex = this->GetFixedValueIndex( fixedImageValue );
  unsigned int movingIndex =this->GetMovingValueIndex( movingImageValue );

  double dJPDF = 0;
  double dFmPDF = 0;
  double dMmPDF = 0;
  double jointPDFValue = 0;
  double fixedImagePDFValue = 0;
  double movingImagePDFValue = 0;

  double movingImageParzenWindowTerm = movingImageValue 
                       / static_cast<double>( this->m_MovingImageBinSize ) 
                       - static_cast<double>( this->m_MovingImageNormalizedMin );
  unsigned int movingImageParzenWindowIndex = 
    static_cast<unsigned int>( floor( movingImageParzenWindowTerm ) );

  double fixedImageParzenWindowTerm = fixedImageValue 
                       / static_cast<double>( this->m_FixedImageBinSize )
                       - static_cast<double>( this->m_FixedImageNormalizedMin );
  unsigned int fixedImageParzenWindowIndex = 
    static_cast<unsigned int>( floor( fixedImageParzenWindowTerm ) );

  if ( movingImageParzenWindowTerm < 2 )
    {
    movingImageParzenWindowTerm = 2;
    }
  else if ( movingImageParzenWindowTerm > this->m_NumberOfHistogramBins - 3 )
    {
    movingImageParzenWindowTerm = this->m_NumberOfHistogramBins - 3;
    }
  if ( fixedImageParzenWindowTerm < 2 )
    {
    fixedImageParzenWindowTerm = 2;
    }
  else if ( fixedImageParzenWindowTerm > this->m_NumberOfHistogramBins - 3 )
    {
    fixedImageParzenWindowTerm = this->m_NumberOfHistogramBins - 3;
    }

  typename JointPDFType::IndexType pdfind2;
  pdfind2[1] = fixedIndex;
  pdfind2[0] = movingIndex;

  typename JointPDFType::PointType pdfind;
  pdfind[1] = fixedImageParzenWindowTerm;
  pdfind[0] = movingImageParzenWindowTerm;
  jointPDFValue=pdfinterpolator->Evaluate(pdfind);
  dJPDF = (1.0)*(pdfinterpolator->EvaluateDerivative( pdfind ))[1];

  typename MarginalPDFType::PointType mind;
  mind[0]=fixedImageParzenWindowTerm;
  dFmPDF = pdfinterpolator2->EvaluateDerivative( mind )[0];
  fixedImagePDFValue = pdfinterpolator2->Evaluate( mind );  
  typename MarginalPDFType::IndexType mind2;
  mind2[0] = fixedIndex;
  mind[0]=movingImageParzenWindowTerm;

  double term1 = 0;
  double term2 = 0;
  if( jointPDFValue > 1e-16 && fixedImagePDFValue > 1e-8 )
    {
    term1 = dJPDF/jointPDFValue;
    term2 = dFmPDF/fixedImagePDFValue;
    value =  vnl_math_abs( term1*(-1.0) + term2);
    if ( value > 1 ) 
      {
      term1 = 0;
      term2 = 0;
      }
    //value =  vnl_math_abs( term1*(-1.0) + term2 );
    }  // end if-block to check non-zero bin contribution
  else 
    {
    value = 0;
    }
  return (term1*(-1.0) + term2);
}

template <class TFixedImage, class TMovingImage, class TDeformationField>
double
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::GetValueAndDerivative2( IndexType oindex,
                          MeasureType& value,
                          DerivativeType& derivative ) 
{
  DerivativeType zero( ImageDimension );
  zero.Fill( 0 );

  double fixedImageValue = static_cast<double>( this->m_FixedImage->GetPixel( oindex ) );
  double movingImageValue0 = static_cast<double>( this->m_MovingImage->GetPixel( oindex ) );
  double movingImageValue = this->GetMovingImageValue( oindex, derivative );
  unsigned int fixedIndex = this->GetFixedValueIndex( fixedImageValue );
  unsigned int movingIndex = this->GetMovingValueIndex( movingImageValue );

  MarginalPDFIndexType mind;
  mind[0]=fixedIndex;
  double fixedImagePDFValue = m_FixedImageMarginalPDF->GetPixel( mind );  
  mind[0] = movingIndex;
  double movingImagePDFValue = m_MovingImageMarginalPDF->GetPixel( mind );

  JointPDFValueType *pdfPtr = m_JointPDF->GetBufferPointer() +
    ( fixedIndex*this->m_NumberOfHistogramBins);
  int pdfMovingIndex = static_cast<int>( movingIndex );
  pdfPtr += pdfMovingIndex;
  double jointPDFValue = *(pdfPtr); 	

  typename JointPDFDerivativesType::IndexType jointPDFDerivIndex;
  double denom =  movingImagePDFValue*fixedImagePDFValue;
  double MIval=0;
  if( jointPDFValue > 1e-16 &&  denom > 1e-16 )
    {
      double entropynorm = 1.0 ;//
      double pRatio =( 1.0 + log( jointPDFValue / denom ));
      MIval = pRatio;
    }
  value = static_cast<MeasureType>( MIval );
  return MIval;
}

template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::GetDerivative( const ParametersType& parameters, DerivativeType & derivative ) const
{
  MeasureType value;
  this->GetValueAndDerivative( parameters, value, derivative );
}

template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::ComputeImageDerivatives( 
  const MovingImagePointType& mappedPoint, 
  ImageDerivativesType& gradient ) const
{
  if( this->m_InterpolatorIsBSpline )
    {
    gradient = this->m_BSplineInterpolator->EvaluateDerivative( mappedPoint );
    }
  else
    {
   gradient = this->m_DerivativeCalculator->Evaluate( mappedPoint );
    }
}

template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::TransformPoint( 
  unsigned int sampleNumber, 
  MovingImagePointType& mappedPoint,
  bool& sampleOk,
  double& movingImageValue ) const
{
  sampleOk = true;
  movingImageValue = 
    this->m_MovingImage->GetPixel( m_FixedImageSamples[sampleNumber].FixedImageIndex );
  this->m_FixedImage->TransformIndexToPhysicalPoint
    ( this->m_FixedImageSamples[sampleNumber].FixedImageIndex, mappedPoint );
  return;
}

// Method to reinitialize the seed of the random number generator
template <class TFixedImage, class TMovingImage, class TDeformationField> void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::ReinitializeSeed()
{
  // This method should be the same used in the ImageRandomIterator
  //  vnl_sample_reseed();
}

// Method to reinitialize the seed of the random number generator
template <class TFixedImage, class TMovingImage, class TDeformationField> void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::ReinitializeSeed(int seed)
{
  // This method should be the same used in the ImageRandomIterator
  //  vnl_sample_reseed(seed);
}

template <class TFixedImage, class TMovingImage, class TDeformationField> void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::GetValueAndDerivative3(IndexType oindex,
			 MeasureType& value,
			 DerivativeType& derivative) 
{
  value = NumericTraits< MeasureType >::Zero;
  derivative.Fill( NumericTraits< MeasureType >::Zero );

  double fixedImageValue = static_cast<double>( this->m_FixedImage->GetPixel( oindex ) );
  double movingImageValue = static_cast<double>( this->m_MovingImage->GetPixel( oindex ) );
  unsigned int fixedIndex = this->GetFixedValueIndex( fixedImageValue );
  unsigned int movingIndex =this->GetMovingValueIndex( movingImageValue );

  MarginalPDFIndexType mind;
  mind[0] = fixedIndex;
  double fixedImagePDFValue = this->m_FixedImageMarginalPDF->GetPixel( mind );  
  mind[0] = movingIndex;
  double movingImagePDFValue = this->m_MovingImageMarginalPDF->GetPixel( mind );

  JointPDFValueType *pdfPtr = this->m_JointPDF->GetBufferPointer() +
    ( fixedIndex*this->m_NumberOfHistogramBins);
  int pdfMovingIndex = static_cast<int>( movingIndex );
  pdfPtr += pdfMovingIndex;
  double jointPDFValue = *(pdfPtr); 	

  double pRatio = 0;
  if( jointPDFValue > 1e-16 &&  fixedImagePDFValue > 1e-16 )
    {
    pRatio = log( jointPDFValue / fixedImagePDFValue );
    for ( unsigned int mu = 0; mu < m_NumberOfParameters; mu++)
     	{
      typename JointPDFDerivativesType::IndexType dpdfind;
      dpdfind[0] = mu;
      dpdfind[1] = fixedIndex;
      dpdfind[2] = movingIndex;
      JointPDFDerivativesType::PixelType val = 
            this->m_JointPDFDerivatives->GetPixel( dpdfind );
      derivative[mu] = val;
     	}
    }  
  value = pRatio;
  return;
}

template <class TFixedImage, class TMovingImage, class TDeformationField> void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::ComputePDFDerivatives( int pdfFixedIndex,
                         int pdfMovingIndex,
                         CovariantVectorType& fixedImageGradientValue,
                         double cubicBSplineDerivativeValue ) const
{
  for ( unsigned int mu = 0; mu < this->m_NumberOfParameters; mu++ )
    {
    typename JointPDFDerivativesType::IndexType dpdfind;
    dpdfind[0]=mu;
    dpdfind[1]=pdfFixedIndex;
    dpdfind[2]=pdfMovingIndex;
    
    JointPDFDerivativesType::PixelType val = this->m_JointPDFDerivatives->GetPixel( dpdfind );
    val -= fixedImageGradientValue[mu] * cubicBSplineDerivativeValue;
    this->m_JointPDFDerivatives->SetPixel( dpdfind,val );
    }
}

template <class TFixedImage, class TMovingImage, class TDeformationField>
void
AvantsMutualInformationRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::GetProbabilitiesLocal( FixedImageIndexType centerindex, float radius ) 
{
  typedef ImageRegionConstIteratorWithIndex<FixedImageType> RandomIterator;
  RandomIterator randIter( this->m_FixedImage, 
       this->m_FixedImage->GetLargestPossibleRegion() );

  for ( unsigned int j = 0; j < this->m_NumberOfHistogramBins; j++ )
    {
    MarginalPDFIndexType mind;
    mind[0] = j;
    this->m_FixedImageMarginalPDF->SetPixel( mind, 0 );
    this->m_MovingImageMarginalPDF->SetPixel( mind, 0 );
    }
  this->m_JointPDF->FillBuffer( 0.0 );
  this->m_JointHist->FillBuffer( 0.0 );

  unsigned long nSamples = 0;
  unsigned long nFixedImageSamples = 0;
  RandomIterator iter( this->m_FixedImage, 
       this->m_FixedImage->GetLargestPossibleRegion() );  

  for( iter.GoToBegin(); !iter.IsAtEnd(); ++iter )
    {
    FixedImageIndexType index = iter.GetIndex();

    float dist = 0;
    for ( unsigned int i = 0; i < ImageDimension; i++) 
      {
      dist += ( index[i]-centerindex[i] )*( index[i]-centerindex[i] );
      }
    dist = sqrt( dist );

    if ( dist < radius )
      {
      double fixedImageValue = iter.Get();
      CovariantVectorType fixedGradient;
      fixedGradient.Fill( 0 );
      SizeType imagesize = this->m_FixedImage->GetLargestPossibleRegion().GetSize();
      
      bool inimage = true;
      for ( unsigned int dd = 0; dd < ImageDimension; dd++ )
       	{
        if ( index[dd] < 1 || index[dd] >= imagesize[dd] - 1  ) 
          {
          inimage = false;
          }
       	}    
      if ( inimage ) 
        {
        fixedGradient = this->m_FixedImageGradientCalculator->EvaluateAtIndex( index );
        }
      MovingImagePointType mappedPoint;
      double movingImageValue;
      movingImageValue = this->m_MovingImage->GetPixel( index );

      double movingImageParzenWindowTerm = movingImageValue 
                           / static_cast<double>( this->m_MovingImageBinSize )
                           - static_cast<double>( this->m_MovingImageNormalizedMin );
      unsigned int movingImageParzenWindowIndex = 
        static_cast<unsigned int>( floor( movingImageParzenWindowTerm ) );

      if ( movingImageParzenWindowIndex < 2 )
        {
        movingImageParzenWindowIndex = 2;
        }
      else if ( movingImageParzenWindowIndex > this->m_NumberOfHistogramBins - 3 )
        {
        movingImageParzenWindowIndex = this->m_NumberOfHistogramBins - 3;
        }

      double fixedImageParzenWindowTerm = fixedImageValue 
                    / static_cast<double>( this->m_FixedImageBinSize )
                    - static_cast<double>( this->m_FixedImageNormalizedMin );
      unsigned int fixedImageParzenWindowIndex = 
        static_cast<unsigned int>( floor( fixedImageParzenWindowTerm ) );

      // Make sure the extreme values are in valid bins
      if ( fixedImageParzenWindowIndex < 2 )
        {
        fixedImageParzenWindowIndex = 2;
        }
      else if ( fixedImageParzenWindowIndex > this->m_NumberOfHistogramBins - 3 )
        {
        fixedImageParzenWindowIndex = this->m_NumberOfHistogramBins - 3;
        }

      typename JointPDFType::IndexType pdfind;

      JointPDFValueType *pdfPtr = this->m_JointPDF->GetBufferPointer() +
        	( fixedImageParzenWindowIndex*this->m_NumberOfHistogramBins );
      // Move the pointer to the first affected bin
      int pdfMovingIndex = static_cast<int>( movingImageParzenWindowIndex );
      pdfPtr += pdfMovingIndex;
      *(pdfPtr) += static_cast<PDFValueType>( 1 );
      ++nSamples;
     	}
    }

  /**
   * Normalize the PDFs, compute moving image marginal PDF
   *
   */
  typedef ImageRegionIterator<JointPDFType> JointPDFIteratorType;
  JointPDFIteratorType jointPDFIterator ( this->m_JointPDF, 
            this->m_JointPDF->GetBufferedRegion() );

  // Compute joint PDF normalization factor (to ensure joint PDF sum adds to 1.0)
  double jointPDFSum = 0.0;

  jointPDFIterator.GoToBegin();
  while( !jointPDFIterator.IsAtEnd() )
    {
    float temp = jointPDFIterator.Get();
    jointPDFIterator.Set(temp);
    jointPDFSum += temp;
    ++jointPDFIterator;
    }

  if ( jointPDFSum == 0.0 )
    {
    itkExceptionMacro( "Joint PDF summed to zero" );
    }

  // Normalize the PDF bins
  jointPDFIterator.GoToEnd();
  while( !jointPDFIterator.IsAtBegin() )
    {
    --jointPDFIterator;
    jointPDFIterator.Value() /= static_cast<PDFValueType>( jointPDFSum );
    }

  typedef DiscreteGaussianImageFilter<JointPDFType,JointPDFType> dgtype;
  typename dgtype::Pointer dg = dgtype::New();
  dg->SetInput( this->m_JointPDF );
  dg->SetVariance( 1.0 );
  dg->SetUseImageSpacingOff();
  dg->SetMaximumError( .01f );
  dg->Update();
  this->m_JointPDF = dg->GetOutput();
 
  // Compute moving image marginal PDF by summing over fixed image bins.
  typedef ImageLinearIteratorWithIndex<JointPDFType> JointPDFLinearIterator;
  JointPDFLinearIterator linearIter( 
    this->m_JointPDF, this->m_JointPDF->GetBufferedRegion() );
  linearIter.SetDirection( 0 );
  linearIter.GoToBegin();
  unsigned int fixedIndex = 0;
  while( !linearIter.IsAtEnd() )
    {      
    double sum = 0.0;
    while ( !linearIter.IsAtEndOfLine() )
	     {
	     sum += linearIter.Get();
      ++linearIter;
	     }
      
    MarginalPDFIndexType mind;
    mind[0] = fixedIndex;
    this->m_FixedImageMarginalPDF->SetPixel( mind, static_cast<PDFValueType>( sum ) );
    linearIter.NextLine();
    ++fixedIndex;
    }

  linearIter.SetDirection( 1 );
  linearIter.GoToBegin();
  unsigned int movingIndex = 0;

  while( !linearIter.IsAtEnd() )
    {
    double sum = 0.0;
    while( !linearIter.IsAtEndOfLine() )
     	{
	     sum += linearIter.Get();
	     ++linearIter;
     	}
    MarginalPDFIndexType mind;
    mind[0]=movingIndex;
    this->m_MovingImageMarginalPDF->SetPixel(mind,static_cast<PDFValueType>(sum));

    
    linearIter.NextLine();
    ++movingIndex;
    }
}

} // end namespace itk


#endif


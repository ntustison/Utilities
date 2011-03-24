/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSectionMutualInformationRegistrationFunction.cxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:13:43 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkSectionMutualInformationRegistrationFunction_txx
#define _itkSectionMutualInformationRegistrationFunction_txx

#include "itkSectionMutualInformationRegistrationFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkCovariantVector.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageIterator.h"
#include "vnl/vnl_math.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkBSplineDeformableTransform.h"
#include "itkExtractImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
//#include "vnl/vnl_sample_reseed.h"
namespace itk
{


/**
 * Constructor
 */
template < class TFixedImage, class TMovingImage , class TDeformationField >
SectionMutualInformationRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::SectionMutualInformationRegistrationFunction()
{

  m_ZeroInZ=false;
  this-> Superclass::m_NormalizeGradient=true;
  this->m_NumberOfSpatialSamples = 500;
  this->m_NumberOfHistogramBins = 50;

  //  this->SetComputeGradient(false); // don't use the default gradient for now

  this->m_InterpolatorIsBSpline = false;
  this->m_TransformIsBSpline    = false;

  // Initialize PDFs to NULL
  m_JointPDF = NULL;
  m_JointPDFDerivatives = NULL;

  m_OpticalFlow=false;
  typename TransformType::Pointer transformer = TransformType::New();
  this->SetTransform(transformer);

  //  typename BSplineInterpolatorType::Pointer interpolator = BSplineInterpolatorType::New();
  InterpolatorPointer interpolator = InterpolatorType::New();
  this->SetInterpolator (interpolator);


  m_FixedImageMask=NULL;
  m_MovingImageMask=NULL;

  // Initialize memory
  m_CubicBSplineDerivativeKernel = NULL;
  m_BSplineInterpolator = NULL;
  m_DerivativeCalculator = NULL;
  m_NumParametersPerDim = 0;
  m_NumBSplineWeights = 0;
  m_BSplineTransform = NULL;
  m_NumberOfParameters = 0;

  m_FixedImageGradientCalculator = GradientCalculatorType::New();
  m_MovingImageGradientCalculator = GradientCalculatorType::New();


  typename DefaultInterpolatorType::Pointer interp =  DefaultInterpolatorType::New();
  typename DefaultInterpolatorType::Pointer interp2 = DefaultInterpolatorType::New();

  m_MovingImageInterpolator = static_cast<InterpolatorType*>(
    interp.GetPointer() );
  m_FixedImageInterpolator = static_cast<InterpolatorType*>(
    interp2.GetPointer() );
  m_Interpolator = static_cast<InterpolatorType*>(
    interp.GetPointer() );

  //std::cout << " done declaring " << std::endl;

}


/**
 * Print out internal information about this class
 */
template < class TFixedImage, class TMovingImage  , class TDeformationField>
void
SectionMutualInformationRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::PrintSelf(std::ostream& os, Indent indent) const
{
  
  Superclass::PrintSelf(os, indent);

  os << indent << "NumberOfSpatialSamples: ";
  os << m_NumberOfSpatialSamples << std::endl;
  os << indent << "NumberOfHistogramBins: ";
  os << m_NumberOfHistogramBins << std::endl;

  // Debugging information
  os << indent << "NumberOfParameters: ";
  os << m_NumberOfParameters << std::endl;
  os << indent << "FixedImageNormalizedMin: ";
  os << m_FixedImageNormalizedMin << std::endl;
  os << indent << "MovingImageNormalizedMin: ";
  os << m_MovingImageNormalizedMin << std::endl;
  os << indent << "FixedImageBinSize: "; 
  os << m_FixedImageBinSize << std::endl;
  os << indent << "MovingImageBinSize: ";
  os << m_MovingImageBinSize << std::endl;
  os << indent << "InterpolatorIsBSpline: ";
  os << m_InterpolatorIsBSpline << std::endl;
  os << indent << "TransformIsBSpline: ";
  os << m_TransformIsBSpline << std::endl;
  
}


/**
 * Initialize
 */
  template <class TFixedImage, class TMovingImage, class TDeformationField> 
void
SectionMutualInformationRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::InitializeIteration()
{

  this->ComputeMetricImage();

  if (this->m_FixedImage )
    m_NumberOfSlices = this->Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[ImageDimension-1];

  m_FixedImageGradientCalculator->SetInputImage( Superclass::m_FixedImage );
  m_MovingImageGradientCalculator->SetInputImage( Superclass::m_MovingImage );
  m_FixedImageInterpolator->SetInputImage( Superclass::m_FixedImage );
  m_Interpolator->SetInputImage( Superclass::m_MovingImage );

  m_FixedImageSpacing    = Superclass::m_FixedImage->GetSpacing();
  m_FixedImageOrigin     = Superclass::m_FixedImage->GetOrigin();
  m_Normalizer      = 0.0;
  unsigned long nsam = 1;
  for( unsigned int k = 0; k < ImageDimension; k++ )
    {
    m_Normalizer += m_FixedImageSpacing[k] * m_FixedImageSpacing[k];
    nsam*=this->Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[k];
    }
  m_Normalizer /= static_cast<double>( ImageDimension );
  m_NumberOfSpatialSamples=(unsigned long)((double)nsam*0.2);
  
  /**
   * Allocate memory for the fixed image sample container.
   */
  if (m_NumberOfSpatialSamples > 20000) m_NumberOfSpatialSamples=20000;
  if (m_NumberOfHistogramBins > 2000) m_NumberOfSpatialSamples=2000;

  m_FixedImageSamples.resize( m_NumberOfSpatialSamples);
    
  for (int i=0; i<m_NumberOfSlices; i++)
    {
      m_FixedImageMarginalPDF[i].resize(m_NumberOfHistogramBins, 0.0 );
      m_MovingImageMarginalPDF[i].resize( m_NumberOfHistogramBins, 0.0 );
    }


  /**
   * Allocate memory for the joint PDF and joint PDF derivatives.
   * The joint PDF and joint PDF derivatives are store as itk::Image.
   */
  m_JointPDF = JointPDFType::New();
  m_JointPDFDerivatives = JointPDFDerivativesType::New();

  // Instantiate a region, index, size
  JointPDFRegionType            jointPDFRegion;
  JointPDFIndexType              jointPDFIndex;
  JointPDFSizeType              jointPDFSize;

  JointPDFDerivativesRegionType  jointPDFDerivativesRegion;
  JointPDFDerivativesIndexType  jointPDFDerivativesIndex;
  JointPDFDerivativesSizeType    jointPDFDerivativesSize;

  // For the joint PDF define a region starting from {0,0} 
  // with size {m_NumberOfHistogramBins, m_NumberOfHistogramBins}.
  // The dimension represents fixed image parzen window index
  // and moving image parzen window index, respectively.
  jointPDFIndex.Fill( 0 ); 
  jointPDFSize.Fill( m_NumberOfHistogramBins ); 

  unsigned int zsize = this->Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[ImageDimension-1];
  jointPDFSize[2]=zsize;

  m_MovingImageNormalizedMin.set_size(zsize);
  m_MovingImageNormalizedMin.fill(1.e16);
  m_FixedImageNormalizedMin.set_size(zsize);
  m_FixedImageNormalizedMin.fill(1.e16);
  m_MovingImageBinSize.set_size(zsize);
  m_MovingImageBinSize.fill(0.1);
  m_FixedImageBinSize.set_size(zsize);
  m_FixedImageBinSize.fill(0.1);

  
  jointPDFRegion.SetIndex( jointPDFIndex );
  jointPDFRegion.SetSize( jointPDFSize );

  // Set the regions and allocate
  m_JointPDF->SetRegions( jointPDFRegion );
  m_JointPDF->Allocate();

  // For the derivatives of the joint PDF define a region starting from {0,0,0} 
  // with size {m_NumberOfParameters,m_NumberOfHistogramBins, 
  // m_NumberOfHistogramBins}. The dimension represents transform parameters,
  // fixed image parzen window index and moving image parzen window index,
  // respectively. 
  jointPDFDerivativesIndex.Fill( 0 ); 
  jointPDFDerivativesSize[0] = m_NumberOfParameters;
  jointPDFDerivativesSize[1] = m_NumberOfHistogramBins;
  jointPDFDerivativesSize[2] = m_NumberOfHistogramBins;

  jointPDFDerivativesRegion.SetIndex( jointPDFDerivativesIndex );
  jointPDFDerivativesRegion.SetSize( jointPDFDerivativesSize );

  //std::cout << " 9 ";
  // Set the regions and allocate
  m_JointPDFDerivatives->SetRegions( jointPDFDerivativesRegion );
  m_JointPDFDerivatives->Allocate();


  /**
   * Setup the kernels used for the Parzen windows.
   */
  m_CubicBSplineKernel = CubicBSplineFunctionType::New();
  m_CubicBSplineDerivativeKernel = CubicBSplineDerivativeFunctionType::New();    


  //std::cout << " 10 ";
  /** 
   * Uniformly sample the fixed image (within the fixed image region)
   * to create the sample points list.
   */
  //std::cout << " 11 ";

  this->SampleFixedImageDomain( m_FixedImageSamples );

  /**
   * Pre-compute the fixed image parzen window index for 
   * each point of the fixed image sample points list.
   */
  //std::cout << " 12 ";
  
  /**
   * Check if the interpolator is of type BSplineInterpolateImageFunction.
   * If so, we can make use of its EvaluateDerivatives method.
   * Otherwise, we instantiate an external central difference
   * derivative calculator.
   *
   * TODO: Also add it the possibility of using the default gradient
   * provided by the superclass.
   *
   */
  m_InterpolatorIsBSpline = true;

  BSplineInterpolatorType * testPtr = dynamic_cast<BSplineInterpolatorType *>(
    this->m_Interpolator.GetPointer() );
  if ( !testPtr )
    {
    m_InterpolatorIsBSpline = false;

    m_DerivativeCalculator = DerivativeFunctionType::New();
    m_DerivativeCalculator->SetInputImage( this->m_MovingImage );

    m_BSplineInterpolator = NULL;
    //    itkDebugMacro( "Interpolator is not BSpline" );
    } 
  else
    {
    m_BSplineInterpolator = testPtr;
    m_DerivativeCalculator = NULL;
    // itkDebugMacro( "Interpolator is BSpline" );
    }

  m_TransformIsBSpline = false;
  m_BSplineTransform = NULL;
  
  m_NormalizeMetric=1.0;
  for (int i=0; i<ImageDimension; i++)
    m_NormalizeMetric*=this->Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[i];
  
  this->GetProbabilities();
  this-> Superclass::m_Energy=0;
  
}



/**
 * Get the both Value and Derivative Measure
 */
template < class TFixedImage, class TMovingImage  , class TDeformationField>
void
SectionMutualInformationRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::GetProbabilities() 
{

  typedef itk::ExtractImageFilter<TFixedImage,ImageSliceType>  extractorType;  
  typedef typename itk::ExtractImageFilter<TFixedImage,ImageSliceType>::Pointer  extractorPointer;
  typedef typename extractorType::InputImageRegionType inRegionType;
  typename ImageSliceType::SpacingType spc2d;
  typename TFixedImage::IndexType index;
  typename TFixedImage::SizeType exsize = Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize();
  typename TFixedImage::SizeType outsize = Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize();
  typename TFixedImage::SpacingType outspacing = Superclass::m_FixedImage->GetSpacing();
  // Reset the joint pdfs to zero
  m_JointPDF->FillBuffer( 0.0 );
  m_JointPDFDerivatives->FillBuffer( 0.0 );

  for (int i=0; i<ImageDimension; i++) 
    {
      index[i]=0;
    }
  for (int i=0; i<ImageDimension; i++) 
    {
      spc2d[i]=outspacing[i];
    }
  exsize[ImageDimension-1]=1;
  
  inRegionType  region;
  region.SetSize(exsize);

  // for each slice  
  unsigned int slc = 0;
  for (slc=0;  slc <  Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[ImageDimension-1]; slc++)
    {
  index[ImageDimension-1]=slc;
  region.SetIndex(index);
	    

  extractorPointer ex1=extractorType::New();      
  ex1->SetInput(this->m_FixedImage);
  ex1->SetExtractionRegion(region);
  ex1->Update();
  m_FixedImageSlice=ex1->GetOutput();     
  m_FixedImageSlice->SetSpacing(spc2d);
  this->ReinitializeSeed();
  this->SampleFixedImageSlice( m_FixedImageSamples , slc);
  this->ComputeFixedImageParzenWindowIndices( m_FixedImageSamples , slc);

  // Reset marginal pdf to all zeros.
  // Assumed the size has already been set to NumberOfHistogramBins
  // in Initialize().
  for ( unsigned int i = 0; i < m_NumberOfSlices; i++ )
  for ( unsigned int j = 0; j < m_NumberOfHistogramBins; j++ )
    {
    m_FixedImageMarginalPDF[i][j]  = 0.0;
    m_MovingImageMarginalPDF[i][j] = 0.0;
    }

  // Declare iterators for iteration over the sample container
  typename FixedImageSpatialSampleContainer::const_iterator fiter;
  typename FixedImageSpatialSampleContainer::const_iterator fend = 
    m_FixedImageSamples.end();

  unsigned long nSamples=0;
  unsigned long nFixedImageSamples=0;

  for ( fiter = m_FixedImageSamples.begin(); fiter != fend; ++fiter )
    {

      bool sampleOk=true;
      double movingImageValue;
      
      typename TFixedImage::IndexType mvind = m_FixedImageSamples[nFixedImageSamples].FixedImageIndex;
      mvind[ImageDimension-1]=slc;
      movingImageValue = this->Superclass::m_MovingImage->GetPixel(mvind);


    if( sampleOk )
      {
	++nSamples; 
      	// Determine parzen window arguments (see eqn 6 of Avants paper [2]).    
	double movingImageParzenWindowTerm =
	  movingImageValue / m_MovingImageBinSize[slc] - m_MovingImageNormalizedMin[slc];
	unsigned int movingImageParzenWindowIndex = 
	  static_cast<unsigned int>( floor( movingImageParzenWindowTerm ) );
	
	// Make sure the extreme values are in valid bins     
	if ( movingImageParzenWindowIndex < 2 )
	  {
	    movingImageParzenWindowIndex = 2;
	  }
	else if ( movingImageParzenWindowIndex > (m_NumberOfHistogramBins - 3) )
	  {
	    movingImageParzenWindowIndex = m_NumberOfHistogramBins - 3;
	  }
	m_FixedImageMarginalPDF[slc][(*fiter).FixedImageParzenWindowIndex] +=
	  static_cast<PDFValueType>( 1 );
	m_MovingImageMarginalPDF[slc][movingImageParzenWindowIndex] +=
	  static_cast<PDFValueType>( 1 );
       
	JointPDFValueType *pdfPtr = m_JointPDF->GetBufferPointer() +
	  ( slc*m_NumberOfHistogramBins*m_NumberOfHistogramBins+
	   (*fiter).FixedImageParzenWindowIndex * m_NumberOfHistogramBins );
	// Move the pointer to the first affected bin
	int pdfMovingIndex = static_cast<int>( movingImageParzenWindowIndex );
	pdfPtr += pdfMovingIndex;
	*(pdfPtr) += static_cast<PDFValueType>( 1 );
	
      } //end if-block check sampleOk

    ++nFixedImageSamples;

    } // end iterating over fixed image spatial sample container for loop

  }
  
  bool smoothjh=true;
  if (smoothjh)
    {
      float sig=1.;
      typedef RecursiveGaussianImageFilter<JointPDFType,JointPDFType>  filterType;      
      {
	typename filterType::Pointer filter = filterType::New();
	filter->SetInput( this->m_JointPDF ); 
	filter->SetDirection( 0 );  
	filter->SetSigma(sig);
	filter->Update();
	this->m_JointPDF =filter->GetOutput();
      }
      {
	typename filterType::Pointer filter = filterType::New();
	filter->SetInput( this->m_JointPDF ); 
	filter->SetDirection( 1 );  
	filter->SetSigma(sig);
	filter->Update();
	this->m_JointPDF =filter->GetOutput();
      }
    }
  
  typedef ImageRegionIterator<JointPDFType> JointPDFIteratorType;
  JointPDFIteratorType jointPDFIterator ( m_JointPDF, m_JointPDF->GetBufferedRegion() );
  jointPDFIterator.GoToBegin();

  for (slc=0;  slc <  Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[ImageDimension-1]; slc++)
    {      
      double jointPDFSum = 0.0; 
      for (int i=0; i<m_NumberOfHistogramBins; i++)
	for (int j=0; j<m_NumberOfHistogramBins; j++)
	  {
	    typename JointPDFType::IndexType ind;
	    ind[0]=i;
	    ind[1]=j;
	    ind[2]=slc;
	    jointPDFSum += m_JointPDF->GetPixel(ind);
	  }      
      //      std::cout << " jpfSum " << jointPDFSum << " slc " << slc <<  std::endl;
      if (jointPDFSum == 0.0) jointPDFSum = 1.0;
      for (int i=0; i<m_NumberOfHistogramBins; i++)
	for (int j=0; j<m_NumberOfHistogramBins; j++)
	  {
	    typename JointPDFType::IndexType ind;
	    ind[0]=i;
	    ind[1]=j;
	    ind[2]=slc;
	    m_JointPDF->SetPixel(ind,m_JointPDF->GetPixel(ind)/jointPDFSum);
	  }      

      for (int i=0; i<m_NumberOfHistogramBins; i++)
	{
	  double sum = 0.0;
	  for (int j=0; j<m_NumberOfHistogramBins; j++)
	    {
	      typename JointPDFType::IndexType ind;
	      ind[0]=i;
	      ind[1]=j;
	      ind[2]=slc;
	      sum+=m_JointPDF->GetPixel(ind);
	    }      
	  m_MovingImageMarginalPDF[slc][i] = static_cast<PDFValueType>(sum);
	}

      for (int i=0; i<m_NumberOfHistogramBins; i++)
	{
	  double sum = 0.0;
	  for (int j=0; j<m_NumberOfHistogramBins; j++)
	    {
	      typename JointPDFType::IndexType ind;
	      ind[0]=j;
	      ind[1]=i;
	      ind[2]=slc;
	      sum+=m_JointPDF->GetPixel(ind);
	    }      
	  m_FixedImageMarginalPDF[slc][i] = static_cast<PDFValueType>(sum);
	}
           

    }


}


/**
 * Uniformly sample the fixed image domain using a random walk
 */
template < class TFixedImage, class TMovingImage , class TDeformationField >
void
SectionMutualInformationRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::SampleFixedImageDomain( FixedImageSpatialSampleContainer& samples )
{
 
  // Set up a random interator within the user specified fixed image region.
  typedef ImageRandomConstIteratorWithIndex<FixedImageType> RandomIterator;
  RandomIterator randIter( this->m_FixedImage,this->GetFixedImage()->GetLargestPossibleRegion() );

  randIter.SetNumberOfSamples( m_NumberOfSpatialSamples );
  randIter.GoToBegin();

  typename FixedImageSpatialSampleContainer::iterator iter;
  typename FixedImageSpatialSampleContainer::const_iterator end=samples.end();

  if( this->m_FixedImageMask )
    {

    }
  else
    {
    for( iter=samples.begin(); iter != end; ++iter )
      {
      // Get sampled index
      FixedImageIndexType index = randIter.GetIndex();
      // Get sampled fixed image value
      (*iter).FixedImageValue = randIter.Get();
      // Translate index to point
      (*iter).FixedImageIndex = index;
      this->Superclass::m_FixedImage->TransformIndexToPhysicalPoint( index,
                                                   (*iter).FixedImagePointValue );
      // Jump to random position
      ++randIter;

      }
    }
}


/**
 * Uniformly sample the fixed image domain using a random walk
 */
template < class TFixedImage, class TMovingImage , class TDeformationField >
void
SectionMutualInformationRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::SampleFixedImageSlice( FixedImageSpatialSampleContainer& samples , unsigned int slc)
{
 
  // Set up a random interator within the user specified fixed image region.
  typedef ImageRandomConstIteratorWithIndex<ImageSliceType> RandomIterator;
  RandomIterator randIter( this->m_FixedImageSlice,this->m_FixedImageSlice->GetLargestPossibleRegion() );

  randIter.SetNumberOfSamples( m_NumberOfSpatialSamples );
  randIter.GoToBegin();

  typename FixedImageSpatialSampleContainer::iterator iter;
  typename FixedImageSpatialSampleContainer::const_iterator end=samples.end();

  double movingImageMin = NumericTraits<double>::max();
  double movingImageMax = NumericTraits<double>::NonpositiveMin();
  double fixedImageMin = NumericTraits<double>::max();
  double fixedImageMax = NumericTraits<double>::NonpositiveMin();

  for( iter=samples.begin(); iter != end; ++iter )
    {
      // Get sampled index
      SliceIndexType index = randIter.GetIndex();
      // Get sampled fixed image value
      (*iter).FixedImageValue = randIter.Get();
      // Translate index to point
      (*iter).FixedImageIndex = index;
      this->m_FixedImageSlice->TransformIndexToPhysicalPoint( index,
			      (*iter).FixedImagePointValue );

      double sample = static_cast<double>(  randIter.Get() );
      
      if ( sample < fixedImageMin )
	{
	  fixedImageMin = sample;
	}
      
      if ( sample > fixedImageMax )
	{
	  fixedImageMax = sample;
	}
      
      sample = static_cast<double>( this->Superclass::m_MovingImage->GetPixel(index) );
      
      if ( sample < movingImageMin )
	{
	  movingImageMin = sample;
	}
      
      if ( sample > movingImageMax )
	{
	  movingImageMax = sample;
	}
      // Jump to random position
      ++randIter;
      
    }
  
  const int padding = 2;  // this will pad by 2 bins

  m_FixedImageBinSize[slc] = ( fixedImageMax - fixedImageMin ) /
    static_cast<double>( m_NumberOfHistogramBins - 2 * padding );
  m_FixedImageNormalizedMin[slc] = fixedImageMin / m_FixedImageBinSize[slc] - 
    static_cast<double>( padding );

  m_MovingImageBinSize[slc] = ( movingImageMax - movingImageMin ) /
    static_cast<double>( m_NumberOfHistogramBins - 2 * padding );
  m_MovingImageNormalizedMin[slc] = movingImageMin / m_MovingImageBinSize[slc] -
    static_cast<double>( padding );

}




/**
 * Uniformly sample the fixed image domain using a random walk
 */
template < class TFixedImage, class TMovingImage , class TDeformationField >
void
SectionMutualInformationRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::ComputeFixedImageParzenWindowIndices( FixedImageSpatialSampleContainer& samples , unsigned int slc)
{

  typename FixedImageSpatialSampleContainer::iterator iter;
  typename FixedImageSpatialSampleContainer::const_iterator end=samples.end();

  for( iter=samples.begin(); iter != end; ++iter )
    {

    // Determine parzen window arguments (see eqn 6 of Avants paper [2]).  
    double windowTerm =
      static_cast<double>( (*iter).FixedImageValue ) / m_FixedImageBinSize[slc] -
        m_FixedImageNormalizedMin[slc];
    unsigned int pindex = static_cast<unsigned int>( floor( windowTerm ) );

    // Make sure the extreme values are in valid bins
    if ( pindex < 2 )
      {
      pindex = 2;
      }
    else if ( pindex > (m_NumberOfHistogramBins - 3) )
      {
      pindex = m_NumberOfHistogramBins - 3;
      }

    (*iter).FixedImageParzenWindowIndex = pindex;

    }

}


/**
 * Get the both Value and Derivative Measure
 */
template < class TFixedImage, class TMovingImage  , class TDeformationField>
void
SectionMutualInformationRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::GetValueAndDerivative(IndexType oindex,
			MeasureType& value,
			DerivativeType& derivative) 
{
  
  DerivativeType zero(ImageDimension);
  zero.Fill(0);

  int sign = 0;
  // if (derivative[ImageDimension-1] > 0.5) sign = 1;
  //  if (derivative[ImageDimension-1] < -0.5) sign = -1;
  unsigned int slc = oindex[ImageDimension-1]+sign;

  double fixedImageValue = (double)this->Superclass::m_FixedImage->GetPixel(oindex);
  double movingImageValue = this->GetMovingImageValue(oindex,derivative);  
  unsigned int fixedIndex = this->GetFixedValueIndex(fixedImageValue,slc);
  unsigned int movingIndex =this->GetMovingValueIndex(movingImageValue,slc);

  double fixedImagePDFValue = m_FixedImageMarginalPDF[slc][fixedIndex];  
  double movingImagePDFValue = m_MovingImageMarginalPDF[slc][movingIndex];

  JointPDFValueType *pdfPtr = m_JointPDF->GetBufferPointer() +
    ( slc*m_NumberOfHistogramBins*m_NumberOfHistogramBins+fixedIndex* m_NumberOfHistogramBins);
  int pdfMovingIndex = static_cast<int>( movingIndex );
  pdfPtr += pdfMovingIndex;
  double jointPDFValue = *(pdfPtr); 	
	  
  // check for non-zero bin contribution
  typename JointPDFDerivativesType::IndexType jointPDFDerivIndex;
  double denom =  movingImagePDFValue*fixedImagePDFValue;
  double derivMIval=0;
  if( jointPDFValue > 1e-16 &&  denom > 1e-16)
    {
      double entropynorm = 1.0 ;//
      double pRatio =( 1.0 + log( jointPDFValue / denom ));
      derivMIval = pRatio ;
    }  // end if-block to check non-zero bin contribution
  
  
  value = static_cast<MeasureType>( derivMIval );

}


/**
 * Get the match measure derivative
 */
template < class TFixedImage, class TMovingImage  , class TDeformationField>
void
SectionMutualInformationRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::GetDerivative( const ParametersType& parameters, DerivativeType & derivative ) const
{
  MeasureType value;
  // call the combined version
  this->GetValueAndDerivative( parameters, value, derivative );
}


/**
 * Compute image derivatives using a central difference function
 * if we are not using a BSplineInterpolator, which includes
 * derivatives.
 */
template < class TFixedImage, class TMovingImage , class TDeformationField >
void
SectionMutualInformationRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::ComputeImageDerivatives( 
  const MovingImagePointType& mappedPoint, 
  ImageDerivativesType& gradient ) const
{
  
  if( m_InterpolatorIsBSpline )
    {
    // Computed moving image gradient using derivative BSpline kernel.
    gradient = m_BSplineInterpolator->EvaluateDerivative( mappedPoint );
    }
  else
    {
    // For all generic interpolator use central differencing.
   gradient = m_DerivativeCalculator->Evaluate( mappedPoint );
    }

}


/**
 * Transform a point from FixedImage domain to MovingImage domain.
 * This function also checks if mapped point is within support region. 
 */
template < class TFixedImage, class TMovingImage , class TDeformationField >
void
SectionMutualInformationRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::TransformPoint( 
  unsigned int sampleNumber, 
  MovingImagePointType& mappedPoint,
  bool& sampleOk,
  double& movingImageValue ) const
{
  
  sampleOk=true;
  movingImageValue = 
    this->Superclass::m_MovingImage->GetPixel(  m_FixedImageSamples[sampleNumber].FixedImageIndex);
  this->Superclass::m_FixedImage->TransformIndexToPhysicalPoint
    (m_FixedImageSamples[sampleNumber].FixedImageIndex,mappedPoint);
  return;
  
}


/**
 * Compute PDF derivatives contribution for each parameter
 */
template < class TFixedImage, class TMovingImage , class TDeformationField >
void
SectionMutualInformationRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::ComputePDFDerivatives( 
  unsigned int sampleNumber, 
  int pdfMovingIndex,
  const ImageDerivativesType& movingImageGradientValue,
  double cubicBSplineDerivativeValue ) const
{


  typename JointPDFDerivativesType::IndexType jointPDFDerivIndex;
  jointPDFDerivIndex[2]=pdfMovingIndex;
  jointPDFDerivIndex[1]= m_FixedImageSamples[sampleNumber].FixedImageParzenWindowIndex;
  
  /**
   * Generic version which works for all transforms.
   */
  // Compute the transform Jacobian.
  typedef typename TransformType::JacobianType JacobianType;
  const JacobianType& jacobian =  this->m_Transform->GetJacobian
    ( m_FixedImageSamples[sampleNumber].FixedImagePointValue );
  
  
  for ( unsigned int mu = 0; mu < m_NumberOfParameters; mu++)//, derivPtr++ )
    {
      double innerProduct = 0.0;
      for ( unsigned int dim = 0; dim < ImageDimension; dim++ )
        {
	  innerProduct += jacobian[dim][mu] * 
	    movingImageGradientValue[mu];
	}
      jointPDFDerivIndex[0]=mu;
      float val=m_JointPDFDerivatives->GetPixel(jointPDFDerivIndex);
      val -= innerProduct * cubicBSplineDerivativeValue;
      m_JointPDFDerivatives->SetPixel(jointPDFDerivIndex,val);
      
    }
  

}


// Method to reinitialize the seed of the random number generator
template < class TFixedImage, class TMovingImage  , class TDeformationField> void
SectionMutualInformationRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::ReinitializeSeed()
{
  // This method should be the same used in the ImageRandomIterator
  //vnl_sample_reseed();
}

// Method to reinitialize the seed of the random number generator
template < class TFixedImage, class TMovingImage  , class TDeformationField> void
SectionMutualInformationRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::ReinitializeSeed(int seed)
{
  // This method should be the same used in the ImageRandomIterator
 // vnl_sample_reseed(seed);
}

//






} // end namespace itk


#endif


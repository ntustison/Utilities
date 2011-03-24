/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSectionMutualInformationRegistrationFunction.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:13:44 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkSectionMutualInformationRegistrationFunction_h
#define __itkSectionMutualInformationRegistrationFunction_h

#include "itkImageToImageMetric.h"
#include "itkAvantsPDEDeformableRegistrationFunction.h"
#include "itkCovariantVector.h"
#include "itkPoint.h"
#include "itkIndex.h"
#include "itkBSplineKernelFunction.h"
#include "itkBSplineDerivativeKernelFunction.h"
#include "itkCentralDifferenceImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkBSplineDeformableTransform.h"
#include "itkTranslationTransform.h"
#include "itkArray2D.h"
#include "itkImageBase.h"
#include "itkTransform.h"
#include "itkInterpolateImageFunction.h"
#include "itkSingleValuedCostFunction.h"
#include "itkExceptionObject.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkSpatialObject.h"
#include "itkConstNeighborhoodIterator.h"

namespace itk
{

/** \class SectionMutualInformationRegistrationFunction
 * \brief Computes the mutual information between two images to be 
 * registered using the method of Avants et al.
 *
 * SectionMutualInformationRegistrationFunction computes the mutual 
 * information between a fixed and moving image to be registered.
 *
 * This class is templated over the FixedImage type and the MovingImage 
 * type.
 *
 * The fixed and moving images are set via methods SetFixedImage() and
 * SetMovingImage(). This metric makes use of user specified Transform and
 * Interpolator. The Transform is used to map points from the fixed image to
 * the moving image domain. The Interpolator is used to evaluate the image
 * intensity at user specified geometric points in the moving image.
 * The Transform and Interpolator are set via methods SetTransform() and
 * SetInterpolator().
 *
 * If a BSplineInterpolationFunction is used, this class obtain
 * image derivatives from the BSpline interpolator. Otherwise, 
 * image derivatives are computed using central differencing.
 *
 * \warning This metric assumes that the moving image has already been
 * connected to the interpolator outside of this class. 
 *
 * The method GetValue() computes of the mutual information
 * while method GetValueAndDerivative() computes
 * both the mutual information and its derivatives with respect to the
 * transform parameters.
 *
 * The calculations are based on the method of Avants et al [1,2]
 * where the probability density distribution are estimated using
 * Parzen histograms. Since the fixed image PDF does not contribute
 * to the derivatives, it does not need to be smooth. Hence, 
 * a zero order (box car) BSpline kernel is used
 * for the fixed image intensity PDF. On the other hand, to ensure
 * smoothness a third order BSpline kernel is used for the 
 * moving image intensity PDF.
 *
 * On Initialize(), the FixedImage is uniformly sampled within
 * the FixedImageRegion. The number of samples used can be set
 * via SetNumberOfSpatialSamples(). Typically, the number of
 * spatial samples used should increase with the image size.
 *
 * During each call of GetValue(), GetDerivatives(),
 * GetValueAndDerivatives(), marginal and joint intensity PDF's
 * values are estimated at discrete position or bins. 
 * The number of bins used can be set via SetNumberOfHistogramBins().
 * To handle data with arbitray magnitude and dynamic range, 
 * the image intensity is scale such that any contribution to the
 * histogram will fall into a valid bin.
 *
 * One the PDF's have been contructed, the mutual information
 * is obtained by doubling summing over the discrete PDF values.
 *
 *
 * Notes: 
 * 1. This class returns the negative mutual information value.
 * 2. This class in not thread safe due the private data structures
 *     used to the store the sampled points and the marginal and joint pdfs.
 *
 * References:
 * [1] "Nonrigid multimodality image registration"
 *      D. Avants, D. R. Haynor, H. Vesselle, T. Lewellen and W. Eubank
 *      Medical Imaging 2001: Image Processing, 2001, pp. 1609-1620.
 * [2] "PET-CT Image Registration in the Chest Using Free-form Deformations"
 *      D. Avants, D. R. Haynor, H. Vesselle, T. Lewellen and W. Eubank
 *      IEEE Transactions in Medical Imaging. Vol.22, No.1, 
        January 2003. pp.120-128.
 * [3] "Optimization of Mutual Information for MultiResolution Image
 *      Registration"
 *      P. Thevenaz and M. Unser
 *      IEEE Transactions in Image Processing, 9(12) December 2000.
 *
 * \ingroup RegistrationMetrics
 * \ingroup ThreadUnSafe
 */
  template <class TFixedImage,class TMovingImage , class TDeformationField>
class ITK_EXPORT SectionMutualInformationRegistrationFunction :
  public AvantsPDEDeformableRegistrationFunction< TFixedImage, TMovingImage , TDeformationField>
{
public:
  /** Standard class typedefs. */
  typedef SectionMutualInformationRegistrationFunction    Self;
  typedef AvantsPDEDeformableRegistrationFunction< TFixedImage,
    TMovingImage, TDeformationField >    Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( SectionMutualInformationRegistrationFunction, 
    AvantsPDEDeformableRegistrationFunction );

  /** Inherit some enums from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

  /** MovingImage image type. */
  typedef typename Superclass::MovingImageType     MovingImageType;
  typedef typename Superclass::MovingImagePointer  MovingImagePointer;

  /** FixedImage image type. */
  typedef typename Superclass::FixedImageType     FixedImageType;
  typedef typename Superclass::FixedImagePointer  FixedImagePointer;
  typedef typename FixedImageType::IndexType      IndexType;
  typedef typename FixedImageType::SizeType       SizeType;
  typedef typename FixedImageType::SpacingType    SpacingType;
  
  typedef Image<typename TFixedImage::PixelType,ImageDimension> ImageSliceType;
  typedef typename ImageSliceType::IndexType      SliceIndexType;
  typedef typename ImageSliceType::PointType      SlicePointType;

  /** Deformation field type. */
  typedef typename Superclass::VectorType VectorType;
  typedef typename Superclass::DeformationFieldType    DeformationFieldType;
  typedef typename Superclass::DeformationFieldTypePointer   
    DeformationFieldTypePointer;


  /** Inherit some enums from the superclass. */
  typedef typename Superclass::PixelType     PixelType;
  typedef typename Superclass::RadiusType    RadiusType;
  typedef typename Superclass::NeighborhoodType    NeighborhoodType;
  typedef typename Superclass::FloatOffsetType  FloatOffsetType;
  typedef typename Superclass::TimeStepType TimeStepType;


  /** Interpolator type. */
  typedef double CoordRepType;
  //typedef NearestNeighborInterpolateImageFunction<MovingImageType,CoordRepType>
  typedef LinearInterpolateImageFunction<MovingImageType,CoordRepType>
    ///    BSplineInterpolateImageFunction<MovingImageType,CoordRepType> 
    InterpolatorType;
  typedef typename InterpolatorType::Pointer         InterpolatorPointer;
  typedef typename InterpolatorType::PointType       PointType;
  typedef InterpolatorType DefaultInterpolatorType;
  //  typedef LinearInterpolateImageFunction<MovingImageType,CoordRepType>
  //DefaultInterpolatorType;

  /** Covariant vector type. */
  typedef CovariantVector<double,itkGetStaticConstMacro(ImageDimension)> CovariantVectorType;

  /** Gradient calculator type. */
  typedef CentralDifferenceImageFunction<FixedImageType> GradientCalculatorType;
  typedef typename GradientCalculatorType::Pointer   GradientCalculatorPointer;

  /** Set the moving image interpolator. */
  void SetMovingImageInterpolator( InterpolatorType * ptr )
  { m_MovingImageInterpolator = ptr; }
  
  /** Get the moving image interpolator. */
  InterpolatorType * GetMovingImageInterpolator(void)
    { return m_MovingImageInterpolator; }
  
  /** This class uses a constant timestep of 1. */
  virtual TimeStepType ComputeGlobalTimeStep(void *itkNotUsed(GlobalData)) const
    { return 1; }

  /** Return a pointer to a global data structure that is passed to
   * this object from the solver at each calculation.  */
  virtual void *GetGlobalDataPointer() const
  {
        GlobalDataStruct *global = new GlobalDataStruct();
    //    global->m_SumOfSquaredDifference  = 0.0;
    /// global->m_NumberOfPixelsProcessed = 0L;
    // global->m_SumOfSquaredChange      = 0;
       return global;
  }

  
  /** Release memory for global data structure. */
  virtual void ReleaseGlobalDataPointer( void *GlobalData ) const
  { 
     delete (GlobalDataStruct *) GlobalData;  
  }

  /** Set the object's state before each iteration. */
  virtual void InitializeIteration();


  typedef double CoordinateRepresentationType;

  /** Types inherited from Superclass. */
  typedef TranslationTransform<CoordinateRepresentationType, 
    //                    itkGetStaticConstMacro(ImageDimension),
                    itkGetStaticConstMacro(ImageDimension)> TransformType;

  typedef ImageToImageMetric< TFixedImage, TMovingImage > Metricclass;

  typedef typename TransformType::Pointer            TransformPointer;
  typedef typename Metricclass::TransformJacobianType    TransformJacobianType;
  //  typedef typename Metricclass::InterpolatorType         InterpolatorType;
  typedef typename Metricclass::MeasureType              MeasureType;
  typedef typename Metricclass::DerivativeType           DerivativeType;
  typedef typename TransformType::ParametersType           ParametersType;
  typedef typename Metricclass::FixedImageConstPointer   FixedImageConstPointer;
  typedef typename Metricclass::MovingImageConstPointer  MovingImageCosntPointer;
  // typedef typename TransformType::CoordinateRepresentationType  CoordinateRepresentationType;

  /** Index and Point typedef support. */
  typedef typename ImageSliceType::IndexType            FixedImageIndexType;
  typedef typename ImageSliceType::IndexValueType  FixedImageIndexValueType;
  typedef typename MovingImageType::IndexType           MovingImageIndexType;
  typedef typename TransformType::InputPointType        FixedImagePointType;
  typedef typename TransformType::OutputPointType       MovingImagePointType;



  /** Get the derivatives of the match measure. */
  void GetDerivative( 
    const ParametersType& parameters,
    DerivativeType & Derivative ) const;


  /**  Get the value and derivatives for single valued optimizers. */
  void GetValueAndDerivative( IndexType index, 
                              MeasureType& Value, DerivativeType& Derivative ) ;

  /** Number of spatial samples to used to compute metric */
  //  itkSetClampMacro( NumberOfSpatialSamples, unsigned long,
  //                1, NumericTraits<unsigned long>::max() );
  //itkGetConstReferenceMacro( NumberOfSpatialSamples, unsigned long); 

  /** Number of bins to used in the histogram. Typical value is 50. */
  //  itkSetClampMacro( NumberOfHistogramBins, unsigned long,
  //                1, NumericTraits<unsigned long>::max() );
  // itkGetConstReferenceMacro( NumberOfHistogramBins, unsigned long);   
  void SetNumberOfHistogramBins(unsigned long nhb) { m_NumberOfHistogramBins=nhb;}
  unsigned long GetNumberOfHistogramBins() {return m_NumberOfHistogramBins;}
  void SetNumberOfSpatialSamples(unsigned long nhb) { m_NumberOfSpatialSamples=nhb;}
  unsigned long GetNumberOfSpatialSamples() {return m_NumberOfSpatialSamples;}


  /** Provide API to reinitialize the seed of the random number generator */
  static void ReinitializeSeed();
  static void ReinitializeSeed(int);  


  void SetTransform(TransformPointer t){m_Transform=t;}
  TransformPointer GetTransform(){return m_Transform;}
  void SetInterpolator(InterpolatorPointer t){m_Interpolator=t;}
  InterpolatorPointer GetInterpolator(){ return m_Interpolator; }

  void GetProbabilities();

  virtual CovariantVectorType  OpticalFlowUpdate(const NeighborhoodType &neighborhood)
  {
    // Get fixed image related information
    IndexType index=neighborhood.GetIndex();    
    typename TDeformationField::PixelType vec;
    if ( Superclass::m_DeformationField )
      {
      vec = Superclass::m_DeformationField->GetPixel(index);
      }
    else
      {
      vec.Fill( 0 );
      } 
    CovariantVectorType update;
    double fixedValue;
    CovariantVectorType fixedGradient;
    double fixedGradientSquaredMagnitude = 0;
    fixedValue = (double) Superclass::m_FixedImage->GetPixel( index );
    fixedGradient = m_FixedImageGradientCalculator->EvaluateAtIndex( index );
    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
	fixedGradientSquaredMagnitude += vnl_math_sqr( fixedGradient[j] );
      } 
    double movingValue;
    int j;
    PointType mappedPoint;
    for( j = 0; j < ImageDimension; j++ )
      {
	mappedPoint[j] = double( index[j] ) * m_FixedImageSpacing[j] + 
	  m_FixedImageOrigin[j];
	mappedPoint[j] += vec[j];
      }
    if( m_MovingImageInterpolator->IsInsideBuffer( mappedPoint ) )
      {
	movingValue = m_MovingImageInterpolator->Evaluate( mappedPoint );
      }
    else
      {
	for( j = 0; j < ImageDimension; j++ )
	  {
	    update[j] = 0.0;
	  }
	return update;
      }
    double speedValue = fixedValue - movingValue;
    double denominator = vnl_math_sqr( speedValue ) / m_Normalizer + 
      fixedGradientSquaredMagnitude;
    double m_DenominatorThreshold = 1e-9;
    double m_IntensityDifferenceThreshold = 0.001;  
    if ( vnl_math_abs(speedValue) < m_IntensityDifferenceThreshold || 
	 denominator < m_DenominatorThreshold )
      {
	for( j = 0; j < ImageDimension; j++ )
	  {
	    update[j] = 0.0;
	  }
	return update;
      }
    for( j = 0; j < ImageDimension; j++ )
      {
	update[j] = speedValue * fixedGradient[j] / denominator;
      }
    return update;
  }


  virtual VectorType ComputeUpdate(const NeighborhoodType &neighborhood,
                                   void *globalData,
                                   const FloatOffsetType &offset = FloatOffsetType(0.0))
  {
    VectorType update;
    update.Fill(0.0);   
    IndexType oindex = neighborhood.GetIndex();
    
    FixedImageType* img =const_cast<FixedImageType *>(this->Superclass::m_FixedImage.GetPointer()); 
    if (!img) return update;
    typename FixedImageType::SpacingType spacing=img->GetSpacing();
    typename FixedImageType::SizeType imagesize=img->GetLargestPossibleRegion().GetSize();
    bool inimage=true;
    
    for (unsigned int dd=0; dd<ImageDimension; dd++)
      {
	if ( oindex[dd] < 1 || 
	     oindex[dd] >= static_cast<typename IndexType::IndexValueType>(imagesize[dd]-1) ) 
	  return update;
      }    
    
    CovariantVectorType fixedGradient;
    //    ImageDerivativesType fixedGradient;
    CovariantVectorType fixedGradientNorm;
    //std::cout << " grad " << std::endl;

    double loce=0.0;
    double nccp1=0,nccm1=0;
    ParametersType fdvec1(ImageDimension);
    ParametersType fdvec2(ImageDimension);    
    
    if (m_OpticalFlow) 
      {
	fixedGradientNorm = this->OpticalFlowUpdate(neighborhood);		
	for (int imd=0; imd<ImageDimension; imd++)
	  {
	    fdvec1[imd]=fixedGradientNorm[imd]*spacing[imd]*0.6;
	    fdvec2[imd]=fixedGradientNorm[imd]*(-1.)*spacing[imd]*0.6;
	  }
	if (m_ZeroInZ) { fdvec1[ImageDimension-1]=0; fdvec2[ImageDimension-1]=0; }
	this->GetValueAndDerivative(oindex,nccp1,fdvec1);
	this->GetValueAndDerivative(oindex,nccm1,fdvec2);
	loce+=(nccp1+nccm1);
	this-> Superclass::m_Energy+=loce/m_NormalizeMetric;
	float sign=nccp1-nccm1;
	if (sign < 0 ) sign=(0.0);  
	if (sign > 0.0 ) sign=1.0;
	for (int imd=0; imd<ImageDimension; imd++) update[imd]=sign*fixedGradientNorm[imd];
      }
    else
      {
	fixedGradient = m_FixedImageGradientCalculator->EvaluateAtIndex( oindex ); 
	//      fixedGradient = m_MovingImageGradientCalculator->EvaluateAtIndex( oindex ); 
	//PointType pt;  Superclass::m_MovingImage->TransformIndexToPhysicalPoint(oindex,pt);
	//this->ComputeImageDerivatives(pt,fixedGradient);
	float mag1=0;
	for (int imd=0; imd<ImageDimension; imd++) 
	  {
	    mag1+=fixedGradient[imd]*fixedGradient[imd];
	    fixedGradientNorm[imd]=fixedGradient[imd];
	  }
	mag1=sqrt(mag1);
	if (mag1 > 1.e-5 ) fixedGradientNorm/=mag1;  
	
	for (int imd=0; imd<ImageDimension; imd++)
	  {
	    fdvec1[imd]=fixedGradientNorm[imd]*spacing[imd]*0.6;
	    fdvec2[imd]=fixedGradientNorm[imd]*(-1.)*spacing[imd]*0.6;
	  }
	if (m_ZeroInZ) { fdvec1[ImageDimension-1]=0; fdvec2[ImageDimension-1]=0; }
	this->GetValueAndDerivative(oindex,nccp1,fdvec1);
	this->GetValueAndDerivative(oindex,nccm1,fdvec2);
	//nccp1=this->GetValue(fdvec1,oindex);
	//nccm1=this->GetValue(fdvec2,oindex);
	float sign=nccp1-nccm1; 
        ///double denominator = mag1;
	double denominator = vnl_math_sqr( sign ) +  mag1;// from demons
	if (denominator < 1.e-12 || ! Superclass::m_NormalizeGradient )  denominator = (1.0);
	//	for (int imd=0; imd<ImageDimension; imd++) update[imd]=sign*fixedGradientNorm[imd];
	for (int imd=0; imd<ImageDimension; imd++) update[imd]=sign*fixedGradient[imd]/denominator*spacing[imd];
	loce+=(nccp1+nccm1);
      }

    if (m_ZeroInZ) { update[ImageDimension-1]=0.0; }
    loce/=(2.0*(float)ImageDimension);
    this-> Superclass::m_Energy+=loce/m_NormalizeMetric;
    //    std::cout << " Loce " << loce << std::endl;
    
    if (this->m_MetricImage) this->Superclass::m_MetricImage->SetPixel(oindex,loce);
    
    if (ImageDimension == 2)
      {
	
	if (this->m_MetricImage &&  
	    oindex[0] == this->Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[0]-5 && 
	    oindex[1] == this->Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[1]-5 )
	  { 
	    this->WriteImages();
	  }
      }
    else if (ImageDimension == 3)
      {
	
	if (this->m_MetricImage &&  
	    oindex[0] == this->Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[0]-5 && 
	    oindex[1] == this->Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[1]-5 && 
	    oindex[2] == this->Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[2]-5 )
	  { 
	    this->WriteImages();
	  }
      }

    //    for (int imd=0; imd<ImageDimension; imd++) mag1+=update[imd]*update[imd];
    //  mag1=sqrt(mag1);
    // if (mag1 > 1.e-5) update=update*(1.0/mag1); 

    return update*this->Superclass:: Superclass::m_GradientStep;
  }

  void WriteImages()
  {
  }
  

  void SetOpticalFlow(bool b){ m_OpticalFlow = b; }

  void SetZeroInZ(bool b) {m_ZeroInZ=b;}

protected:

  SectionMutualInformationRegistrationFunction();
  virtual ~SectionMutualInformationRegistrationFunction() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** FixedImage image neighborhood iterator type. */
  typedef ConstNeighborhoodIterator<FixedImageType> FixedImageNeighborhoodIteratorType;
  
  /** A global data type for this class of equation. Used to store
   * iterators for the fixed image. */
  struct GlobalDataStruct
   {
   FixedImageNeighborhoodIteratorType   m_FixedImageIterator;
   };


  /**
   * A fixed image spatial sample consists of the fixed domain point 
   * and the fixed image value at that point. */
  class FixedImageSpatialSample
  {
  public:
    FixedImageSpatialSample():FixedImageValue(0.0)
    { FixedImagePointValue.Fill(0.0); }
    ~FixedImageSpatialSample() {};

    IndexType                     FixedImageIndex;
    FixedImagePointType           FixedImagePointValue;
    double                        FixedImageValue;
    unsigned int                  FixedImageParzenWindowIndex;
  };

  /** FixedImageSpatialSample typedef support. */
  typedef std::vector<FixedImageSpatialSample>  FixedImageSpatialSampleContainer;

  /** Container to store a set of points and fixed image values. */
  FixedImageSpatialSampleContainer    m_FixedImageSamples;

  /** Uniformly select a sample set from the fixed image domain. */
  virtual void SampleFixedImageDomain( 
    FixedImageSpatialSampleContainer& samples);

  virtual void SampleFixedImageSlice(  FixedImageSpatialSampleContainer& samples, unsigned int );

  void SampleFixedImageDomainLocal( 
    FixedImageSpatialSampleContainer& samples, typename TFixedImage::IndexType index);

  /** Transform a point from FixedImage domain to MovingImage domain.
   * This function also checks if mapped point is within support region. */
  virtual void TransformPoint( unsigned int sampleNumber,
                               MovingImagePointType& mappedPoint, bool& sampleWithinSupportRegion,
                               double& movingImageValue ) const;

  unsigned int GetFixedValueIndex(double fixedImageValue, unsigned int slc)
  {

    // Determine parzen window arguments (see eqn 6 of Avants paper [2]).    
    double fixedImageParzenWindowTerm =
      fixedImageValue / m_FixedImageBinSize[slc] - m_FixedImageNormalizedMin[slc];
    unsigned int fixedImageParzenWindowIndex = 
      static_cast<unsigned int>( floor( fixedImageParzenWindowTerm ) );
    // Make sure the extreme values are in valid bins     
    if ( fixedImageParzenWindowIndex < 2 )
      {
	fixedImageParzenWindowIndex = 2;
      }
    else if ( fixedImageParzenWindowIndex > (m_NumberOfHistogramBins - 3) )
      {
	fixedImageParzenWindowIndex = m_NumberOfHistogramBins - 3;
      }

    return fixedImageParzenWindowIndex;
  }


  double GetFixedImageValue(IndexType oindex, DerivativeType derivative)
  {
    double fixedImageValue=0;
    PointType mappedPoint;
    for(int j = 0; j < ImageDimension; j++ )
      {
	mappedPoint[j] = double( oindex[j] ) * this->Superclass::m_FixedImage->GetSpacing()[j] + 
	  this->Superclass::m_FixedImage->GetOrigin()[j];
	mappedPoint[j] += derivative[j];
      }
    if( m_FixedImageInterpolator->IsInsideBuffer( mappedPoint ) )
      {
	fixedImageValue = m_FixedImageInterpolator->Evaluate( mappedPoint );
      }
    
    return fixedImageValue;
  }
  
  double GetMovingImageValue(IndexType oindex, DerivativeType derivative)
  {
    double movingImageValue=0;
    PointType mappedPoint;
    for(int j = 0; j < ImageDimension; j++ )
      {
	mappedPoint[j] = double( oindex[j] ) * this->Superclass::m_FixedImage->GetSpacing()[j] + 
	  this->Superclass::m_FixedImage->GetOrigin()[j];
	float temp = derivative[j];
	if (j == ( ImageDimension-1 ) && !m_ZeroInZ) 
	  {
	    if (temp > 0.5) temp=1.0;
	    if (temp < -0.5) temp=-1.0;
	  }
	mappedPoint[j] += temp;
      }
    if( m_MovingImageInterpolator->IsInsideBuffer( mappedPoint ) )
      {
	movingImageValue = m_MovingImageInterpolator->Evaluate( mappedPoint );
      }
    
    return movingImageValue;
  }

  unsigned int GetMovingValueIndex(double movingImageValue, unsigned int slc)
  {
    
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
    return movingImageParzenWindowIndex;
  }


private:

  SectionMutualInformationRegistrationFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented


  /** The marginal PDFs are stored as std::vector. */
  typedef float PDFValueType;
  typedef std::vector<PDFValueType> MarginalPDFType;

  /** The fixed image marginal PDF. */
  mutable MarginalPDFType m_FixedImageMarginalPDF[1000];

  /** The moving image marginal PDF. */
  mutable MarginalPDFType m_MovingImageMarginalPDF[1000];

  bool m_ZeroInZ;

  /** Typedef for the joint PDF and PDF derivatives are stored as ITK Images. */
  typedef Image<PDFValueType,3> JointPDFType;
  typedef Image<PDFValueType,3> JointPDFDerivativesType;
  typedef JointPDFType::IndexType                JointPDFIndexType;
  typedef JointPDFType::PixelType                JointPDFValueType;
  typedef JointPDFType::RegionType              JointPDFRegionType;
  typedef JointPDFType::SizeType                JointPDFSizeType;
  typedef JointPDFDerivativesType::IndexType    JointPDFDerivativesIndexType;
  typedef JointPDFDerivativesType::PixelType    JointPDFDerivativesValueType;
  typedef JointPDFDerivativesType::RegionType    JointPDFDerivativesRegionType;
  typedef JointPDFDerivativesType::SizeType      JointPDFDerivativesSizeType;

  /** The joint PDF and PDF derivatives. */
  typename JointPDFType::Pointer m_JointPDF;
  typename JointPDFDerivativesType::Pointer m_JointPDFDerivatives;

  unsigned long m_NumberOfSpatialSamples;
  unsigned long m_NumberOfParameters;

  /** Variables to define the marginal and joint histograms. */
  unsigned long m_NumberOfHistogramBins;
  vnl_vector<double> m_MovingImageNormalizedMin;
  vnl_vector<double> m_FixedImageNormalizedMin;
  vnl_vector<double> m_FixedImageBinSize;
  vnl_vector<double> m_MovingImageBinSize;

  /** Typedefs for BSpline kernel and derivative functions. */
  typedef BSplineKernelFunction<3> CubicBSplineFunctionType;
  typedef BSplineDerivativeKernelFunction<3> 
  CubicBSplineDerivativeFunctionType;

  /** Cubic BSpline kernel for computing Parzen histograms. */
  typename CubicBSplineFunctionType::Pointer m_CubicBSplineKernel;
  typename CubicBSplineDerivativeFunctionType::Pointer m_CubicBSplineDerivativeKernel;

  /** Precompute fixed image parzen window indices. */
  virtual void ComputeFixedImageParzenWindowIndices( FixedImageSpatialSampleContainer& samples , unsigned int);

  /**
   * Types and variables related to image derivative calculations.
   * If a BSplineInterpolationFunction is used, this class obtain
   * image derivatives from the BSpline interpolator. Otherwise, 
   * image derivatives are computed using central differencing.
   */
  typedef CovariantVector< double,
                           itkGetStaticConstMacro(ImageDimension) > ImageDerivativesType;

  /** Compute image derivatives at a point. */
  virtual void ComputeImageDerivatives( const MovingImagePointType& mappedPoint,
                                        ImageDerivativesType& gradient ) const;

  /** Boolean to indicate if the interpolator BSpline. */
  bool m_InterpolatorIsBSpline;

  // boolean to determine if we use mono-modality assumption
  bool m_OpticalFlow;

  /** Typedefs for using BSpline interpolator. */
  typedef  BSplineInterpolateImageFunction<MovingImageType,
                   CoordinateRepresentationType> BSplineInterpolatorType;

  /** Pointer to BSplineInterpolator. */
  typename BSplineInterpolatorType::Pointer m_BSplineInterpolator;

  /** Typedefs for using central difference calculator. */
  typedef CentralDifferenceImageFunction<MovingImageType,
                                         CoordinateRepresentationType> DerivativeFunctionType;

  /** Pointer to central difference calculator. */
  typename DerivativeFunctionType::Pointer m_DerivativeCalculator;


  /** Compute PDF derivative contribution for each parameter. */
  virtual void ComputePDFDerivatives( unsigned int sampleNumber,
                                      int movingImageParzenWindowIndex,
                                      const ImageDerivativesType& movingImageGradientValue,
                                      double cubicBSplineDerivativeValue ) const;

  /**
   * Types and variables related to BSpline deformable transforms.
   * If the transform is of type third order BSplineDeformableTransform,
   * then we can speed up the metric derivative calculation by
   * only inspecting the parameters within the support region
   * of a mapped point.
   */

  /** Boolean to indicate if the transform is BSpline deformable. */
  bool m_TransformIsBSpline;

  /** The number of BSpline parameters per image dimension. */
  long m_NumParametersPerDim;

  /** 
   * The number of BSpline transform weights is the number of
   * of parameter in the support region (per dimension ). */   
  unsigned long m_NumBSplineWeights;


  /** 
   * Enum of the deformabtion field spline order. 
   */
  enum { DeformationSplineOrder = 3 };

  /**
   * Typedefs for the BSplineDeformableTransform.
   */
  typedef BSplineDeformableTransform<
    CoordinateRepresentationType,
    ::itk::GetImageDimension<FixedImageType>::ImageDimension,
    DeformationSplineOrder> BSplineTransformType;
  typedef typename BSplineTransformType::WeightsType  BSplineTransformWeightsType;
  typedef typename BSplineTransformType::ParameterIndexArrayType  BSplineTransformIndexArrayType;

  /**
   * Variables used when transform is of type BSpline deformable.
   */
  typename BSplineTransformType::Pointer m_BSplineTransform;

  /**
   * Cache pre-transformed points, weights, indices and 
   * within support region flag.
   */
  typedef typename BSplineTransformWeightsType::ValueType WeightsValueType;
  typedef          Array2D<WeightsValueType> BSplineTransformWeightsArrayType;
  typedef typename BSplineTransformIndexArrayType::ValueType IndexValueType;
  typedef          Array2D<IndexValueType> BSplineTransformIndicesArrayType;
  typedef          std::vector<MovingImagePointType> MovingImagePointArrayType;
  typedef          std::vector<bool> BooleanArrayType;

  BSplineTransformWeightsArrayType      m_BSplineTransformWeightsArray;
  BSplineTransformIndicesArrayType      m_BSplineTransformIndicesArray;
  MovingImagePointArrayType             m_PreTransformPointsArray;
  BooleanArrayType                      m_WithinSupportRegionArray;
  
  typename TFixedImage::SpacingType                  m_FixedImageSpacing;
  typename TFixedImage::PointType                  m_FixedImageOrigin;

  typedef FixedArray<unsigned long, 
    ::itk::GetImageDimension<FixedImageType>::ImageDimension> ParametersOffsetType;
  ParametersOffsetType                  m_ParametersOffset;

  mutable TransformPointer    m_Transform;
  InterpolatorPointer         m_Interpolator;
  InterpolatorPointer             m_FixedImageInterpolator;
  InterpolatorPointer             m_MovingImageInterpolator;

  GradientCalculatorPointer       m_FixedImageGradientCalculator;
  GradientCalculatorPointer       m_MovingImageGradientCalculator;

  FixedImagePointer m_FixedImageMask;
  MovingImagePointer m_MovingImageMask;

  double m_NormalizeMetric;
  float m_Normalizer;

  typename ImageSliceType::Pointer m_FixedImageSlice;
  typename ImageSliceType::Pointer m_MovingImageSlice;

  unsigned int m_NumberOfSlices;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSectionMutualInformationRegistrationFunction.cxx"
#endif

#endif


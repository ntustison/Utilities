#ifndef MZ_itkESMImageToImageMetric_H_
#define MZ_itkESMImageToImageMetric_H_

#include "itkESMTransform.h"

#include <itkImageBase.h>
#include <itkInterpolateImageFunction.h>
#include <itkImageToImageMetric.h>
#include <itkExceptionObject.h>
#include <itkGradientRecursiveGaussianImageFilter.h>
#include <itkSpatialObject.h>
#include <itkCentralDifferenceImageFunction.h>
#include <itkCovariantVector.h>

#include <itkMultiThreader.h>

namespace itk
{

/** \class ESMImageToImageMetric
 * \brief Computes similarity between regions of two images.
 *
 * This Class is templated over the type of the two input images.
 * It expects a Transform and an Interpolator to be plugged in.
 * This particular class is the base class for a hierarchy of
 * similarity metrics.
 *
 * This class computes a value that measures the similarity
 * between the Fixed image and the transformed Moving image.
 * The Interpolator is used to compute intensity values on
 * non-grid positions resulting from mapping points through
 * the Transform.
 *
 *
 * \ingroup RegistrationMetrics
 *
 */

template <class TFixedImage,  class TMovingImage>
class ITK_EXPORT ESMImageToImageMetric
   : public ImageToImageMetric<TFixedImage,TMovingImage>
{
public:
   /** Standard class typedefs. */
   typedef ESMImageToImageMetric        Self;
   typedef ImageToImageMetric<TFixedImage,
      TMovingImage>                     Superclass;
   typedef SmartPointer<Self>           Pointer;
   typedef SmartPointer<const Self>     ConstPointer;

   /** Type used for representing point components  */
   typedef typename Superclass::CoordinateRepresentationType CoordinateRepresentationType;

   /** Run-time type information (and related methods). */
   itkTypeMacro(ESMImageToImageMetric, ImageToImageMetric);

   /**  Type of the moving Image. */
   typedef typename Superclass::MovingImageType          MovingImageType;

   /**  Type of the fixed Image. */
   typedef typename Superclass::FixedImageType           FixedImageType;

   /** Constants for the image dimensions */
   itkStaticConstMacro(ImageDimension, unsigned int, TFixedImage::ImageDimension);

   typedef itk::ImageRegion<ImageDimension>              ImageRegionType;

   /**  Type of the Transform Base class */
   typedef ESMTransform<CoordinateRepresentationType,
      itkGetStaticConstMacro(ImageDimension)>            ESMTransformType;

   typedef typename ESMTransformType::Pointer            ESMTransformPointer;
   typedef typename ESMTransformType::PointType          PointType;
   typedef typename ESMTransformType::ParametersType     TransformParametersType;
   typedef typename ESMTransformType::JacobianType       TransformJacobianType;

   /** Index and Point typedef support. */
   typedef itk::Index<ImageDimension>                    ImageIndexType;

   typedef std::vector<ImageIndexType>                   ImageIndexContainer;

   /**  Type of the Interpolator Base class */
   typedef InterpolateImageFunction< MovingImageType,
      CoordinateRepresentationType >                     InterpolatorType;

   /** Gaussian filter to compute the gradient of the Moving Image */
   typedef typename Superclass::GradientImageType        GradientImageType;

   /**  Type for the mask of the fixed image. Only pixels that are "inside"
        this mask will be considered for the computation of the metric */
   typedef typename Superclass::FixedImageMaskType       FixedImageMaskType;
   typedef typename Superclass::MovingImageMaskType      MovingImageMaskType;

   /**  Type of the measure. */
   typedef typename Superclass::MeasureType                    MeasureType;

   /**  Type of the derivative. */
   typedef typename Superclass::DerivativeType                 DerivativeType;

   /**  Type of the Hessian. */
   typedef vnl_matrix<double>                                  HessianType;

   /**  Type of the parameters. */
   typedef typename Superclass::ParametersType                 ParametersType;

   /** This method returns the evaluation of the cost function corresponding
    * to the specified parameters.   */
   virtual void GetValueDerivativeAndHessian( const ParametersType & parameters,
                                              MeasureType & value,
                                              DerivativeType & derivative,
                                              HessianType & hessian ) const = 0;

   /** Connect the Transform. */
   itkSetObjectMacro( ESMTransform, ESMTransformType );

   /** Get a pointer to the Transform.  */
   itkGetConstObjectMacro( ESMTransform, ESMTransformType );

   /** Get the number of pixels considered in the computation. */
   unsigned long GetNumberOfMovingImageSamples( void )
   {
      return this->GetNumberOfPixelsCounted();
   }
   itkGetConstReferenceMacro( NumberOfPixelsCounted, unsigned long );

   /** Set the region over which the metric will be computed */
   void SetFixedImageRegion( const ImageRegionType reg );

   /** Set/Get the moving image mask. */
   using Superclass::SetMovingImageMask;

   void SetMovingImageMask( const MovingImageMaskType* mask )
   {
      // This is a work-around for the fact that the non-optimized
      // itk registration imaplementation is not const-correct
      // w.r.t. to the masks
      this->Superclass::SetMovingImageMask(const_cast<MovingImageMaskType*>(mask));
   }

   /** Set/Get the fixed image mask. */
   using Superclass::SetFixedImageMask;

   void SetFixedImageMask( const FixedImageMaskType* mask )
   {
      // This is a work-around for the fact that the non-optimized
      // itk registration imaplementation is not const-correct
      // w.r.t. to the masks
      this->Superclass::SetFixedImageMask(const_cast<FixedImageMaskType*>(mask));
   }

   /** Set the fixed image indexes to be used as the samples when
    *   computing the match metric */
   void SetFixedImageIndexes( const ImageIndexContainer & indexes );
   void SetUseFixedImageIndexes( bool useIndex );
   itkGetConstReferenceMacro( UseFixedImageIndexes, bool );

   /** Set/Get number of threads to use for computations. */
   void SetNumberOfThreads( unsigned int numberOfThreads );
   itkGetConstReferenceMacro( NumberOfThreads, unsigned int );

   /** Computes the gradient image and assigns it to m_GradientImage */
   virtual void ComputeGradient( void );

   /** Set the parameters defining the Transform. */
   void SetTransformParameters( const ParametersType & parameters ) const;

   void IncrementalESMTransformUpdate( const ParametersType & update ) const;

   /** Return the number of parameters required by the Transform */
   unsigned int GetNumberOfParameters( void ) const
   {
      return m_ESMTransform->GetNumberOfParameters();
   }

   /** Initialize the Metric by making sure that all the components
    *  are present and plugged together correctly     */
   virtual void Initialize( void ) throw ( ExceptionObject );

   /** Initialize the components related to supporting multiple threads */
   virtual void MultiThreadingInitialize( void ) throw ( ExceptionObject );

   /** Number of spatial samples to used to compute metric
    *   This sets the number of samples.  */
   virtual void SetNumberOfFixedImageSamples( unsigned long numSamples );
   itkGetConstReferenceMacro( NumberOfFixedImageSamples, unsigned long );

   /** Number of spatial samples to used to compute metric
    *   This sets the number of samples.  */
   void SetNumberOfSpatialSamples( unsigned long num )
   {
      this->SetNumberOfFixedImageSamples( num );
   }
   unsigned long GetNumberOfSpatialSamples( void )
   {
      return this->GetNumberOfFixedImageSamples();
   }

   /** Select whether the metric will be computed using all the pixels on the
    * fixed image region, or only using a set of randomly selected pixels.
    * This value override IntensityThreshold, Masks, and SequentialSampling. */
   void SetUseAllPixels( bool useAllPixels );
   void UseAllPixelsOn( void )
   {
      this->SetUseAllPixels( true );
   }
   void UseAllPixelsOff( void )
   {
      this->SetUseAllPixels( false );
   }
   itkGetConstReferenceMacro( UseAllPixels, bool );

   /** If set to true, then every pixel in the fixed image will be scanned to
    * determine if it should be used in registration metric computation.  A
    * pixel will be chosen if it meets any mask or threshold limits set.  If
    * set to false, then UseAllPixels will be set to false. */
   void SetUseSequentialSampling( bool sequentialSampling );
   itkGetConstReferenceMacro( UseSequentialSampling, bool );

   /** Reinitialize the seed of the random number generator that selects the
    * sample of pixels used for estimating the image histograms and the joint
    * histogram. By nature, this metric is not deterministic, since at each run
    * it may select a different set of pixels. By initializing the random number
    * generator seed to the same value you can restore determinism. On the other
    * hand, calling the method ReinitializeSeed() without arguments will use the
    * clock from your machine in order to have a very random initialization of
    * the seed. This will indeed increase the non-deterministic behavior of the
    * metric. */
   void ReinitializeSeed();
   void ReinitializeSeed( int seed );

#ifdef ITK_USE_CONCEPT_CHECKING
   /** Begin concept checking */
   itkConceptMacro(SameDimensionCheck1,
                   (Concept::SameDimension<TFixedImage::ImageDimension,TMovingImage::ImageDimension>));
   /** End concept checking */
#endif

protected:
   ESMImageToImageMetric();
   virtual ~ESMImageToImageMetric();

   void PrintSelf(std::ostream& os, Indent indent) const;

#ifdef ITK_USE_OPTIMIZED_REGISTRATION_METHODS
   typedef typename Superclass::FixedImageSamplePoint FixedImageSamplePoint;
#else
   /** \class FixedImageSamplePoint
    * A fixed image spatial sample consists of the fixed domain point
    * and the fixed image value at that point. */
   /// @cond
   class FixedImageSamplePoint
   {
   public:
      FixedImageSamplePoint()
      {
         point.Fill(0.0);
         value = 0;
         valueIndex = 0;
      }
      ~FixedImageSamplePoint() {};

   public:
      PointType                     point;
      double                        value;
      unsigned int                  valueIndex;
   };
   /// @endcond
#endif

   bool                      m_UseFixedImageIndexes;
   ImageIndexContainer       m_FixedImageIndexes;

   /** FixedImageSamplePoint typedef support. */
   typedef std::vector<FixedImageSamplePoint> FixedImageSampleContainer;

   /** Uniformly select a sample set from the fixed image domain. */
   virtual void SampleFixedImageDomain( FixedImageSampleContainer & samples) const;

   virtual void SampleFixedImageIndexes( FixedImageSampleContainer & samples) const;

   /** Gather all the pixels from the fixed image domain. */
   virtual void SampleFullFixedImageDomain( FixedImageSampleContainer & samples) const;

   /** Container to store a set of points and fixed image values. */
   FixedImageSampleContainer   m_FixedImageSamples;

   unsigned long               m_NumberOfParameters;
   mutable ParametersType      m_Parameters;

   unsigned long               m_NumberOfFixedImageSamples;
   //m_NumberOfPixelsCounted must be mutable because the const
   //thread consolidation functions merge each threads valus
   //onto this accumulator variable.
   mutable unsigned long       m_NumberOfPixelsCounted;

   /** Main transform to be used in thread = 0 */
   ESMTransformPointer         m_ESMTransform;

   /** Copies of Transform helpers per thread (N-1 of them, since m_Transform
    * will do the work for thread=0. */
   ESMTransformPointer       * m_ThreaderESMTransform;

   unsigned int                m_NumberOfThreads;

   bool                        m_UseAllPixels;
   bool                        m_UseSequentialSampling;

   bool                        m_ReseedIterator;

   int                         m_RandomSeed;

   typename GradientImageType::Pointer  m_FixedGradientImage;

   typedef          std::vector<PointType>                     PointArrayType;
   typedef          std::vector<bool>                          BooleanArrayType;

   /** Typedefs for using central difference calculator. */
   typedef CentralDifferenceImageFunction<MovingImageType,
      CoordinateRepresentationType>                            DerivativeFunctionType;
   typedef CovariantVector< double,
      itkGetStaticConstMacro(ImageDimension) >                 ImageDerivativesType;
   typedef CentralDifferenceImageFunction<FixedImageType,
      CoordinateRepresentationType>                            FixedDerivativeFunctionType;

   /** Transform a point from FixedImage domain to MovingImage domain.
    * This function also checks if mapped point is within support region. */
   virtual void TransformPoint( unsigned int sampleNumber,
                                PointType& mappedPoint,
                                bool& sampleWithinSupportRegion,
                                double& movingImageValue,
                                unsigned int threadID ) const;

#ifdef ITK_USE_OPTIMIZED_REGISTRATION_METHODS
   using Superclass::TransformPointWithDerivatives;
#endif
   virtual void TransformPointWithDerivatives( unsigned int sampleNumber,
                                               PointType& mappedPoint,
                                               bool& sampleWithinSupportRegion,
                                               double& movingImageValue,
                                               ImageDerivativesType & movingGradient,
                                               ImageDerivativesType & fixedGradient,
                                               unsigned int threadID ) const;

   /** Pointer to central difference calculator. */
   typename DerivativeFunctionType::Pointer             m_DerivativeCalculator;
   typename FixedDerivativeFunctionType::Pointer        m_FixedDerivativeCalculator;

   /** Compute image derivatives at a point. */
#ifdef ITK_USE_OPTIMIZED_REGISTRATION_METHODS
   using Superclass::ComputeImageDerivatives;
#endif
   virtual void ComputeImageDerivatives( const PointType & mappedPoint,
                                         ImageDerivativesType & gradient ) const;

   virtual void ComputeFixedImageDerivatives( const PointType & fixedPoint,
                                              ImageDerivativesType & gradient ) const;


   /**
    * Types and variables related to multi-threading
    */

   typedef MultiThreader               MultiThreaderType;

   struct MultiThreaderParameterType
   {
      ESMImageToImageMetric               * metric;
   };

   MultiThreaderType::Pointer               m_Threader;
   MultiThreaderParameterType               m_ThreaderParameter;
   mutable unsigned int                   * m_ThreaderNumberOfMovingImageSamples;
   bool                                     m_WithinThreadPreProcess;
   bool                                     m_WithinThreadPostProcess;

   // threaded value comp stuff
   void                           GetValueMultiThreadedPreProcessInitiate( void ) const;
   void                           GetValueMultiThreadedInitiate( void ) const;
   void                           GetValueMultiThreadedPostProcessInitiate( void ) const;
   static ITK_THREAD_RETURN_TYPE  GetValueMultiThreadedPreProcess( void * arg );
   static ITK_THREAD_RETURN_TYPE  GetValueMultiThreaded( void * arg );
   static ITK_THREAD_RETURN_TYPE  GetValueMultiThreadedPostProcess( void * arg );

   virtual inline void       GetValueThread( unsigned int threadID ) const;
   virtual inline void       GetValueThreadPreProcess(
      unsigned int itkNotUsed(threadID),
      bool itkNotUsed(withinSampleThread) ) const
   { };
   virtual inline bool       GetValueThreadProcessSample(
      unsigned int itkNotUsed(threadID),
      unsigned long itkNotUsed(fixedImageSample),
      const PointType & itkNotUsed(mappedPoint),
      double itkNotUsed(movingImageValue)) const
   {
      return false;
   }

   virtual inline void       GetValueThreadPostProcess(
      unsigned int itkNotUsed(threadID), bool itkNotUsed(withinSampleThread) ) const {}

   // threaded value and deriv comp stuff
   void                          GetValueAndDerivativeMultiThreadedPreProcessInitiate( void ) const;
   void                          GetValueAndDerivativeMultiThreadedInitiate( void ) const;
   void                          GetValueAndDerivativeMultiThreadedPostProcessInitiate( void ) const;
   static ITK_THREAD_RETURN_TYPE GetValueAndDerivativeMultiThreadedPreProcess( void * arg );

   static ITK_THREAD_RETURN_TYPE GetValueAndDerivativeMultiThreaded( void * arg );

   static ITK_THREAD_RETURN_TYPE GetValueAndDerivativeMultiThreadedPostProcess( void * arg );

   virtual inline void  GetValueAndDerivativeThread(unsigned int threadID) const;
   virtual inline void  GetValueAndDerivativeThreadPreProcess(
      unsigned int itkNotUsed(threadID), bool itkNotUsed(withinSampleThread)) const {}
#ifdef ITK_USE_OPTIMIZED_REGISTRATION_METHODS
   using Superclass::GetValueAndDerivativeThreadProcessSample;
#endif
   virtual inline bool  GetValueAndDerivativeThreadProcessSample(
      unsigned int itkNotUsed(threadID),
      unsigned long itkNotUsed(fixedImageSample),
      const PointType & itkNotUsed(mappedPoint),
      double itkNotUsed(movingImageValue),
      const ImageDerivativesType & itkNotUsed(movingImageGradientValue),
      const ImageDerivativesType & itkNotUsed(fixedImageGradientValue) ) const
   {
      return false;
   }

   virtual inline void  GetValueAndDerivativeThreadPostProcess(
      unsigned int itkNotUsed(threadID), bool itkNotUsed(withinSampleThread) ) const {}

   // threaded value and deriv comp stuff
   void                GetValueDerivativeAndHessianMultiThreadedPreProcessInitiate( void) const;
   void                GetValueDerivativeAndHessianMultiThreadedInitiate( void ) const;
   void                GetValueDerivativeAndHessianMultiThreadedPostProcessInitiate( void) const;
   static ITK_THREAD_RETURN_TYPE GetValueDerivativeAndHessianMultiThreadedPreProcess(void * arg);

   static ITK_THREAD_RETURN_TYPE GetValueDerivativeAndHessianMultiThreaded(void * arg);

   static ITK_THREAD_RETURN_TYPE GetValueDerivativeAndHessianMultiThreadedPostProcess(void * arg);

   virtual inline void  GetValueDerivativeAndHessianThread(unsigned int threadID) const;
   virtual inline void  GetValueDerivativeAndHessianThreadPreProcess(
      unsigned int itkNotUsed(threadID), bool itkNotUsed(withinSampleThread)) const {}
   virtual inline bool  GetValueDerivativeAndHessianThreadProcessSample(
      unsigned int itkNotUsed(threadID),
      unsigned long itkNotUsed(fixedImageSample),
      const PointType & itkNotUsed(mappedPoint),
      double itkNotUsed(movingImageValue),
      const ImageDerivativesType & itkNotUsed(movingImageGradientValue),
      const ImageDerivativesType & itkNotUsed(fixedImageGradientValue) ) const
   {
      return false;
   }

   virtual inline void  GetValueDerivativeAndHessianThreadPostProcess(
      unsigned int itkNotUsed(threadID), bool itkNotUsed(withinSampleThread) ) const {}

   /** Synchronizes the threader transforms with the transform
    *   member variable.
    */
   void SynchronizeTransforms() const;

private:
   ESMImageToImageMetric(const Self&); //purposely not implemented
   void operator=(const Self&); //purposely not implemented

   // m_ESMTransform should be used in lieu of m_Transform m_Transform
   virtual void SetTransform(typename Superclass::TransformType*){assert(false);};
};

} // end namespace itk

#include "itkESMImageToImageMetric.txx"

#endif

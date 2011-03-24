#ifndef MZ_ESMImageToImageMetric_TXX_
#define MZ_ESMImageToImageMetric_TXX_


#include "itkESMImageToImageMetric.h"

#include <itkImageRandomConstIteratorWithIndex.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>

namespace itk
{

/**
 * Constructor
 */
template <class TFixedImage, class TMovingImage>
ESMImageToImageMetric<TFixedImage,TMovingImage>
::ESMImageToImageMetric()
   :Superclass()
{
   m_NumberOfFixedImageSamples = 50000;
   m_UseAllPixels = false;
   m_UseSequentialSampling = false;
   m_UseFixedImageIndexes = false;
   m_ReseedIterator = false;
   m_RandomSeed = -1;

   m_Threader = MultiThreaderType::New();
   m_ThreaderParameter.metric = this;
   m_ThreaderNumberOfMovingImageSamples = NULL;
   m_WithinThreadPreProcess = false;
   m_WithinThreadPostProcess = false;

   m_ESMTransform         = NULL; // has to be provided by the user.
   m_ThreaderESMTransform = NULL; // constructed at initialization.

   // For convenience
   this->m_Interpolator  = itk::LinearInterpolateImageFunction<TMovingImage>::New();

   m_DerivativeCalculator = NULL;
   m_FixedDerivativeCalculator = NULL;

   m_NumberOfThreads = m_Threader->GetNumberOfThreads();
}

template <class TFixedImage, class TMovingImage>
ESMImageToImageMetric<TFixedImage,TMovingImage>
::~ESMImageToImageMetric()
{

   if(m_ThreaderNumberOfMovingImageSamples != NULL)
   {
      delete [] m_ThreaderNumberOfMovingImageSamples;
   }
   m_ThreaderNumberOfMovingImageSamples = NULL;

   if(m_ThreaderESMTransform != NULL)
   {
      delete [] m_ThreaderESMTransform;
   }
   m_ThreaderESMTransform = NULL;
}

/**
 * Set the number of threads. This will be clamped by the
 * multithreader, so we must check to see if it is accepted.
 */
template <class TFixedImage, class TMovingImage>
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::SetNumberOfThreads( unsigned int numberOfThreads )
{
   m_Threader->SetNumberOfThreads( numberOfThreads);
   m_NumberOfThreads = m_Threader->GetNumberOfThreads();
}

/**
 * Set the parameters that define a unique transform
 */
template <class TFixedImage, class TMovingImage>
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::SetTransformParameters( const ParametersType & parameters ) const
{
   if( !m_ESMTransform )
   {
      itkExceptionMacro(<<"Transform has not been assigned");
   }
   m_ESMTransform->SetParameters( parameters );

   m_Parameters = parameters;
}

template <class TFixedImage, class TMovingImage>
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::SetNumberOfFixedImageSamples( unsigned long numSamples )
{
   if( numSamples != m_NumberOfFixedImageSamples )
   {
      m_NumberOfFixedImageSamples = numSamples;
      if( m_NumberOfFixedImageSamples != this->GetFixedImageRegion().GetNumberOfPixels() )
      {
         this->SetUseAllPixels( false );
      }
      this->Modified();
   }
}

template <class TFixedImage, class TMovingImage>
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::SetFixedImageIndexes( const ImageIndexContainer & indexes )
{
   this->SetUseFixedImageIndexes( true );
   m_NumberOfFixedImageSamples = indexes.size();
   m_FixedImageIndexes.resize( m_NumberOfFixedImageSamples );
   for(unsigned int i=0; i<m_NumberOfFixedImageSamples; i++)
   {
      m_FixedImageIndexes[i] = indexes[i];
   }
}

template <class TFixedImage, class TMovingImage>
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::SetUseFixedImageIndexes( bool useIndexes )
{
   if( useIndexes != m_UseFixedImageIndexes )
   {
      m_UseFixedImageIndexes = useIndexes;
      if( m_UseFixedImageIndexes )
      {
         this->SetUseAllPixels( false );
      }
      else
      {
         this->Modified();
      }
   }
}

template <class TFixedImage, class TMovingImage>
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::SetFixedImageRegion( const ImageRegionType reg )
{
   if( reg != this->GetFixedImageRegion() )
   {
      Superclass::SetFixedImageRegion(reg);
      if( this->GetUseAllPixels() )
      {
         this->SetNumberOfFixedImageSamples( this->GetFixedImageRegion().GetNumberOfPixels() );
      }
   }
}

template <class TFixedImage, class TMovingImage>
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::SetUseAllPixels( bool useAllPixels )
{
  if( useAllPixels != m_UseAllPixels )
    {
    m_UseAllPixels = useAllPixels;
    if( m_UseAllPixels )
      {
      this->SetNumberOfFixedImageSamples( this->GetFixedImageRegion().GetNumberOfPixels() );
      this->SetUseSequentialSampling( true );
      }
    else
      {
      this->SetUseSequentialSampling( false );
      this->Modified();
      }
    }
}


template <class TFixedImage, class TMovingImage>
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::SetUseSequentialSampling( bool useSequential )
{
  if( useSequential != m_UseSequentialSampling )
    {
    m_UseSequentialSampling = useSequential;
    if( !m_UseSequentialSampling )
      {
      this->SetUseAllPixels( false );
      }
    else
      {
      this->Modified();
      }
    }
}

template <class TFixedImage, class TMovingImage>
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::IncrementalESMTransformUpdate( const ParametersType & update ) const
{
   if( !m_ESMTransform )
   {
      itkExceptionMacro(<<"Transform has not been assigned");
   }
   m_ESMTransform->IncrementalUpdate( update );

   m_Parameters = m_ESMTransform->GetParameters();
}

/**
 * Initialize
 */
template <class TFixedImage, class TMovingImage>
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::Initialize(void) throw ( ExceptionObject )
{

   if( !this->m_ESMTransform )
   {
      itkExceptionMacro(<<"Transform is not present");
   }
   this->m_NumberOfParameters = this->m_ESMTransform->GetNumberOfParameters();

   if( !this->m_Interpolator )
   {
      itkExceptionMacro(<<"Interpolator is not present");
   }

   if( !this->m_MovingImage )
   {
      itkExceptionMacro(<<"MovingImage is not present");
   }

   if( !this->m_FixedImage )
   {
      itkExceptionMacro(<<"FixedImage is not present");
   }


   if( this->GetFixedImageRegion().GetNumberOfPixels() == 0 )
   {
      itkExceptionMacro(<<"FixedImageRegion is empty");
   }

   // If the image is provided by a source, update the source.
   if( this->m_MovingImage->GetSource() )
   {
      this->m_MovingImage->GetSource()->Update();
   }

   // If the image is provided by a source, update the source.
   if( this->m_FixedImage->GetSource() )
   {
      this->m_FixedImage->GetSource()->Update();
   }

   // Make sure the FixedImageRegion is within the FixedImage buffered region
   if ( !this->m_FixedImage->GetBufferedRegion().IsInside( this->GetFixedImageRegion() ) )
   {
      typename Superclass::FixedImageRegionType croppedRegion( this->GetFixedImageRegion() );
      if ( !croppedRegion.Crop( this->m_FixedImage->GetBufferedRegion() ) )
      {
         itkExceptionMacro(
            <<"FixedImageRegion does not overlap the fixed image buffered region" );
      }
      this->SetFixedImageRegion( croppedRegion );
   }

   this->m_Interpolator->SetInputImage( this->m_MovingImage );

   if ( this->m_ComputeGradient )
   {
      this->ComputeGradient();
   }

   // If there are any observers on the metric, call them to give the
   // user code a chance to set parameters on the metric
   this->InvokeEvent( InitializeEvent() );

}


/**
 * MultiThreading Initialize
 */
template <class TFixedImage, class TMovingImage>
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::MultiThreadingInitialize(void) throw ( ExceptionObject )
{

   m_Threader->SetNumberOfThreads( m_NumberOfThreads );

   if(m_ThreaderNumberOfMovingImageSamples != NULL)
   {
      delete [] m_ThreaderNumberOfMovingImageSamples;
   }
   m_ThreaderNumberOfMovingImageSamples = new unsigned int[m_NumberOfThreads-1];

   // Allocate the array of transform clones to be used in every thread
   if(m_ThreaderESMTransform != NULL)
   {
      delete [] m_ThreaderESMTransform;
   }
   m_ThreaderESMTransform = new ESMTransformPointer[m_NumberOfThreads-1];
   for( unsigned int ithread=0; ithread < m_NumberOfThreads-1; ++ithread)
   {
      // Create a copy of the main transform to be used in this thread.
      LightObject::Pointer anotherTransform = this->m_ESMTransform->CreateAnother();
      // This static_cast should always work since the pointer was created by
      // CreateAnother() called from the transform itself.
      ESMTransformType * transformCopy = static_cast< ESMTransformType * >( anotherTransform.GetPointer() );
      /** Set the fixed parameters first. Some transforms have parameters which depend on
          the values of the fixed parameters. For instance, the BSplineDeformableTransform
          checks the grid size (part of the fixed parameters) before setting the parameters. */
      transformCopy->SetFixedParameters( this->m_ESMTransform->GetFixedParameters() );
      transformCopy->SetParameters( this->m_ESMTransform->GetParameters() );
      this->m_ThreaderESMTransform[ithread] = transformCopy;
   }

   m_FixedImageSamples.resize( m_NumberOfFixedImageSamples );
   if( m_UseSequentialSampling )
   {
      //
      // Take all the pixels within the fixed image region)
      // to create the sample points list.
      //
      SampleFullFixedImageDomain( m_FixedImageSamples );
   }
   else
   {
      if( m_UseFixedImageIndexes )
      {
         //
         //  Use the list of indexes passed to the SetFixedImageIndexes
         //  member function .
         //
         SampleFixedImageIndexes( m_FixedImageSamples );
      }
      else
      {
         //
         // Uniformly sample the fixed image (within the fixed image region)
         // to create the sample points list.
         //
         SampleFixedImageDomain( m_FixedImageSamples );
      }
   }

   m_DerivativeCalculator = DerivativeFunctionType::New();
   m_FixedDerivativeCalculator = FixedDerivativeFunctionType::New();

#ifdef ITK_USE_ORIENTED_IMAGE_DIRECTION
   m_DerivativeCalculator->UseImageDirectionOn();
   m_FixedDerivativeCalculator->UseImageDirectionOn();
#endif

   m_DerivativeCalculator->SetInputImage( this->m_MovingImage );
   m_FixedDerivativeCalculator->SetInputImage( this->m_FixedImage );
}


/**
 * Uniformly sample the fixed image domain using a random walk
 */
template < class TFixedImage, class TMovingImage >
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::SampleFixedImageIndexes( FixedImageSampleContainer & samples ) const
{
   typename FixedImageSampleContainer::iterator iter;

   unsigned long len = m_FixedImageIndexes.size();
   if( len != m_NumberOfFixedImageSamples
       || samples.size() != m_NumberOfFixedImageSamples )
   {
      throw ExceptionObject(__FILE__, __LINE__,
                            "Index list size does not match desired number of samples" );
   }

   iter=samples.begin();
   for(unsigned long i=0; i<len; i++)
   {
      // Get sampled index
      ImageIndexType index = this->m_FixedImageIndexes[i];
      // Translate index to point
      this->m_FixedImage->TransformIndexToPhysicalPoint( index, (*iter).point );

      // Get sampled fixed image value
      (*iter).value = this->m_FixedImage->GetPixel( index );
      (*iter).valueIndex = 0;

      ++iter;
   }
}

/**
 * Sample the fixed image using a random walk
 */
template < class TFixedImage, class TMovingImage >
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::SampleFixedImageDomain( FixedImageSampleContainer & samples ) const
{
   if( samples.size() != m_NumberOfFixedImageSamples )
   {
      throw ExceptionObject(__FILE__, __LINE__,
                            "Sample size does not match desired number of samples" );
   }

   // Set up a random interator within the user specified fixed image region.
   typedef ImageRandomConstIteratorWithIndex<FixedImageType> RandomIterator;
   RandomIterator randIter( this->m_FixedImage, this->GetFixedImageRegion() );

   typename FixedImageSampleContainer::iterator iter;
   typename FixedImageSampleContainer::const_iterator end=samples.end();

   if( this->m_FixedImageMask.IsNotNull() )
   {
      PointType inputPoint;

      iter=samples.begin();
      unsigned long int samplesFound = 0;
      randIter.SetNumberOfSamples( m_NumberOfFixedImageSamples * 1000 );
      randIter.GoToBegin();
      while( iter != end )
      {
         if( randIter.IsAtEnd() )
         {
            // Must be a small mask since after many random samples we don't
            // have enough to fill the desired number.   So, we will replicate
            // the samples we've found so far to fill-in the desired number
            // of samples
            unsigned long int count = 0;
            while( iter != end )
            {
               (*iter).point = samples[count].point;
               (*iter).value = samples[count].value;
               (*iter).valueIndex = 0;
               ++count;
               if(count >= samplesFound)
               {
                  count = 0;
               }
               ++iter;
            }
            break;
         }

         // Get sampled index
         ImageIndexType index = randIter.GetIndex();
         // Check if the Index is inside the mask, translate index to point
         this->m_FixedImage->TransformIndexToPhysicalPoint( index, inputPoint );

         if( this->m_FixedImageMask.IsNotNull() )
         {
            double val;
            if( this->m_FixedImageMask->ValueAt( inputPoint, val ) )
            {
               if( val == 0 )
               {
                  ++randIter; // jump to another random position
                  continue;
               }
            }
            else
            {
               ++randIter; // jump to another random position
               continue;
            }
         }

         // Translate index to point
         (*iter).point = inputPoint;
         // Get sampled fixed image value
         (*iter).value = randIter.Get();
         (*iter).valueIndex = 0;

         ++samplesFound;
         ++randIter;
         ++iter;
      }
   }
   else
   {
      randIter.SetNumberOfSamples( m_NumberOfFixedImageSamples );
      randIter.GoToBegin();
      for( iter=samples.begin(); iter != end; ++iter )
      {
         // Get sampled index
         ImageIndexType index = randIter.GetIndex();
         // Translate index to point
         this->m_FixedImage->TransformIndexToPhysicalPoint( index,
                                                            (*iter).point );
         // Get sampled fixed image value
         (*iter).value = randIter.Get();
         (*iter).valueIndex = 0;

         // Jump to random position
         ++randIter;
      }
   }
}

/**
 * Sample the fixed image domain using all pixels in the Fixed image region
 */
template < class TFixedImage, class TMovingImage >
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::SampleFullFixedImageDomain( FixedImageSampleContainer& samples ) const
{

   if( samples.size() != m_NumberOfFixedImageSamples )
   {
      throw ExceptionObject(__FILE__, __LINE__,
                            "Sample size does not match desired number of samples" );
   }

   // Set up a region interator within the user specified fixed image region.
   typedef ImageRegionConstIteratorWithIndex<FixedImageType> RegionIterator;
   RegionIterator regionIter( this->m_FixedImage, this->GetFixedImageRegion() );

   regionIter.GoToBegin();

   typename FixedImageSampleContainer::iterator iter;
   typename FixedImageSampleContainer::const_iterator end=samples.end();

   if( this->m_FixedImageMask.IsNotNull() )
   {
      PointType inputPoint;

      // repeat until we get enough samples to fill the array
      iter=samples.begin();
      while( iter != end )
      {
         // Get sampled index
         ImageIndexType index = regionIter.GetIndex();
         // Check if the Index is inside the mask, translate index to point
         this->m_FixedImage->TransformIndexToPhysicalPoint( index, inputPoint );

         if( this->m_FixedImageMask.IsNotNull() )
         {
            // If not inside the mask, ignore the point
            if( !this->m_FixedImageMask->IsInside( inputPoint ) )
            {
               ++regionIter; // jump to next pixel
               if( regionIter.IsAtEnd() )
               {
                  regionIter.GoToBegin();
               }
               continue;
            }
         }

         // Translate index to point
         (*iter).point = inputPoint;
         // Get sampled fixed image value
         (*iter).value = regionIter.Get();
         (*iter).valueIndex = 0;

         ++regionIter;
         if( regionIter.IsAtEnd() )
         {
            regionIter.GoToBegin();
         }
         ++iter;
      }
   }
   else // not restricting sample throwing to a mask
   {
      for( iter=samples.begin(); iter != end; ++iter )
      {
         // Get sampled index
         ImageIndexType index = regionIter.GetIndex();

         // Translate index to point
         this->m_FixedImage->TransformIndexToPhysicalPoint( index,
                                                            (*iter).point );
         // Get sampled fixed image value
         (*iter).value = regionIter.Get();
         (*iter).valueIndex = 0;

         ++regionIter;
         if( regionIter.IsAtEnd() )
         {
            regionIter.GoToBegin();
         }
      }
   }
}

/**
 * Compute the gradient image and assign it to m_GradientImage.
 */
template <class TFixedImage, class TMovingImage>
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::ComputeGradient()
{
   {
   // Moving image gradient
   typedef GradientRecursiveGaussianImageFilter< MovingImageType,
         GradientImageType >  GradientImageFilterType;
   typename GradientImageFilterType::Pointer gradientFilter = GradientImageFilterType::New();

   gradientFilter->SetInput( this->m_MovingImage );

   const typename MovingImageType::SpacingType & spacing =
      this->m_MovingImage->GetSpacing();
   double maximumSpacing=0.0;
   for(unsigned int i=0; i<ImageDimension; i++)
   {
      if( spacing[i] > maximumSpacing )
      {
         maximumSpacing = spacing[i];
      }
   }
   gradientFilter->SetSigma( maximumSpacing );
   gradientFilter->SetNormalizeAcrossScale( true );
   gradientFilter->SetNumberOfThreads( m_NumberOfThreads );

#ifdef ITK_USE_ORIENTED_IMAGE_DIRECTION
   gradientFilter->SetUseImageDirection( true );
#endif

   gradientFilter->Update();

   this->m_GradientImage = gradientFilter->GetOutput();
   }

   {
   // Fixed image gradient
   typedef GradientRecursiveGaussianImageFilter< FixedImageType,
         GradientImageType > FixedGradientImageFilterType;
   typename FixedGradientImageFilterType::Pointer fixedGradientFilter
      = FixedGradientImageFilterType::New();

   fixedGradientFilter->SetInput( this->m_FixedImage );

   const typename FixedImageType::SpacingType & fixedSpacing =
      this->m_FixedImage->GetSpacing();
   double maximumFixedSpacing=0.0;
   for(unsigned int i=0; i<ImageDimension; i++)
   {
      if( fixedSpacing[i] > maximumFixedSpacing )
      {
         maximumFixedSpacing = fixedSpacing[i];
      }
   }
   fixedGradientFilter->SetSigma( maximumFixedSpacing );
   fixedGradientFilter->SetNormalizeAcrossScale( true );
   fixedGradientFilter->SetNumberOfThreads( m_NumberOfThreads );

#ifdef ITK_USE_ORIENTED_IMAGE_DIRECTION
   fixedGradientFilter->SetUseImageDirection( true );
#endif

   fixedGradientFilter->Update();

   this->m_FixedGradientImage = fixedGradientFilter->GetOutput();
   }
}

// Method to reinitialize the seed of the random number generator
template < class TFixedImage, class TMovingImage  > void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::ReinitializeSeed()
{
   Statistics::MersenneTwisterRandomVariateGenerator::GetInstance()->SetSeed();
}

// Method to reinitialize the seed of the random number generator
template < class TFixedImage, class TMovingImage  > void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::ReinitializeSeed(int seed)
{
   Statistics::MersenneTwisterRandomVariateGenerator::GetInstance()->SetSeed(
      seed);
}


/**
 * Transform a point from FixedImage domain to MovingImage domain.
 * This function also checks if mapped point is within support region.
 */
template < class TFixedImage, class TMovingImage >
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::TransformPoint( unsigned int sampleNumber,
                  PointType& mappedPoint,
                  bool& sampleOk,
                  double& movingImageValue,
                  unsigned int threadID ) const
{
   sampleOk = true;
   ESMTransformType * transform;

   if( threadID > 0 )
   {
      transform = this->m_ThreaderESMTransform[threadID-1];
   }
   else
   {
      transform = this->m_ESMTransform;
   }

   // Use generic transform to compute mapped position
   mappedPoint = transform->TransformPoint( this->m_FixedImageSamples[sampleNumber].point );
   sampleOk = true;

   if(sampleOk)
   {
      // If user provided a mask over the Moving image
      if ( this->m_MovingImageMask )
      {
         // Check if mapped point is within the support region of the moving image
         // mask
         sampleOk = sampleOk && this->m_MovingImageMask->IsInside( mappedPoint );
      }


      // Check if mapped point inside image buffer
      sampleOk = sampleOk && this->m_Interpolator->IsInsideBuffer( mappedPoint );
      if( sampleOk )
      {
         movingImageValue = this->m_Interpolator->Evaluate( mappedPoint );
      }
   }
}


/**
 * Transform a point from FixedImage domain to MovingImage domain.
 * This function also checks if mapped point is within support region.
 */
template < class TFixedImage, class TMovingImage >
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::TransformPointWithDerivatives( unsigned int sampleNumber,
                                 PointType& mappedPoint,
                                 bool& sampleOk,
                                 double& movingImageValue,
                                 ImageDerivativesType & movingImageGradient,
                                 ImageDerivativesType & fixedImageGradient,
                                 unsigned int threadID ) const
{
   ESMTransformType * transform;
   sampleOk = true;

   if( threadID > 0 )
   {
      transform = this->m_ThreaderESMTransform[threadID-1];
   }
   else
   {
      transform = this->m_ESMTransform;
   }

   // Use generic transform to compute mapped position
   mappedPoint = transform->TransformPoint( this->m_FixedImageSamples[sampleNumber].point );
   sampleOk = true;

   if(sampleOk)
   {
      // If user provided a mask over the Moving image
      if ( this->m_MovingImageMask )
      {
         // Check if mapped point is within the support region of the moving image
         // mask
         sampleOk = sampleOk && this->m_MovingImageMask->IsInside( mappedPoint );
      }

      // Check if mapped point inside image buffer
      sampleOk = sampleOk && this->m_Interpolator->IsInsideBuffer( mappedPoint );
      if( sampleOk )
      {
         this->ComputeImageDerivatives( mappedPoint, movingImageGradient );
         movingImageValue = this->m_Interpolator->Evaluate( mappedPoint );

         this->ComputeFixedImageDerivatives(
            this->m_FixedImageSamples[sampleNumber].point, fixedImageGradient );
      }
   }
}

/**
 * Compute image derivatives
 */
template < class TFixedImage, class TMovingImage >
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::ComputeImageDerivatives( const PointType & mappedPoint,
                           ImageDerivativesType & gradient ) const
{
   if ( this->m_ComputeGradient )
   {
      ContinuousIndex<double, ImageDimension> tempIndex;
      this->m_MovingImage->TransformPhysicalPointToContinuousIndex( mappedPoint,
                                                                    tempIndex );
      ImageIndexType mappedIndex;
      mappedIndex.CopyWithRound( tempIndex );
      gradient = this->m_GradientImage->GetPixel( mappedIndex );
   }
   else
   {
      // if not using the gradient image
      gradient = this->m_DerivativeCalculator->Evaluate( mappedPoint );
   }
}

/**
 * Compute fixed image derivatives
 */
template < class TFixedImage, class TMovingImage >
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::ComputeFixedImageDerivatives( const PointType & fixedPoint,
                                ImageDerivativesType & gradient ) const
{
   if ( this->m_ComputeGradient )
   {
      ContinuousIndex<double, ImageDimension> tempIndex;
      this->m_FixedImage->TransformPhysicalPointToContinuousIndex( fixedPoint,
                                                                   tempIndex );
      ImageIndexType fixedIndex;
      fixedIndex.CopyWithRound( tempIndex );
      gradient = this->m_FixedGradientImage->GetPixel( fixedIndex );
   }
   else
   {
      // if not using the gradient image
      gradient = this->m_FixedDerivativeCalculator->Evaluate( fixedPoint );
   }
}

template < class TFixedImage, class TMovingImage  >
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::GetValueMultiThreadedPreProcessInitiate( void ) const
{
   this->SynchronizeTransforms();

   m_Threader->SetSingleMethod(GetValueMultiThreadedPreProcess,
                               (void *)(&m_ThreaderParameter));
   m_Threader->SingleMethodExecute();
}

template < class TFixedImage, class TMovingImage  >
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::GetValueMultiThreadedInitiate( void ) const
{

   this->SynchronizeTransforms();

   m_Threader->SetSingleMethod(GetValueMultiThreaded,
                               (void *)(&m_ThreaderParameter));
   m_Threader->SingleMethodExecute();

   for( unsigned int threadID = 0; threadID<m_NumberOfThreads-1; threadID++ )
   {
      this->m_NumberOfPixelsCounted += m_ThreaderNumberOfMovingImageSamples[threadID];
   }
}

template < class TFixedImage, class TMovingImage  >
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::GetValueMultiThreadedPostProcessInitiate( void ) const
{
   m_Threader->SetSingleMethod(GetValueMultiThreadedPostProcess,
                               (void *)(&m_ThreaderParameter));
   m_Threader->SingleMethodExecute();
}

/**
 * Get the match Measure
 */
template < class TFixedImage, class TMovingImage  >
ITK_THREAD_RETURN_TYPE
ESMImageToImageMetric<TFixedImage,TMovingImage>
::GetValueMultiThreadedPreProcess( void * arg )
{
   int threadID;
   MultiThreaderParameterType * mtParam;

   threadID = ((MultiThreaderType::ThreadInfoStruct *)(arg))->ThreadID;

   mtParam = (MultiThreaderParameterType *)
      (((MultiThreaderType::ThreadInfoStruct *)(arg))->UserData);

   mtParam->metric->GetValueThreadPreProcess(threadID, false);

   return ITK_THREAD_RETURN_VALUE;
}


/**
 * Get the match Measure
 */
template < class TFixedImage, class TMovingImage  >
ITK_THREAD_RETURN_TYPE
ESMImageToImageMetric<TFixedImage,TMovingImage>
::GetValueMultiThreaded( void * arg )
{
   int threadID;
   MultiThreaderParameterType * mtParam;

   threadID = ((MultiThreaderType::ThreadInfoStruct *)(arg))->ThreadID;

   mtParam = (MultiThreaderParameterType *)
      (((MultiThreaderType::ThreadInfoStruct *)(arg))->UserData);

   mtParam->metric->GetValueThread(threadID);

   return ITK_THREAD_RETURN_VALUE;
}
/**
 * Get the match Measure
 */
template < class TFixedImage, class TMovingImage  >
ITK_THREAD_RETURN_TYPE
ESMImageToImageMetric<TFixedImage,TMovingImage>
::GetValueMultiThreadedPostProcess( void * arg )
{
   int threadID;
   MultiThreaderParameterType * mtParam;

   threadID = ((MultiThreaderType::ThreadInfoStruct *)(arg))->ThreadID;

   mtParam = (MultiThreaderParameterType *)
      (((MultiThreaderType::ThreadInfoStruct *)(arg))->UserData);

   mtParam->metric->GetValueThreadPostProcess(threadID, false);

   return ITK_THREAD_RETURN_VALUE;
}


template < class TFixedImage, class TMovingImage  >
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::GetValueThread( unsigned int threadID ) const
{
   // Figure out how many samples to process
   int chunkSize = m_NumberOfFixedImageSamples / m_NumberOfThreads;

   // Skip to this thread's samples to process
   unsigned int fixedImageSample = threadID * chunkSize;

   if(threadID == m_NumberOfThreads - 1)
   {
      chunkSize = m_NumberOfFixedImageSamples
         - ((m_NumberOfThreads-1)
            * chunkSize);
   }

   int numSamples = 0;

   if(m_WithinThreadPreProcess)
   {
      this->GetValueThreadPreProcess(threadID, true);
   }

   // Process the samples
   PointType mappedPoint;
   bool sampleOk;
   double movingImageValue;
   for( int count=0; count < chunkSize; ++count, ++fixedImageSample )
   {
      // Get moving image value
      this->TransformPoint( fixedImageSample, mappedPoint, sampleOk, movingImageValue,
                            threadID );

      if( sampleOk )
      {
         // CALL USER FUNCTION
         if(GetValueThreadProcessSample(threadID, fixedImageSample,
                                        mappedPoint, movingImageValue))
         {
            ++numSamples;
         }
      }
   }

   if(threadID > 0)
   {
      m_ThreaderNumberOfMovingImageSamples[threadID-1] = numSamples;
   }
   else
   {
      m_NumberOfPixelsCounted = numSamples;
   }

   if(m_WithinThreadPostProcess)
   {
      this->GetValueThreadPostProcess(threadID, true);
   }
}

template < class TFixedImage, class TMovingImage  >
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::GetValueAndDerivativeMultiThreadedPreProcessInitiate( void ) const
{
   this->SynchronizeTransforms();

   m_Threader->SetSingleMethod(GetValueAndDerivativeMultiThreadedPreProcess,
                               (void *)(&m_ThreaderParameter));
   m_Threader->SingleMethodExecute();
}

template < class TFixedImage, class TMovingImage  >
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::GetValueAndDerivativeMultiThreadedInitiate( void ) const
{
   this->SynchronizeTransforms();

   m_Threader->SetSingleMethod(GetValueAndDerivativeMultiThreaded,
                               (void *)(&m_ThreaderParameter));
   m_Threader->SingleMethodExecute();

   for( unsigned int threadID = 0; threadID<m_NumberOfThreads-1; threadID++ )
   {
      this->m_NumberOfPixelsCounted += m_ThreaderNumberOfMovingImageSamples[threadID];
   }
}

template < class TFixedImage, class TMovingImage  >
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::GetValueAndDerivativeMultiThreadedPostProcessInitiate( void ) const
{
   m_Threader->SetSingleMethod(GetValueAndDerivativeMultiThreadedPostProcess,
                               (void *)(&m_ThreaderParameter));
   m_Threader->SingleMethodExecute();
}

/**
 * Get the match Measure
 */
template < class TFixedImage, class TMovingImage  >
ITK_THREAD_RETURN_TYPE
ESMImageToImageMetric<TFixedImage,TMovingImage>
::GetValueAndDerivativeMultiThreadedPreProcess( void * arg )
{
   int threadID;
   MultiThreaderParameterType * mtParam;

   threadID = ((MultiThreaderType::ThreadInfoStruct *)(arg))->ThreadID;

   mtParam = (MultiThreaderParameterType *)
      (((MultiThreaderType::ThreadInfoStruct *)(arg))->UserData);

   mtParam->metric->GetValueAndDerivativeThreadPreProcess(threadID, false);

   return ITK_THREAD_RETURN_VALUE;
}

/**
 * Get the match Measure
 */
template < class TFixedImage, class TMovingImage  >
ITK_THREAD_RETURN_TYPE
ESMImageToImageMetric<TFixedImage,TMovingImage>
::GetValueAndDerivativeMultiThreaded( void * arg )
{
   int threadID;
   MultiThreaderParameterType * mtParam;

   threadID = ((MultiThreaderType::ThreadInfoStruct *)(arg))->ThreadID;

   mtParam = (MultiThreaderParameterType *)
      (((MultiThreaderType::ThreadInfoStruct *)(arg))->UserData);

   mtParam->metric->GetValueAndDerivativeThread(threadID);

   return ITK_THREAD_RETURN_VALUE;
}

/**
 * Get the match Measure
 */
template < class TFixedImage, class TMovingImage  >
ITK_THREAD_RETURN_TYPE
ESMImageToImageMetric<TFixedImage,TMovingImage>
::GetValueAndDerivativeMultiThreadedPostProcess( void * arg )
{
   int threadID;
   MultiThreaderParameterType * mtParam;

   threadID = ((MultiThreaderType::ThreadInfoStruct *)(arg))->ThreadID;

   mtParam = (MultiThreaderParameterType *)
      (((MultiThreaderType::ThreadInfoStruct *)(arg))->UserData);

   mtParam->metric->GetValueAndDerivativeThreadPostProcess(threadID, false);

   return ITK_THREAD_RETURN_VALUE;
}

template < class TFixedImage, class TMovingImage  >
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::GetValueAndDerivativeThread( unsigned int threadID ) const
{
   // Figure out how many samples to process
   int chunkSize = m_NumberOfFixedImageSamples / m_NumberOfThreads;

   // Skip to this thread's samples to process
   unsigned int fixedImageSample = threadID * chunkSize;

   if(threadID == m_NumberOfThreads - 1)
   {
      chunkSize = m_NumberOfFixedImageSamples
         - ((m_NumberOfThreads-1)
            * chunkSize);
   }

   int numSamples = 0;

   if(m_WithinThreadPreProcess)
   {
      this->GetValueAndDerivativeThreadPreProcess(threadID, true);
   }

   // Process the samples
   PointType mappedPoint;
   bool sampleOk;
   double movingImageValue;
   ImageDerivativesType movingImageGradientValue;
   ImageDerivativesType fixedImageGradientValue;
   for( int count=0; count < chunkSize; ++count, ++fixedImageSample )
   {
      // Get moving image value
      TransformPointWithDerivatives( fixedImageSample, mappedPoint, sampleOk,
                                     movingImageValue, movingImageGradientValue,
                                     fixedImageGradientValue, threadID );

      if( sampleOk )
      {
         // CALL USER FUNCTION
         if( this->GetValueAndDerivativeThreadProcessSample(
                threadID,
                fixedImageSample,
                mappedPoint,
                movingImageValue,
                movingImageGradientValue,
                fixedImageGradientValue ))
         {
            ++numSamples;
         }
      }
   }

   if(threadID > 0)
   {
      m_ThreaderNumberOfMovingImageSamples[threadID-1] = numSamples;
   }
   else
   {
      m_NumberOfPixelsCounted = numSamples;
   }

   if(m_WithinThreadPostProcess)
   {
      this->GetValueAndDerivativeThreadPostProcess(threadID, true);
   }
}


template < class TFixedImage, class TMovingImage  >
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::GetValueDerivativeAndHessianMultiThreadedPreProcessInitiate( void ) const
{
   this->SynchronizeTransforms();

   m_Threader->SetSingleMethod(GetValueDerivativeAndHessianMultiThreadedPreProcess,
                               (void *)(&m_ThreaderParameter));
   m_Threader->SingleMethodExecute();
}

template < class TFixedImage, class TMovingImage  >
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::GetValueDerivativeAndHessianMultiThreadedInitiate( void ) const
{
   this->SynchronizeTransforms();

   m_Threader->SetSingleMethod(GetValueDerivativeAndHessianMultiThreaded,
                               (void *)(&m_ThreaderParameter));
   m_Threader->SingleMethodExecute();

   for( unsigned int threadID = 0; threadID<m_NumberOfThreads-1; threadID++ )
   {
      this->m_NumberOfPixelsCounted += m_ThreaderNumberOfMovingImageSamples[threadID];
   }
}

template < class TFixedImage, class TMovingImage  >
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::GetValueDerivativeAndHessianMultiThreadedPostProcessInitiate( void ) const
{
   m_Threader->SetSingleMethod(GetValueDerivativeAndHessianMultiThreadedPostProcess,
                               (void *)(&m_ThreaderParameter));
   m_Threader->SingleMethodExecute();
}

/**
 * Get the match Measure
 */
template < class TFixedImage, class TMovingImage  >
ITK_THREAD_RETURN_TYPE
ESMImageToImageMetric<TFixedImage,TMovingImage>
::GetValueDerivativeAndHessianMultiThreadedPreProcess( void * arg )
{
   int threadID;
   MultiThreaderParameterType * mtParam;

   threadID = ((MultiThreaderType::ThreadInfoStruct *)(arg))->ThreadID;

   mtParam = (MultiThreaderParameterType *)
      (((MultiThreaderType::ThreadInfoStruct *)(arg))->UserData);

   mtParam->metric->GetValueDerivativeAndHessianThreadPreProcess(threadID, false);

   return ITK_THREAD_RETURN_VALUE;
}

/**
 * Get the match Measure
 */
template < class TFixedImage, class TMovingImage  >
ITK_THREAD_RETURN_TYPE
ESMImageToImageMetric<TFixedImage,TMovingImage>
::GetValueDerivativeAndHessianMultiThreaded( void * arg )
{
   int threadID;
   MultiThreaderParameterType * mtParam;

   threadID = ((MultiThreaderType::ThreadInfoStruct *)(arg))->ThreadID;

   mtParam = (MultiThreaderParameterType *)
      (((MultiThreaderType::ThreadInfoStruct *)(arg))->UserData);

   mtParam->metric->GetValueDerivativeAndHessianThread(threadID);

   return ITK_THREAD_RETURN_VALUE;
}

/**
 * Get the match Measure
 */
template < class TFixedImage, class TMovingImage  >
ITK_THREAD_RETURN_TYPE
ESMImageToImageMetric<TFixedImage,TMovingImage>
::GetValueDerivativeAndHessianMultiThreadedPostProcess( void * arg )
{
   int threadID;
   MultiThreaderParameterType * mtParam;

   threadID = ((MultiThreaderType::ThreadInfoStruct *)(arg))->ThreadID;

   mtParam = (MultiThreaderParameterType *)
      (((MultiThreaderType::ThreadInfoStruct *)(arg))->UserData);

   mtParam->metric->GetValueDerivativeAndHessianThreadPostProcess(threadID, false);

   return ITK_THREAD_RETURN_VALUE;
}

template < class TFixedImage, class TMovingImage  >
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::GetValueDerivativeAndHessianThread( unsigned int threadID ) const
{
   // Figure out how many samples to process
   int chunkSize = m_NumberOfFixedImageSamples / m_NumberOfThreads;

   // Skip to this thread's samples to process
   unsigned int fixedImageSample = threadID * chunkSize;

   if(threadID == m_NumberOfThreads - 1)
   {
      chunkSize = m_NumberOfFixedImageSamples
         - ((m_NumberOfThreads-1)
            * chunkSize);
   }

   int numSamples = 0;

   if(m_WithinThreadPreProcess)
   {
      this->GetValueDerivativeAndHessianThreadPreProcess(threadID, true);
   }

   // Process the samples
   PointType mappedPoint;
   bool sampleOk;
   double movingImageValue;
   ImageDerivativesType movingImageGradientValue;
   ImageDerivativesType fixedImageGradientValue;
   for( int count=0; count < chunkSize; ++count, ++fixedImageSample )
   {
      // Get moving image value
      TransformPointWithDerivatives( fixedImageSample, mappedPoint, sampleOk,
                                     movingImageValue, movingImageGradientValue,
                                     fixedImageGradientValue, threadID );

      if( sampleOk )
      {
         // CALL USER FUNCTION
         if( this->GetValueDerivativeAndHessianThreadProcessSample(
                threadID,
                fixedImageSample,
                mappedPoint,
                movingImageValue,
                movingImageGradientValue,
                fixedImageGradientValue ))
         {
            ++numSamples;
         }
      }
   }

   if(threadID > 0)
   {
      m_ThreaderNumberOfMovingImageSamples[threadID-1] = numSamples;
   }
   else
   {
      this->m_NumberOfPixelsCounted = numSamples;
   }

   if(m_WithinThreadPostProcess)
   {
      this->GetValueDerivativeAndHessianThreadPostProcess(threadID, true);
   }
}


/**
 * PrintSelf
 */
template <class TFixedImage, class TMovingImage>
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
   Superclass::PrintSelf( os, indent );

   os << indent << "NumberOfFixedImageSamples: ";
   os << m_NumberOfFixedImageSamples << std::endl;

   if( m_UseFixedImageIndexes )
   {
      os << indent << "Use Fixed Image Indexes: True" << std::endl;
      os << indent << "Number of Fixed Image Indexes = "
         << m_FixedImageIndexes.size() << std::endl;
   }
   else
   {
      os << indent << "Use Fixed Image Indexes: False" << std::endl;
   }

   if( m_UseSequentialSampling )
   {
      os << indent << "Use Sequential Sampling: True" << std::endl;
   }
   else
   {
      os << indent << "Use Sequential Sampling: False" << std::endl;
   }

   os << indent << "UseAllPixels: ";
   os << m_UseAllPixels << std::endl;

   os << indent << "Threader: " << m_Threader << std::endl;
   os << indent << "Number of Threads: " << m_NumberOfThreads << std::endl;
   os << indent << "ThreaderParameter: " << std::endl;
   os << indent << "ThreaderNumberOfMovingImageSamples: " << std::endl;
   if( m_ThreaderNumberOfMovingImageSamples )
   {
      for(unsigned int i=0; i<m_NumberOfThreads-1; i++)
      {
         os << "  Thread[" << i << "]= " << (unsigned int)m_ThreaderNumberOfMovingImageSamples[i] << std::endl;
      }
   }

   os << indent << "Transform:    " << m_ESMTransform.GetPointer()    << std::endl;
   os << indent << "Number of Moving Image Samples: " << m_NumberOfPixelsCounted
      << std::endl;
}

/** This method can be const because we are not altering the m_ThreaderTransform
 *  pointer. We are altering the object that m_ThreaderTransform[idx] points at.
 *  This is allowed under C++ const rules.
 */
template <class TFixedImage, class TMovingImage>
void
ESMImageToImageMetric<TFixedImage,TMovingImage>
::SynchronizeTransforms() const
{
   for( unsigned int threadID = 0; threadID<m_NumberOfThreads-1; threadID++ )
   {
      /** Set the fixed parameters first. Some transforms have parameters which depend on
          the values of the fixed parameters. For instance, the BSplineDeformableTransform
          checks the grid size (part of the fixed parameters) before setting the parameters. */
      this->m_ThreaderESMTransform[threadID]->SetFixedParameters( this->m_ESMTransform->GetFixedParameters() );
      this->m_ThreaderESMTransform[threadID]->SetParameters( this->m_ESMTransform->GetParameters() );
   }
}

} // end namespace itk

#endif

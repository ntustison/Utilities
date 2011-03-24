#ifndef MZ_ESMMeanSquaresImageToImageMetric_TXX_
#define MZ_ESMMeanSquaresImageToImageMetric_TXX_

#include "itkESMMeanSquaresImageToImageMetric.h"

#include <itkCovariantVector.h>
#include <itkImageRandomConstIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkImageIterator.h>
#include <vnl/vnl_math.h>

namespace itk
{

/**
 * Constructor
 */
template < class TFixedImage, class TMovingImage >
ESMMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::ESMMeanSquaresImageToImageMetric()
   :Superclass()
{
   this->SetComputeGradient(true);

   m_ThreaderMSE = NULL;
   m_ThreaderMSEDerivatives = NULL;
   m_ThreaderMSEHessians = NULL;
   this->m_WithinThreadPreProcess = false;
   this->m_WithinThreadPostProcess = false;

   //  For backward compatibility, the default behavior is to use all the pixels
   //  in the fixed image.
   this->SetUseAllPixels( true );
}

template < class TFixedImage, class TMovingImage >
ESMMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::~ESMMeanSquaresImageToImageMetric()
{
   if(m_ThreaderMSE != NULL)
   {
      delete [] m_ThreaderMSE;
   }
   m_ThreaderMSE = NULL;

   if(m_ThreaderMSEDerivatives != NULL)
   {
      delete [] m_ThreaderMSEDerivatives;
   }
   m_ThreaderMSEDerivatives = NULL;

   if(m_ThreaderMSEHessians != NULL)
   {
      delete [] m_ThreaderMSEHessians;
   }
   m_ThreaderMSEHessians = NULL;
}

/**
 * Print out internal information about this class
 */
template < class TFixedImage, class TMovingImage  >
void
ESMMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::PrintSelf(std::ostream& os, Indent indent) const
{

   Superclass::PrintSelf(os, indent);

}


/**
 * Initialize
 */
template <class TFixedImage, class TMovingImage>
void
ESMMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::Initialize(void) throw ( ExceptionObject )
{

   this->Superclass::Initialize();
   this->Superclass::MultiThreadingInitialize();

   if(m_ThreaderMSE != NULL)
   {
      delete [] m_ThreaderMSE;
   }
   m_ThreaderMSE = new double[this->m_NumberOfThreads];

   if(m_ThreaderMSEDerivatives != NULL)
   {
      delete [] m_ThreaderMSEDerivatives;
   }
   m_ThreaderMSEDerivatives = new DerivativeType[this->m_NumberOfThreads];
   for(unsigned int threadID=0; threadID<this->m_NumberOfThreads; ++threadID)
   {
      m_ThreaderMSEDerivatives[threadID].SetSize( this->m_NumberOfParameters );
   }

   if(m_ThreaderMSEHessians != NULL)
   {
      delete [] m_ThreaderMSEHessians;
   }
   m_ThreaderMSEHessians = new HessianType[this->m_NumberOfThreads];
   for(unsigned int threadID=0; threadID<this->m_NumberOfThreads; ++threadID)
   {
      m_ThreaderMSEHessians[threadID].set_size( this->m_NumberOfParameters, this->m_NumberOfParameters );
   }
}

template < class TFixedImage, class TMovingImage  >
inline bool
ESMMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::GetValueThreadProcessSample(
   unsigned int threadID,
   unsigned long fixedImageSample,
   const PointType & itkNotUsed(mappedPoint),
   double movingImageValue) const
{
   double diff = movingImageValue - this->m_FixedImageSamples[fixedImageSample].value;

   m_ThreaderMSE[threadID] += diff*diff;

   return true;
}

template < class TFixedImage, class TMovingImage  >
typename ESMMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::MeasureType
ESMMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::GetValue( const ParametersType & parameters ) const
{
   if( !this->m_FixedImage )
   {
      itkExceptionMacro( << "Fixed image has not been assigned" );
   }

   // Set output measures values to zero
   memset( m_ThreaderMSE,
           0,
           this->m_NumberOfThreads * sizeof(MeasureType) );

   // Set up the parameters in the transform
   this->SetTransformParameters( parameters );


   // MUST BE CALLED TO INITIATE PROCESSING
   this->GetValueMultiThreadedInitiate();

   itkDebugMacro( "Ratio of voxels mapping into moving image buffer: "
                  << this->m_NumberOfPixelsCounted << " / "
                  << this->m_NumberOfFixedImageSamples
                  << std::endl );

   if( this->m_NumberOfPixelsCounted <
       this->m_NumberOfFixedImageSamples / 4 )
   {
      /*itkExceptionMacro( "Too many samples map outside moving image buffer: "
                         << this->m_NumberOfPixelsCounted << " / "
                         << this->m_NumberOfFixedImageSamples
                         << std::endl );*/
      return std::numeric_limits<MeasureType>::max();
   }

   double mse = m_ThreaderMSE[0];
   for(unsigned int t=1; t<this->m_NumberOfThreads; ++t)
   {
      mse += m_ThreaderMSE[t];
   }
   mse /= this->m_NumberOfPixelsCounted;

   return mse;
}


template < class TFixedImage, class TMovingImage  >
inline bool
ESMMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::GetValueAndDerivativeThreadProcessSample(
   unsigned int threadID,
   unsigned long fixedImageSample,
   const PointType & itkNotUsed(mappedPoint),
   double movingImageValue,
   const ImageDerivativesType & movingImageGradientValue,
   const ImageDerivativesType & fixedImageGradientValue ) const
{
   const double minusDiff = movingImageValue
      - this->m_FixedImageSamples[fixedImageSample].value;

   m_ThreaderMSE[threadID] += minusDiff*minusDiff;

   const PointType & fixedImagePoint = this->m_FixedImageSamples[fixedImageSample].point;

   // Need to use one of the threader transforms if we're
   // not in thread 0.
   //
   // Use a raw pointer here to avoid the overhead of smart pointers.
   // For instance, Register and UnRegister have mutex locks around
   // the reference counts.
   ESMTransformType* transform;

   if (threadID > 0)
   {
      transform = this->m_ThreaderESMTransform[threadID - 1];
   }
   else
   {
      transform = this->m_ESMTransform;
   }

   // Jacobian is evaluated at the unmapped (fixed image) point.
   const TransformJacobianType & spatial_jac = transform->GetSpatialJacobian( fixedImagePoint );
   const TransformJacobianType & increment_jac = transform->GetIncrementalUpdateJacobian( fixedImagePoint );

   ///\todo replace use of spatial Jacobian by prior resampling of the moving image

   m_ThreaderMSEDerivatives[threadID] +=
      ( (movingImageGradientValue.GetVnlVector()*spatial_jac + fixedImageGradientValue.GetVnlVector())
        * minusDiff ) * increment_jac;

   //(movingImageGradientValue*minusDiff ).GetVnlVector() * (spatial_jac*increment_jac);

   return true;
}

template < class TFixedImage, class TMovingImage  >
inline bool
ESMMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::GetValueDerivativeAndHessianThreadProcessSample(
   unsigned int threadID,
   unsigned long fixedImageSample,
   const PointType & itkNotUsed(mappedPoint),
   double movingImageValue,
   const ImageDerivativesType & movingImageGradientValue,
   const ImageDerivativesType & fixedImageGradientValue ) const
{
   const double minusDiff = movingImageValue - this->m_FixedImageSamples[fixedImageSample].value;

   m_ThreaderMSE[threadID] += minusDiff*minusDiff;

   const PointType & fixedImagePoint = this->m_FixedImageSamples[fixedImageSample].point;

   // Need to use one of the threader transforms if we're
   // not in thread 0.
   //
   // Use a raw pointer here to avoid the overhead of smart pointers.
   // For instance, Register and UnRegister have mutex locks around
   // the reference counts.
   ESMTransformType* transform;

   if (threadID > 0)
   {
      transform = this->m_ThreaderESMTransform[threadID - 1];
   }
   else
   {
      transform = this->m_ESMTransform;
   }

   // Jacobian is evaluated at the unmapped (fixed image) point.
   const TransformJacobianType & spatial_jac = transform->GetSpatialJacobian( fixedImagePoint );
   const TransformJacobianType & increment_jac = transform->GetIncrementalUpdateJacobian( fixedImagePoint );

   ///\todo replace use of spatial Jacobian by prior resampling of the moving image

   //const vnl_vector<double> minusJ( movingImageGradientValue.GetVnlVector()*(spatial_jac*increment_jac) );
   const vnl_vector<double> minusJtimes2(
      ( movingImageGradientValue.GetVnlVector()*spatial_jac
        + fixedImageGradientValue.GetVnlVector() ) * increment_jac );

   //m_ThreaderMSEDerivatives[threadID] += minusJ * minusDiff;
   m_ThreaderMSEDerivatives[threadID] += minusJtimes2 * minusDiff;

   // The following code is equivalent to
   // m_ThreaderMSEHessians[threadID] += outer_product( minusJ, minusJ );
   // but uses the fact that it is twice the same vector
   for(unsigned int i=0; i<this->m_NumberOfParameters; ++i)
   {
      //m_ThreaderMSEHessians[threadID](i,i) += vnl_math_sqr(minusJ[i]);
      m_ThreaderMSEHessians[threadID](i,i) += vnl_math_sqr(minusJtimes2[i]);
      for(unsigned int j=i+1; j<this->m_NumberOfParameters; ++j)
      {
         //const double tmp = minusJ[i]*minusJ[j];
         const double tmp = minusJtimes2[i]*minusJtimes2[j];
         m_ThreaderMSEHessians[threadID](i,j) += tmp;
         m_ThreaderMSEHessians[threadID](j,i) += tmp;
      }
   }

   return true;
}

/**
 * Get the both Value and Derivative Measure
 */
template < class TFixedImage, class TMovingImage  >
void
ESMMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::GetValueAndDerivative( const ParametersType & parameters,
                         MeasureType & value,
                         DerivativeType & derivative) const
{

   if( !this->m_FixedImage )
   {
      itkExceptionMacro( << "Fixed image has not been assigned" );
   }

   // Set up the parameters in the transform
   this->SetTransformParameters( parameters );

   // Set output values to zero
   memset( m_ThreaderMSE,
           0,
           this->m_NumberOfThreads * sizeof(MeasureType) );

   // Set output derivatives values to zero
   if(derivative.GetSize() != this->m_NumberOfParameters)
   {
      derivative = DerivativeType( this->m_NumberOfParameters );
   }
   memset( derivative.data_block(),
           0,
           this->m_NumberOfParameters * sizeof(double) );

   for( unsigned int threadID = 0; threadID<this->m_NumberOfThreads; threadID++ )
   {
      memset( m_ThreaderMSEDerivatives[threadID].data_block(),
              0,
              this->m_NumberOfParameters * sizeof(double) );
   }

   // MUST BE CALLED TO INITIATE PROCESSING
   this->GetValueAndDerivativeMultiThreadedInitiate();

   itkDebugMacro( "Ratio of voxels mapping into moving image buffer: "
                  << this->m_NumberOfPixelsCounted << " / "
                  << this->m_NumberOfFixedImageSamples
                  << std::endl );

   if( this->m_NumberOfPixelsCounted <
       this->m_NumberOfFixedImageSamples / 4 )
   {
      /*itkExceptionMacro( "Too many samples map outside moving image buffer: "
                         << this->m_NumberOfPixelsCounted << " / "
                         << this->m_NumberOfFixedImageSamples
                         << std::endl );*/
      value = std::numeric_limits<MeasureType>::max();
      return;
   }

   value = 0;
   for(unsigned int t=0; t<this->m_NumberOfThreads; ++t)
   {
      value += m_ThreaderMSE[t];
      derivative += m_ThreaderMSEDerivatives[t];
   }

   value /= this->m_NumberOfPixelsCounted;
   //derivative *= 2.0/this->m_NumberOfPixelsCounted;
   derivative *= 1.0/this->m_NumberOfPixelsCounted;
}


template < class TFixedImage, class TMovingImage  >
void
ESMMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::GetValueDerivativeAndHessian( const ParametersType & parameters,
                                MeasureType & value,
                                DerivativeType & derivative,
                                HessianType & hessian ) const
{
   if( !this->m_FixedImage )
   {
      itkExceptionMacro( << "Fixed image has not been assigned" );
   }

   // Set up the parameters in the transform
   this->SetTransformParameters( parameters );

   // Set output measures values to zero
   value = 0;
   memset( m_ThreaderMSE,
           0,
           this->m_NumberOfThreads * sizeof(MeasureType) );

   // Set output derivatives values to zero
   if(derivative.GetSize() != this->m_NumberOfParameters)
   {
      derivative = DerivativeType( this->m_NumberOfParameters );
   }
   memset( derivative.data_block(),
           0,
           this->m_NumberOfParameters * sizeof(double) );

   for( unsigned int threadID = 0; threadID<this->m_NumberOfThreads; ++threadID )
   {
      memset( m_ThreaderMSEDerivatives[threadID].data_block(),
              0,
              this->m_NumberOfParameters * sizeof(double) );
   }

   // Set output Hessians values to zero
   if ( (hessian.rows() != this->m_NumberOfParameters) ||
        (hessian.columns() != this->m_NumberOfParameters) )
   {
      hessian = HessianType( this->m_NumberOfParameters, this->m_NumberOfParameters );
   }
   memset( hessian.data_block(),
           0,
           this->m_NumberOfParameters * this->m_NumberOfParameters * sizeof(double) );

   for( unsigned int threadID = 0; threadID<this->m_NumberOfThreads; ++threadID )
   {
      memset( m_ThreaderMSEHessians[threadID].data_block(),
              0,
              this->m_NumberOfParameters * this->m_NumberOfParameters * sizeof(double) );
   }


   // MUST BE CALLED TO INITIATE PROCESSING
   this->GetValueDerivativeAndHessianMultiThreadedInitiate();

   itkDebugMacro( "Ratio of voxels mapping into moving image buffer: "
                  << this->m_NumberOfPixelsCounted << " / "
                  << this->m_NumberOfFixedImageSamples
                  << std::endl );

   if( this->m_NumberOfPixelsCounted <
       this->m_NumberOfFixedImageSamples / 4 )
   {
      /*itkExceptionMacro( "Too many samples map outside moving image buffer: "
                         << this->m_NumberOfPixelsCounted << " / "
                         << this->m_NumberOfFixedImageSamples
                         << std::endl );*/
      value = std::numeric_limits<MeasureType>::max();
      return;
   }

   // Merge the threaded results
   for(unsigned int t=0; t<this->m_NumberOfThreads; ++t)
   {
      value += m_ThreaderMSE[t];
      derivative += m_ThreaderMSEDerivatives[t];
      hessian += m_ThreaderMSEHessians[t];
   }

   value /= this->m_NumberOfPixelsCounted;
   //derivative *= 2.0/this->m_NumberOfPixelsCounted;
   //hessian *= 2.0/this->m_NumberOfPixelsCounted;
   derivative *= 1.0/this->m_NumberOfPixelsCounted;
   hessian *= 0.5/this->m_NumberOfPixelsCounted;
}


/**
 * Get the match measure derivative
 */
template < class TFixedImage, class TMovingImage  >
void
ESMMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::GetDerivative( const ParametersType & parameters,
                 DerivativeType & derivative ) const
{
   if( !this->m_FixedImage )
   {
      itkExceptionMacro( << "Fixed image has not been assigned" );
   }

   MeasureType value;
   // call the combined version
   this->GetValueAndDerivative( parameters, value, derivative );
}

} // end namespace itk


#endif

#ifndef MZ_ESMMeanSquaresImageToImageMetric_H_
#define MZ_ESMMeanSquaresImageToImageMetric_H_

#include "itkESMImageToImageMetric.h"
#include "itkCovariantVector.h"
#include "itkPoint.h"
#include "itkIndex.h"

#include "itkMultiThreader.h"

namespace itk
{

template <class TFixedImage,class TMovingImage >
class ITK_EXPORT ESMMeanSquaresImageToImageMetric :
    public ESMImageToImageMetric< TFixedImage, TMovingImage >
{
public:

  /** Standard class typedefs. */
  typedef ESMMeanSquaresImageToImageMetric                   Self;
  typedef ESMImageToImageMetric< TFixedImage, TMovingImage > Superclass;
  typedef SmartPointer<Self>                                 Pointer;
  typedef SmartPointer<const Self>                           ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ESMMeanSquaresImageToImageMetric, ESMImageToImageMetric);

  /** Types inherited from Superclass. */
  typedef typename Superclass::ESMTransformType         ESMTransformType;
  typedef typename Superclass::ESMTransformPointer      ESMTransformPointer;
  typedef typename Superclass::TransformJacobianType    TransformJacobianType;
  typedef typename Superclass::InterpolatorType         InterpolatorType;
  typedef typename Superclass::MeasureType              MeasureType;
  typedef typename Superclass::DerivativeType           DerivativeType;
  typedef typename Superclass::HessianType              HessianType;
  typedef typename Superclass::ParametersType           ParametersType;
  typedef typename Superclass::FixedImageType           FixedImageType;
  typedef typename Superclass::MovingImageType          MovingImageType;
  typedef typename Superclass::PointType                PointType;
  typedef typename Superclass::FixedImageConstPointer   FixedImageConstPointer;
  typedef typename Superclass::MovingImageConstPointer  MovingImageConstPointer;
  typedef typename Superclass::CoordinateRepresentationType
                                                   CoordinateRepresentationType;
  typedef typename Superclass::FixedImageSampleContainer
                                                        FixedImageSampleContainer;
  typedef typename Superclass::ImageDerivativesType     ImageDerivativesType;

  /** The image dimension. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       Superclass::ImageDimension );

  /**
   *  Initialize the Metric by
   *  (1) making sure that all the components are present and plugged
   *      together correctly,
   *  (2) uniformly select NumberOfSpatialSamples within
   *      the FixedImageRegion, and
   *  (3) allocate memory for pdf data structures. */
  virtual void Initialize(void) throw ( ExceptionObject );

  /**  Get the value. */
  MeasureType GetValue( const ParametersType & parameters ) const;

  /** Get the derivatives of the match measure. */
  void GetDerivative( const ParametersType & parameters,
                      DerivativeType & Derivative ) const;

  /**  Get the value and derivatives for single valued optimizers. */
  void GetValueAndDerivative( const ParametersType & parameters,
                              MeasureType & value,
                              DerivativeType & derivative ) const;

  /**  Get the value and derivatives for single valued optimizers. */
  void GetValueDerivativeAndHessian( const ParametersType & parameters,
                                     MeasureType & value,
                                     DerivativeType & derivative,
                                     HessianType & hessian ) const;

protected:

  ESMMeanSquaresImageToImageMetric();
  virtual ~ESMMeanSquaresImageToImageMetric();
  void PrintSelf(std::ostream& os, Indent indent) const;

private:

  //purposely not implemented
  ESMMeanSquaresImageToImageMetric(const Self &);
  //purposely not implemented
  void operator=(const Self &);

  inline bool GetValueThreadProcessSample( unsigned int threadID,
                                       unsigned long fixedImageSample,
                                       const PointType & mappedPoint,
                                       double movingImageValue ) const;

#ifdef ITK_USE_OPTIMIZED_REGISTRATION_METHODS
   using Superclass::GetValueAndDerivativeThreadProcessSample;
#endif
  inline bool GetValueAndDerivativeThreadProcessSample( unsigned int threadID,
                                       unsigned long fixedImageSample,
                                       const PointType & mappedPoint,
                                       double movingImageValue,
                                       const ImageDerivativesType & movingImageGradientValue,
                                       const ImageDerivativesType & fixedImageGradientValue ) const;

  inline bool GetValueDerivativeAndHessianThreadProcessSample(
     unsigned int threadID,
     unsigned long fixedImageSample,
     const PointType & mappedPoint,
     double movingImageValue,
     const ImageDerivativesType & movingImageGradientValue,
     const ImageDerivativesType & fixedImageGradientValue ) const;

  MeasureType    * m_ThreaderMSE;
  DerivativeType * m_ThreaderMSEDerivatives;
  HessianType    * m_ThreaderMSEHessians;

};

} // end namespace itk

#include "itkESMMeanSquaresImageToImageMetric.txx"

#endif

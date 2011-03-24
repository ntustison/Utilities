#ifndef __itkSmoothRepresentationFunction_h
#define __itkSmoothRepresentationFunction_h

#include "itkImage.h"
#include "itkImageFunction.h"

namespace itk
{

/** \class SmoothRepresentationFunction
 * \brief Computes a smooth function representation of the image.
 *
 * This class provides a method for computing a C^\inf smooth
 * function that interpolates over the image domain.  
 * 
 * \cite M. Siqueira, "Mesh Generation From Imaging Data", Ph.D. Thesis,
 * University of Pennsylvania, February 2006.
 *
 * \todo Only works for isotropic pixels.  Possible extension to anisotropic.
 * Since the function is C^infinity continuous, need a way to return any
 * user-specified derivative
 */
template <class TInputImage, class TOutput = double, class TCoordRep = double>
class ITK_EXPORT SmoothRepresentationFunction 
: public ImageFunction<TInputImage, TOutput, TCoordRep>
{
public:
  /** Standard class typedefs. */
  typedef SmoothRepresentationFunction                   Self;
  typedef ImageFunction<TInputImage, TOutput, TCoordRep> Superclass;
  typedef SmartPointer<Self>                             Pointer;
  typedef SmartPointer<const Self>                       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SmoothRepresentationFunction, 
               ImageFunction);

  /** Extract the dimension of the image. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TInputImage::ImageDimension );

  typedef double                                         RealType;
  typedef Image<RealType, ImageDimension>                RealImageType;
  
  /** Standard image types within this class. */
  typedef TInputImage                                    ImageType;
  typedef typename Superclass::IndexType                 IndexType;
  typedef typename Superclass::ContinuousIndexType       ContinuousIndexType;     
  typedef typename Superclass::PointType                 PointType;
  
  /** Set the input image.
   * \warning this method caches BufferedRegion information.
   * If the BufferedRegion has changed, user must call
   * SetInputImage again to update cached values. */
  virtual void SetInputImage( const ImageType * );

  /** Evaluate the function at specified Point position. */
  virtual TOutput Evaluate( const PointType & ) const;

  /** Evaluate the function at specified Index position. */
  virtual TOutput EvaluateAtIndex( const IndexType & ) const;

  /** Evaluate the function at specified ContinousIndex position. */
  virtual TOutput EvaluateAtContinuousIndex( const ContinuousIndexType & ) const;

  /** Evaluate the gradient at specified Point position. */
  PointType EvaluateGradient( const PointType & ) const;

  /** Evaluate the gradient at specified Index position. */
  PointType EvaluateGradientAtIndex( const IndexType & ) const;

  /** Evaluate the gradient at specified ContinousIndex position. */
  PointType EvaluateGradientAtContinuousIndex( const ContinuousIndexType & ) const;

protected:
  SmoothRepresentationFunction();
  virtual ~SmoothRepresentationFunction();
  void PrintSelf( std::ostream& os, Indent indent ) const;

private:
  SmoothRepresentationFunction( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  PointType FromGlobalToLocal( const PointType &, const IndexType & ) const;
  RealType LinearInterpolationFunction( const PointType &, const IndexType & ) const;
  PointType LinearInterpolationGradient( const PointType &, const IndexType & ) const;

  RealType HFunction( const RealType & ) const;
  RealType GFunction( const RealType & ) const;
  RealType WeightFunction( const PointType & ) const;

  RealType HFunctionDerivative( const RealType & ) const;
  RealType GFunctionDerivative( const RealType & ) const;
  PointType WeightFunctionDerivative( const PointType & ) const;

  RealType                                  m_d;
};  // class SmoothRepresentationFunction

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSmoothRepresentationFunction.txx"
#endif

#endif /* __itkSmoothRepresentationFunction_h */

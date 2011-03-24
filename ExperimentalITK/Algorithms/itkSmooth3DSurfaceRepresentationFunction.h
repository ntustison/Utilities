#ifndef __itkSmooth3DSurfaceRepresentationFunction_h
#define __itkSmooth3DSurfaceRepresentationFunction_h

#include "itkImageFunction.h"

namespace itk
{

/** \class Smooth3DSurfaceRepresentationFunction
 * \brief Computes a smooth function representation of the 
 * foreground/background boundary in a binary image.
 *
 * This class provides a method for computing an implicit smooth
 * function which approximates the discrete structure of a 
 * binary image.  The function is constructed such that the n-1
 * dimensional manifold defining the boundary of the binary image
 * is the zero level set of f^{-1}(0) of an implicit function.  
 * 
 * There are two requirements for constructing the representing
 * function, viz.  1) The binary image must be "well-composed" and 
 * 2) the Euclidean distance transform must be available. Both
 * requirements are satisfied by calls to the corresponding filters
 * in the ITK filter library (insight/Code/BasicFunctions) called:
 *    1) itkWellComposedImageFunction
 *    2) itkEuclideanDistanceTransformFunction
 *
 * The function is topology preserving in that the implicit function
 * maintains the topology of the well-composed binary image. 
 *
 * \cite M. Siqueira, 
 *
 * \ingroup Operators ?
 *
 * \todo Only works for isotropic pixels.  Possible extension to anisotropic.
 * Since the function is C^infinity continuous, need a way to return any
 * user-specified derivative
 */
template <class TInputImage, class TOutput = double, class TCoordRep = double>
class ITK_EXPORT Smooth3DSurfaceRepresentationFunction 
: public ImageFunction<TInputImage, TOutput, TCoordRep>
{
public:
  /** Standard class typedefs. */
  typedef Smooth3DSurfaceRepresentationFunction          Self;
  typedef ImageFunction<TInputImage, TOutput, TCoordRep> Superclass;
  typedef SmartPointer<Self>                             Pointer;
  typedef SmartPointer<const Self>                       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Smooth3DSurfaceRepresentationFunction, 
               ImageFunction);

  /** Extract the dimension of the image. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

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
  virtual void SetInputImage(const ImageType *);

  /** Evaluate the function at specified Point position. */
  virtual TOutput Evaluate(const PointType &) const;

  /** Evaluate the function at specified Index position. */
  virtual TOutput EvaluateAtIndex(const IndexType &idx) const
    {
      PointType pt;
      this->m_Image->TransformIndexToPhysicalPoint(idx, pt);
      return this->Evaluate(pt);
    }

  /** Evaluate the function at specified ContinousIndex position. */
  virtual TOutput EvaluateAtContinuousIndex(const ContinuousIndexType &cidx) const
    {
      PointType pt;
      this->m_Image->TransformContinuousIndexToPhysicalPoint(cidx, pt);
      return this->Evaluate(pt);
    }

  /** Evaluate the gradient at specified Point position. */
  PointType EvaluateGradient(const PointType &) const;

  /** Evaluate the gradient at specified Index position. */
  PointType EvaluateGradientAtIndex(const IndexType &idx) const
    {
      PointType pt;
      this->m_Image->TransformIndexToPhysicalPoint(idx, pt);
      return this->EvaluateGradient(pt);
    }

  /** Evaluate the gradient at specified ContinousIndex position. */
  PointType EvaluateGradientAtContinuousIndex(const ContinuousIndexType &cidx) const
    {
      PointType pt;
      this->m_Image->TransformContinuousIndexToPhysicalPoint(cidx, pt);
      return this->EvaluateGradient(pt);
    }

protected:
  Smooth3DSurfaceRepresentationFunction() {};
  virtual ~Smooth3DSurfaceRepresentationFunction() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  Smooth3DSurfaceRepresentationFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  PointType FromGlobalToLocal(const PointType &, const IndexType &) const;
  RealType FunctionLinearInterpolation(const PointType &, const IndexType &) const;
  PointType GradientLinearInterpolation(const PointType &, const IndexType &) const;

  RealType HFunction(const RealType &) const;
  RealType GFunction(const RealType &) const;
  RealType WeightFunction(const PointType &) const;

  RealType HFunctionDerivative(const RealType &) const;
  RealType GFunctionDerivative(const RealType &) const;
  PointType WeightFunctionDerivative(const PointType &) const;

  RealType                                  m_d;

};  // class Smooth3DSurfaceRepresentationFunction

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSmooth3DSurfaceRepresentationFunction.txx"
#endif

#endif /* __itkSmooth3DSurfaceRepresentationFunction_h */

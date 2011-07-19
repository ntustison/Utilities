#ifndef MZ_ESMOptimizer_H_
#define MZ_ESMOptimizer_H_

#include "itkESMOptimizerBase.h"


namespace itk
{

/** \class ESMOptimizer
 * \brief Implements ESM optimization
 *
 * \ingroup Numerics Optimizers
 */
template <class TFixedImage,  class TMovingImage>
class ITK_EXPORT ESMOptimizer :
      public ESMOptimizerBase<TFixedImage,TMovingImage>
{
public:
  /** Standard class typedefs. */
  typedef ESMOptimizer                                Self;
  typedef ESMOptimizerBase<TFixedImage,TMovingImage>  Superclass;
  typedef SmartPointer<Self>                          Pointer;
  typedef SmartPointer<const Self>                    ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( ESMOptimizer, ESMOptimizerBase );

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /**  Parameters type.
   *  It defines a position in the optimization search space. */
  typedef typename Superclass::ParametersType         ParametersType;

  /** Images types */
  typedef typename Superclass::MovingImageType        MovingImageType;
  typedef typename Superclass::FixedImageType         FixedImageType;

  /** Type of the Cost Function   */
  typedef typename Superclass::ESMCostFunctionType    ESMCostFunctionType;
  typedef typename Superclass::ESMCostFunctionPointer ESMCostFunctionPointer;

  /**  Measure type.
   *  It defines a type used to return the cost function value.  */
  typedef typename ESMCostFunctionType::MeasureType    MeasureType;

  /**  Derivative type.
   *  It defines a type used to return the cost function derivative. */
  typedef typename ESMCostFunctionType::DerivativeType DerivativeType;

  /**  Hessian type.
   *  It defines a type used to return the cost function derivative. */
  typedef typename ESMCostFunctionType::HessianType    HessianType;

protected:
  ESMOptimizer();
  virtual ~ESMOptimizer(){};

  virtual void GetUpdate();

  /** Print out internal state */
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  ESMOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk

#include "itkESMOptimizer.hxx"

#endif

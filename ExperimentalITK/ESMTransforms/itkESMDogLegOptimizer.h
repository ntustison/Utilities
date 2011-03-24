#ifndef MZ_ESMDogLegOptimizer_H_
#define MZ_ESMDogLegOptimizer_H_

#include "itkESMOptimizerBase.h"


namespace itk
{

/** \class ESMDogLegOptimizer
 * \brief Implements ESM optimization
 *
 * \ingroup Numerics Optimizers
 */
template <class TFixedImage,  class TMovingImage>
class ITK_EXPORT ESMDogLegOptimizer :
      public ESMOptimizerBase<TFixedImage,TMovingImage>
{
public:
  /** Standard class typedefs. */
  typedef ESMDogLegOptimizer                          Self;
  typedef ESMOptimizerBase<TFixedImage,TMovingImage>  Superclass;
  typedef SmartPointer<Self>                          Pointer;
  typedef SmartPointer<const Self>                    ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( ESMDogLegOptimizer, ESMOptimizerBase );

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

  itkSetMacro( InitialTrustRegionRadius, double );
  itkGetConstReferenceMacro( InitialTrustRegionRadius, double );

protected:
  ESMDogLegOptimizer();
  virtual ~ESMDogLegOptimizer(){};

  /** Print out internal state */
  void PrintSelf(std::ostream& os, Indent indent) const;

  virtual void InitOptimization();
  virtual void GetUpdate();
  virtual void ApplyUpdate(){}; /// does nothing

  typedef enum {
     GaussNewtonStep,
     ScaledToTrustRegionSteepestDescentStep,
     TrueDogLegStep
  } DogLegStepType;

  virtual DogLegStepType GetDogLegStep(ParametersType & step, double & beta);

  ParametersType                   m_GaussNewtonStep;
  ParametersType                   m_SteepestDescentStep;
  ParametersType                   m_DogLegStep;

  double                           m_InitialTrustRegionRadius;
  double                           m_TrustRegionRadius;

private:
  ESMDogLegOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk

#include "itkESMDogLegOptimizer.txx"

#endif

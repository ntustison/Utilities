#ifndef MZ_ESMOptimizerBase_H_
#define MZ_ESMOptimizerBase_H_

#include "itkSingleValuedNonLinearOptimizer.h"
#include "itkSingleValuedVnlCostFunctionAdaptor.h"
#include "itkCommand.h"


namespace itk
{

/** \class ESMOptimizerBase
 * \brief This class is a base for the ESM Optimization methods
 *
 * \ingroup Numerics Optimizers
 */
template <class TFixedImage,  class TMovingImage>
class ITK_EXPORT ESMOptimizerBase :
    public SingleValuedNonLinearOptimizer
{
public:
  /** Standard class typedefs. */
  typedef ESMOptimizerBase                      Self;
  typedef SingleValuedNonLinearOptimizer        Superclass;
  typedef SmartPointer<Self>                    Pointer;
  typedef SmartPointer<const Self>              ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( ESMOptimizerBase, SingleValueNonLinearOptimizer );

  /**  Parameters type.
   *  It defines a position in the optimization search space. */
  typedef Superclass::ParametersType ParametersType;

  /** Images types */
  typedef TMovingImage                               MovingImageType;
  typedef TFixedImage                                FixedImageType;

  /** Type of the Cost Function   */
  typedef ESMImageToImageMetric<FixedImageType,
      MovingImageType>                               ESMCostFunctionType;
  typedef typename ESMCostFunctionType::Pointer      ESMCostFunctionPointer;

  /**  Measure type.
   *  It defines a type used to return the cost function value.  */
  typedef typename ESMCostFunctionType::MeasureType    MeasureType;

  /**  Derivative type.
   *  It defines a type used to return the cost function derivative. */
  typedef typename ESMCostFunctionType::DerivativeType DerivativeType;

  /**  Hessian type.
   *  It defines a type used to return the cost function derivative. */
  typedef typename ESMCostFunctionType::HessianType    HessianType;

  /** Codes of stopping conditions. */
  typedef enum {
     Unknown,
     ImageNotAvailable,
     CostFunctionError,
     MaximumNumberOfIterationsReached,
     RankDeficientHessianMatrix,
     RankDeficientHessianMatrixAndBadGradient,
     SoftStepSizeToleranceReached,
     MeasureToleranceReached,
     GradientToleranceReached,
     BadGainRatio
  } StopConditionType;

  /** Set the cost function. */
  itkSetObjectMacro( ESMCostFunction, ESMCostFunctionType );

  /** Get the cost function. */
  itkGetConstObjectMacro( ESMCostFunction, ESMCostFunctionType );

  /** Get the cost function value at the given parameters. */
  MeasureType GetValue( const ParametersType & parameters) const;


  itkSetMacro( MaximumNumberOfIterations, unsigned long );
  itkGetConstReferenceMacro( MaximumNumberOfIterations, unsigned long );

  itkGetConstMacro( CurrentIteration, unsigned int );

  itkGetConstReferenceMacro( StopCondition, StopConditionType );

  itkGetConstReferenceMacro( Value, MeasureType );
  itkGetConstReferenceMacro( Gradient, DerivativeType );
  itkGetConstReferenceMacro( Hessian, HessianType );

  itkGetConstReferenceMacro( OptimalValue, MeasureType );
  itkGetConstReferenceMacro( OptimalParameters, ParametersType );

  itkSetMacro( SoftStepSizeTolerance, double );
  itkGetConstReferenceMacro( SoftStepSizeTolerance, double );

  itkSetMacro( MeasureTolerance, double );
  itkGetConstReferenceMacro( MeasureTolerance, double );

  itkSetMacro( GradientTolerance, double );
  itkGetConstReferenceMacro( GradientTolerance, double );

  /** Get the reason for termination */
  virtual const std::string GetStopConditionDescription() const;

  /** Start optimization. */
  void StartOptimization();

  /** Stop optimization. */
  void StopOptimization();

protected:
  ESMOptimizerBase();
  virtual ~ESMOptimizerBase(){};

  /** Print out internal state */
  void PrintSelf(std::ostream& os, Indent indent) const;

  virtual void InitOptimization();
  virtual void GetCurrentValueDerivativeAndHessian();
  virtual void GetUpdate() = 0;
  virtual void ApplyUpdate();

  virtual bool GetGaussNewtonStep(ParametersType & step);
  virtual bool GetSteepestDescentStep(ParametersType & step, double & alpha);

  ESMCostFunctionPointer           m_ESMCostFunction;

  unsigned long                    m_MaximumNumberOfIterations;

  bool                             m_Stop;
  unsigned long                    m_CurrentIteration;
  ParametersType                   m_UpdateVector;

  StopConditionType                m_StopCondition;
  std::ostringstream               m_StopConditionDescription;

  MeasureType                      m_Value;
  DerivativeType                   m_Gradient;
  HessianType                      m_Hessian;

  MeasureType                      m_OptimalValue;
  ParametersType                   m_OptimalParameters;

  double                           m_SoftStepSizeTolerance;
  double                           m_MeasureTolerance;
  double                           m_GradientTolerance;

private:
  ESMOptimizerBase(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** This function is set private because it should not be used
   * Please use SetESmCostFunction instead */
  virtual void SetCostFunction(typename Superclass::CostFunctionType *){};
};

} // end namespace itk

#include "itkESMOptimizerBase.txx"

#endif

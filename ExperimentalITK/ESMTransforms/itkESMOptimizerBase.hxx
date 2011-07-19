#ifndef MZ_ESMOptimizerBase_TXX_
#define MZ_ESMOptimizerBase_TXX_

#include "itkESMOptimizerBase.h"

#include <vnl/algo/vnl_cholesky.h>

namespace itk
{

template <class TFixedImage, class TMovingImage>
ESMOptimizerBase<TFixedImage,TMovingImage>
::ESMOptimizerBase()
   :m_ESMCostFunction(0)
   ,m_MaximumNumberOfIterations(100)
   ,m_Stop(false)
   ,m_CurrentIteration(0)
   ,m_Value( NumericTraits<MeasureType>::max() )
   ,m_OptimalValue( NumericTraits<MeasureType>::max() )
   ,m_SoftStepSizeTolerance( 1e-8 )
   ,m_MeasureTolerance( 1e-8 )
   ,m_GradientTolerance( 1e-8 )
{
   m_StopConditionDescription.str("");
}

/**
 * Get the cost function value at the given parameters
 */
template <class TFixedImage, class TMovingImage >
typename ESMOptimizerBase<TFixedImage,TMovingImage>::MeasureType
ESMOptimizerBase<TFixedImage,TMovingImage>
::GetValue( const ParametersType & parameters ) const
{
  itkDebugMacro("Computing CostFunction value at " <<  parameters);

  if(!m_ESMCostFunction)
    {
    ExceptionObject ex;
    ex.SetLocation(__FILE__);
    ex.SetDescription("The cost function must be set prior to calling GetValue");
    throw ex;
    }

  return this->GetESMCostFunction()->GetValue(parameters);
}


template <class TFixedImage, class TMovingImage >
const std::string
ESMOptimizerBase<TFixedImage,TMovingImage>
::GetStopConditionDescription() const
{
  return m_StopConditionDescription.str();
}

template <class TFixedImage, class TMovingImage >
void
ESMOptimizerBase<TFixedImage,TMovingImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  if (m_ESMCostFunction)
  {
     os << indent << "Cost Function: " << m_ESMCostFunction.GetPointer() << std::endl;
  }
}


template <class TFixedImage, class TMovingImage >
void
ESMOptimizerBase<TFixedImage,TMovingImage>
::GetCurrentValueDerivativeAndHessian()
{
   try
   {
      m_ESMCostFunction->GetValueDerivativeAndHessian(
         this->GetCurrentPosition(), m_Value, m_Gradient, m_Hessian );

      if ( this->m_Value < this->m_OptimalValue )
      {
         this->m_OptimalValue = this->m_Value;
         this->m_OptimalParameters = this->GetCurrentPosition();
      }
   }
   catch( ExceptionObject & excp )
   {
      m_StopCondition = CostFunctionError;
      m_StopConditionDescription << "Cost function error after "
                                 << m_CurrentIteration
                                 << " iterations. "
                                 << excp.GetDescription();
      this->StopOptimization();
   }
}


template <class TFixedImage, class TMovingImage >
void
ESMOptimizerBase<TFixedImage,TMovingImage>
::StopOptimization()
{
   itkDebugMacro("StopOptimization - "<<this->GetStopConditionDescription());

   m_Stop = true;
   this->InvokeEvent( EndEvent() );
}


template <class TFixedImage, class TMovingImage >
void
ESMOptimizerBase<TFixedImage,TMovingImage>
::InitOptimization()
{
   itkDebugMacro("StartOptimization");

   this->m_Stop = false;

   this->InvokeEvent( StartEvent() );

   this->m_ESMCostFunction->Initialize();

   const unsigned int spaceDimension = this->m_ESMCostFunction->GetNumberOfParameters();

   this->m_CurrentPosition.set_size( spaceDimension );

   this->m_Gradient.set_size( spaceDimension );
   this->m_Hessian.set_size( spaceDimension, spaceDimension );

   this->m_UpdateVector.set_size( spaceDimension );

   this->SetCurrentPosition( this->GetInitialPosition() );

   this->m_OptimalValue = NumericTraits<MeasureType>::max();
   this->m_OptimalParameters = this->GetInitialPosition();
}


template <class TFixedImage, class TMovingImage >
bool
ESMOptimizerBase<TFixedImage,TMovingImage>
::GetGaussNewtonStep(ParametersType & step)
{
   vnl_cholesky cholesky(this->m_Hessian, vnl_cholesky::quiet);

   if ( cholesky.rank_deficiency() == 0 )
   {
      cholesky.solve(-this->m_Gradient, &(step));
      return true;
   }
   else
   {
      return false;
   }
}


template <class TFixedImage, class TMovingImage >
bool
ESMOptimizerBase<TFixedImage,TMovingImage>
::GetSteepestDescentStep(ParametersType & step, double & alpha)
{
   // Compute ||J.g||^2 = g^T.H.g
   const double jg2 =  bracket(this->m_Gradient, this->m_Hessian, this->m_Gradient);

   // Compute ||g||^2
   const double g2 = this->m_Gradient.squared_magnitude();

   // Compute ||g||^2 / ||J.g||^2 = ||g||^2 / g^T.H.g
   alpha = g2 / jg2;

   // Get the corresponding step
   step = (-alpha) * this->m_Gradient;

   // Check if we didn't do something badly behaved numerically
   if ( g2<=0.0 || jg2 <= std::numeric_limits<double>::epsilon()*g2 )
   {
      return false;
   }

   return true;
}


template <class TFixedImage, class TMovingImage >
void
ESMOptimizerBase<TFixedImage,TMovingImage>
::ApplyUpdate()
{
   // This step might be skipped if we don't tamper with m_ESMCostFunction
   this->m_ESMCostFunction->SetTransformParameters(this->m_CurrentPosition);

   this->m_ESMCostFunction->IncrementalESMTransformUpdate( this->m_UpdateVector );
   this->SetCurrentPosition( this->m_ESMCostFunction->GetESMTransform()->GetParameters() );

   //this->m_CurrentPosition += this->m_UpdateVector;

   //std::cout<<" position: "<<this->GetCurrentPosition()<<std::endl;
}


template <class TFixedImage, class TMovingImage >
void
ESMOptimizerBase<TFixedImage,TMovingImage>
::StartOptimization()
{
   this->InitOptimization();

   while( !m_Stop )
   {
      if( m_CurrentIteration >= m_MaximumNumberOfIterations )
      {
         m_StopCondition = MaximumNumberOfIterationsReached;
         m_StopConditionDescription << "Maximum number of iterations ("
                                    << m_MaximumNumberOfIterations
                                    << ") exceeded.";
         this->StopOptimization();
         break;
      }

      this->GetCurrentValueDerivativeAndHessian();

      if( m_Stop ) break;

      this->GetUpdate();

      if( m_Stop ) break;

      this->ApplyUpdate();

      this->InvokeEvent( IterationEvent() );

      if( m_Stop ) break;

      m_CurrentIteration++;
   }
}

} // end namespace itk

#endif

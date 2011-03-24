#ifndef MZ_ESMDogLegOptimizer_TXX_
#define MZ_ESMDogLegOptimizer_TXX_

#include "itkESMDogLegOptimizer.h"

namespace itk
{

template <class TFixedImage, class TMovingImage>
ESMDogLegOptimizer<TFixedImage,TMovingImage>
::ESMDogLegOptimizer()
   :Superclass()
   ,m_InitialTrustRegionRadius(10.0)
{
}

template <class TFixedImage, class TMovingImage >
void
ESMDogLegOptimizer<TFixedImage,TMovingImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

template <class TFixedImage, class TMovingImage >
void
ESMDogLegOptimizer<TFixedImage,TMovingImage>
::InitOptimization()
{
   this->Superclass::InitOptimization();

   const unsigned int spaceDimension = this->m_ESMCostFunction->GetNumberOfParameters();

   this->m_GaussNewtonStep.set_size( spaceDimension );
   this->m_SteepestDescentStep.set_size( spaceDimension );
   this->m_DogLegStep.set_size( spaceDimension );

   // Initial trust region radius is provided by the user
   this->m_TrustRegionRadius = this->m_InitialTrustRegionRadius;
}

template <class TFixedImage, class TMovingImage >
typename ESMDogLegOptimizer<TFixedImage,TMovingImage>::DogLegStepType
ESMDogLegOptimizer<TFixedImage,TMovingImage>
::GetDogLegStep(ParametersType & step, double & beta)
{
   const double nhgn = this->m_GaussNewtonStep.two_norm();

   if ( nhgn <= this->m_TrustRegionRadius )
   {
      step = this->m_GaussNewtonStep;
      return GaussNewtonStep;
   }

   const double nhsd2 = this->m_SteepestDescentStep.squared_magnitude();
   const double nhsd = std::sqrt( nhsd2 );

   if ( nhsd >= this->m_TrustRegionRadius )
   {
      step = (this->m_TrustRegionRadius/nhsd) * this->m_SteepestDescentStep;
      return ScaledToTrustRegionSteepestDescentStep;
   }

   const vnl_vector<double> sdgndiff( this->m_GaussNewtonStep - this->m_SteepestDescentStep );

   const double c = dot_product( this->m_SteepestDescentStep, sdgndiff );

   const double nsdgndiff2 = sdgndiff.squared_magnitude();

   const double delta2 = vnl_math_sqr(this->m_TrustRegionRadius);

   const double sqrtterm = std::sqrt( c*c + nsdgndiff2*(delta2-nhsd2) );

   if ( c <= 0.0 )
   {
      beta = (-c+sqrtterm) / nsdgndiff2;
   }
   else
   {
      beta = (delta2-nhsd2) / (c+sqrtterm);
   }

   step = this->m_SteepestDescentStep + beta*sdgndiff;

   return TrueDogLegStep;
}

template <class TFixedImage, class TMovingImage >
void
ESMDogLegOptimizer<TFixedImage,TMovingImage>
::GetUpdate()
{
   const bool gnok = this->GetGaussNewtonStep(this->m_GaussNewtonStep);

   double alpha(0.0);
   const bool sdok = this->GetSteepestDescentStep(this->m_SteepestDescentStep,alpha);

   if (!gnok)
   {
      //std::cout<<"bad gn step"<<std::endl;
      if (!sdok)
      {
         //std::cout<<"bad sd step"<<std::endl;
         this->m_StopCondition = Superclass::RankDeficientHessianMatrixAndBadGradient;
         this->m_StopConditionDescription << "The Hessian matrix is rank deficient. H = "
                                          << this->m_Hessian << ". The gradient is also badly scaled g = "
                                          << this->m_Gradient;
         this->StopOptimization();
         return;
      }

      // Hack. This case shouldn't happen often
      this->m_GaussNewtonStep = this->m_SteepestDescentStep;
   }

   if ( !sdok )
   {
      //std::cout<<"bad sd step"<<std::endl;
      // Hack. This case shouldn't happen often
      this->m_SteepestDescentStep = this->m_GaussNewtonStep;
   }

   double beta(0.0);
   const DogLegStepType steptype = this->GetDogLegStep(this->m_DogLegStep, beta);
   //std::cout<<"steptype = "<<steptype<<std::endl;

   // Check if the step size is small compared to eps2^2 or small
   // compared to eps2*x
   const double nh = this->m_DogLegStep.two_norm();
   const double nx = this->GetCurrentPosition().two_norm();
   if ( nh <= this->m_SoftStepSizeTolerance*(nx+this->m_SoftStepSizeTolerance) )
   {
      this->m_StopCondition = Superclass::SoftStepSizeToleranceReached;
      this->m_StopConditionDescription << "The step size is too small. ||h||=" << nh
                                       << " ||x||="<<nx;
      this->StopOptimization();
      return;
   }

   // Get a new proposal for the position
   this->m_ESMCostFunction->IncrementalESMTransformUpdate( this->m_DogLegStep );
   const ParametersType newparams = this->m_ESMCostFunction->GetESMTransform()->GetParameters();

   const MeasureType newvalue = this->GetValue( newparams );

   //std::cout<<"newparams="<<newparams
   //         <<"newvalue="<<newvalue<<std::endl;

   // Check if we accept it by looking at the gain ration
   double Ldiff;
   switch (steptype)
   {
   case GaussNewtonStep:
   {
      Ldiff = this->GetValue();
      break;
   }
   case ScaledToTrustRegionSteepestDescentStep:
   {
      // nhsd could be retrieved from GetDogLegStep
      const double nhsd = this->m_SteepestDescentStep.two_norm();

      Ldiff = this->m_TrustRegionRadius * (2*nhsd-this->m_TrustRegionRadius) / (2*alpha);
      break;
   }
   case TrueDogLegStep:
   {
      // nhsd2 could be retrieved from GetDogLegStep
      const double nhsd2 = this->m_SteepestDescentStep.squared_magnitude();

      Ldiff = vnl_math_sqr(1-beta)*nhsd2/(2*alpha)
         + beta*(2-beta)*this->GetValue();
      break;
   }
   default:
   {
      this->m_StopCondition = Superclass::Unknown;
      this->m_StopConditionDescription << "Unknown dog leg step type "<<steptype;
      this->StopOptimization();
      return;
   }
   }

   if ( Ldiff <= 0.0 )
   {
      this->m_StopCondition = Superclass::BadGainRatio;
      this->m_StopConditionDescription << "Wrong linear model, Ldiff=" << Ldiff;
      this->StopOptimization();
      return;
   }

   // The gain ratio
   const double pho = (this->GetValue() - newvalue) / Ldiff;

   if ( pho > 0.0 )
   {
      // We do get a better step
      this->m_OptimalValue = newvalue;
      this->m_OptimalParameters = newparams;

      this->m_UpdateVector = this->m_DogLegStep;

      // This replaces ApplyUpdate
      this->SetCurrentPosition( newparams );

      // Check if we can stop here
      if ( newvalue <= this->m_MeasureTolerance )
      {
         ///\todo use the infinity norm
         this->m_StopCondition = Superclass::MeasureToleranceReached;
         this->m_StopConditionDescription << "Measure tolerance reached. Value=" << newvalue;
         this->StopOptimization();
         return;
      }

      const double ng_inf = this->m_Gradient.inf_norm();
      if ( ng_inf <= this->m_GradientTolerance )
      {
         ///\todo use the infinity norm
         this->m_StopCondition = Superclass::GradientToleranceReached;
         this->m_StopConditionDescription << "Gradient tolerance reached. Value=" << ng_inf;
         this->StopOptimization();
         return;
      }
   }

   if ( pho > 0.75 )
   {
      // Increase the trust region if needed
      this->m_TrustRegionRadius = std::max(this->m_TrustRegionRadius, 3*nh);
   }
   else if ( pho < 0.25 )
   {
      // Decrease the trust region
      this->m_TrustRegionRadius /= 2.0;

      // Check if the new trust region is small compared to eps2^2 or small
      // compared to eps2*x
      if ( this->m_TrustRegionRadius <= this->m_SoftStepSizeTolerance*(nx+this->m_SoftStepSizeTolerance) )
      {
         this->m_StopCondition = Superclass::SoftStepSizeToleranceReached;
         this->m_StopConditionDescription << "The trust region is too small. Delta="
                                          << this->m_TrustRegionRadius
                                          << " ||x||="<<nx;
         this->StopOptimization();
         return;
      }
   }
}

} // end namespace itk

#endif

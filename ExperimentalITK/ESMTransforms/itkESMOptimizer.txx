#ifndef MZ_ESMOptimizer_TXX_
#define MZ_ESMOptimizer_TXX_

#include "itkESMOptimizer.h"

namespace itk
{

template <class TFixedImage, class TMovingImage>
ESMOptimizer<TFixedImage,TMovingImage>
::ESMOptimizer()
   :Superclass()
{
}

template <class TFixedImage, class TMovingImage >
void
ESMOptimizer<TFixedImage,TMovingImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

template <class TFixedImage, class TMovingImage >
void
ESMOptimizer<TFixedImage,TMovingImage>
::GetUpdate()
{
   if ( !this->GetGaussNewtonStep(this->m_UpdateVector) )
   {
      this->m_StopCondition = Superclass::RankDeficientHessianMatrix;
      this->m_StopConditionDescription << "The Hessian matrix is rank deficient. H = "
                                       << this->m_Hessian;
      this->StopOptimization();
      return;
   }

   // Check if the step size is small compared to eps2^2 or small
   // compared to eps2*x
   const double nh = this->m_UpdateVector.two_norm();
   const double nx = this->GetCurrentPosition().two_norm();
   if ( nh <= this->m_SoftStepSizeTolerance*(nx+this->m_SoftStepSizeTolerance) )
   {
      this->m_StopCondition = Superclass::SoftStepSizeToleranceReached;
      this->m_StopConditionDescription << "The step size is too small. ||h||=" << nh
                                       << " ||x||="<<nx;
      this->StopOptimization();
      return;
   }


   //normalize
   //this->m_UpdateVector = -this->m_Gradient;
   //this->m_UpdateVector[0] *= 1e-4;
   //this->m_UpdateVector.normalize();

   //std::cout<<"  update: "<<this->m_UpdateVector<<std::endl;
}

} // end namespace itk

#endif

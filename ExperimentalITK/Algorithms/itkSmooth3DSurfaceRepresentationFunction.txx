#ifndef __itkSmooth3DSurfaceRepresentationFunction_txx
#define __itkSmooth3DSurfaceRepresentationFunction_txx

#include "itkSmooth3DSurfaceRepresentationFunction.h"
#include "itkBinaryWellComposed3DImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkAbsImageFilter.h"
#include "itkCentralDifferenceImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkCastImageFilter.h"

#include "vnl/vnl_math.h"

namespace itk
{ 

template<class TInputImage, class TOutput, class TCoordRep>
void
Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>
::SetInputImage(const ImageType *ptr)
{
  Superclass::SetInputImage(ptr);

  // Ensure the input image is well-composed
  itkDebugMacro(<< "Making the image well composed.");
  typedef BinaryWellComposed3DImageFilter<ImageType> WellComposedFilterType;
  typename WellComposedFilterType::Pointer wc = WellComposedFilterType::New();
  wc->SetInput(this->m_Image);
  wc->Update();

  // Calculate the distance transform image.
  itkDebugMacro(<< "Calculating the distance transform.");
  typedef SignedMaurerDistanceMapImageFilter<ImageType, ImageType> DistanceFilterType;
  typename DistanceFilterType::Pointer edt = DistanceFilterType::New();
  edt->SetUseImageSpacing(true);
  edt->SetSquaredDistance(false);
  edt->SetInsideIsPositive(true);
  edt->SetInput(wc->GetOutput());
  edt->Update();

  // Re-assign the input image 
  this->m_Image = edt->GetOutput();  

  // Find the minimum value of the absolute distance image.  
  itkDebugMacro(<< "Finding the minimum of the absolute value of the distance image.");
  typedef AbsImageFilter<ImageType, RealImageType> AbsFilterType;
  typename AbsFilterType::Pointer abs = AbsFilterType::New();
  abs->SetInput(edt->GetOutput()); 
  abs->Update();

  typedef StatisticsImageFilter<RealImageType> StatisticsFilterType;
  typename StatisticsFilterType::Pointer stats = StatisticsFilterType::New();
  stats->SetInput(abs->GetOutput());  
  stats->Update();
  RealType fmin = stats->GetMinimum();

  typedef CastImageFilter<ImageType, RealImageType> CastFilterType;
  typename CastFilterType::Pointer caster = CastFilterType::New();
  caster->SetInput(this->m_Image);  

  // Find the maximum gradient component.  
  itkDebugMacro(<< "Finding the maximum gradient component of the distance image.");
  RealType gmax = 0.0;
  typedef CentralDifferenceImageFunction<RealImageType, double> InterpFunction;
  typename InterpFunction::Pointer deriv = InterpFunction::New();
  deriv->SetInputImage(caster->GetOutput());

  ImageRegionIteratorWithIndex<RealImageType> It(caster->GetOutput(), 
        caster->GetOutput()->GetRequestedRegion());
  for (It.GoToBegin(); !It.IsAtEnd(); ++It)
  {
    typename RealImageType::IndexType idx = It.GetIndex();
    typename InterpFunction::ContinuousIndexType cidx;

    for (int i = -1; i <= 1; i+=2)
    {
      for (int j = -1; j <= 1; j+=2)
      {   
        for (int k = -1; k <= 1; k+=2)
        {  
          cidx[0] = idx[0] + static_cast<RealType>(i)*0.4;   
          cidx[1] = idx[1] + static_cast<RealType>(j)*0.4;   
          cidx[2] = idx[2] + static_cast<RealType>(k)*0.4;   

          for (unsigned  m = 0; m < ImageDimension; m++)
          {
            gmax = vnl_math_max(gmax, deriv->EvaluateAtContinuousIndex(cidx)[i]);  
          }  
        }
      }
    }
  }              
  this->m_d = vnl_math_min(0.2, fmin/(2*sqrt(3)*gmax));
}

template<class TInputImage, class TOutput, class TCoordRep>
TOutput
Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>
::Evaluate(const PointType &point) const
{
  itkDebugMacro("Evaluating function at " << point);

  IndexType idx, tmp;
  bool isInside = this->m_Image->TransformPhysicalPointToIndex(point, idx);
  if (!isInside)
  {
    itkExceptionMacro( << "The specified point is outside the function domain."); 
  }

  typedef CastImageFilter<ImageType, RealImageType> CastFilterType;
  typename CastFilterType::Pointer caster = CastFilterType::New();
  caster->SetInput(this->m_Image);  

  RealType num = 0.0;
  RealType den = 0.0;
  for (int i = -1; i <= 1; i++)
  {
    tmp[0] = idx[0] + i;
    for (int j = -1; j <= 1; j++)
    {
      tmp[1] = idx[1] + j;
      for (int k = -1; k <= 1; k++)
      {
        tmp[2] = idx[2] + k;
         
        if (this->m_Image->GetRequestedRegion().IsInside(tmp))
        {
          PointType localPoint = this->FromGlobalToLocal(point, tmp);
          RealType ww = this->WeightFunction(localPoint);
          if (ww > 0.0)
          {
            num += (this->FunctionLinearInterpolation(localPoint, tmp)*ww);
            den += ww;        
          }
        }
      }
    }
  }     
  return (den != 0.0) ? num/den : 0.0;
} 

template<class TInputImage, class TOutput, class TCoordRep>
typename Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>::PointType
Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>
::EvaluateGradient(const PointType &point) const
{
  itkDebugMacro("Evaluating function gradient at " << point);

  IndexType idx, tmp;
  bool isInside = this->GetInput()->TransformPhysicalPointToIndex(point, idx);
  if (!isInside)
  {
    itkExceptionMacro( << "The specified point is outside the function domain."); 
  }
    
  typedef LinearInterpolateImageFunction<RealImageType, double> InterpFunction;
  typename InterpFunction::Pointer interp = InterpFunction::New();
  interp->SetInputImage(this->m_Image);

  PointType grad(0.0);  

  for (int i = -1; i <= 1; i++)
  {
    tmp[0] = idx[0] + i;
    for (int j = -1; j <= 1; j++)
    {
      tmp[1] = idx[1] + j;
      for (int k = -1; k <= 1; k++)
      {
        tmp[2] = idx[2] + k;

        if (this->m_Image->GetRequestedRegion().IsInside(tmp))
        {
          PointType localPoint = this->FromGlobalToLocal(point, tmp);
          RealType ww = this->WeightFunction(localPoint);
          if (ww > 0.0)
          { 
            RealType f = interp->EvaluateAtContinuousIndex(tmp, localPoint);
            PointType dw = this->WeightFunctionDerivative(localPoint);
            PointType df = this->GradientLinearInterpolation(point, tmp);
            grad += (dw*f + df*ww);
          }
        }
      }
    }
  }     
  grad /= this->m_Image->GetSpacing()[0];
  return grad;
}

template<class TInputImage, class TOutput, class TCoordRep>
typename Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>::PointType
Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>
::FromGlobalToLocal(const PointType& pt, const IndexType& idx) const
{
  PointType tmp;
 
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    tmp[i] = (pt[i] - (this->m_Image->GetOrigin()[i] + static_cast<RealType>(idx[i])*this->m_Image->GetSpacing()[i]));
    tmp[i] /= this->m_Image->GetSpacing()[i];
  }
  return tmp;
}

template<class TInputImage, class TOutput, class TCoordRep>
typename Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>::RealType
Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>
::FunctionLinearInterpolation(const PointType& pt, const IndexType& idx) const
{
  typename ImageType::OffsetType offsets[8]; 

  offsets[0][0] = 0;
  offsets[0][1] = 0;
  offsets[0][2] = 0;

  offsets[1][0] = 0;
  offsets[1][1] = 0;
  offsets[1][2] = 1;

  offsets[2][0] = 0;
  offsets[2][1] = 1;
  offsets[2][2] = 0;

  offsets[3][0] = 0;
  offsets[3][1] = 1;
  offsets[3][2] = 1;

  offsets[4][0] = 1;
  offsets[4][1] = 0;
  offsets[4][2] = 0;

  offsets[5][0] = 1;
  offsets[5][1] = 0;
  offsets[5][2] = 1;

  offsets[6][0] = 1;
  offsets[6][1] = 1;
  offsets[6][2] = 0;

  offsets[7][0] = 1;
  offsets[7][1] = 1;
  offsets[7][2] = 1;

  RealType f[8];
  
  for (unsigned int i = 0; i < 8; i++)
  {
    if (this->m_Image->GetRequestedRegion().IsInside(idx+offsets[i]))
    { 
      f[i] = this->m_Image->GetPixel(idx+offsets[i]);
    }
    else
    {
      return 0.0;
    }
  }
  
  return (f[0] + (f[1]-f[0])*pt[2] + (f[2]-f[0])*pt[1] + (f[4]-f[0])*pt[0]
         + (f[0]+f[3]-f[1]-f[2])*pt[1]*pt[2] 
         + (f[0]+f[5]-f[1]-f[4])*pt[0]*pt[2]
         + (f[0]+f[6]-f[2]-f[4])*pt[0]*pt[1]
         + (f[7]+f[4]+f[2]+f[1]-f[6]-f[5]-f[3]-f[0])*pt[0]*pt[1]*pt[2]);
}

template<class TInputImage, class TOutput, class TCoordRep>
typename Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>::PointType
Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>
::GradientLinearInterpolation(const PointType& pt, const IndexType& idx) const
{
  IndexType offsets[8]; 

  offsets[0][0] = 0;
  offsets[0][1] = 0;
  offsets[0][2] = 0;

  offsets[1][0] = 0;
  offsets[1][1] = 0;
  offsets[1][2] = 1;

  offsets[2][0] = 0;
  offsets[2][1] = 1;
  offsets[2][2] = 0;

  offsets[3][0] = 0;
  offsets[3][1] = 1;
  offsets[3][2] = 1;

  offsets[4][0] = 1;
  offsets[4][1] = 0;
  offsets[4][2] = 0;

  offsets[5][0] = 1;
  offsets[5][1] = 0;
  offsets[5][2] = 1;

  offsets[6][0] = 1;
  offsets[6][1] = 1;
  offsets[6][2] = 0;

  offsets[7][0] = 1;
  offsets[7][1] = 1;
  offsets[7][2] = 1;

  RealType f[8];

  PointType g;
  g.Fill(0);
  
  for (unsigned int i = 0; i < 8; i++)
  {
    if (this->m_Image->GetRequestedRegion().IsInside(idx+offsets[i]))
    { 
      f[i] = this->m_Image->GetPixel(idx+offsets[i]);
    }
    else
    {
      return g;
    }
  }  
  
  g[0] = (f[4]-f[0]) + (f[0]+f[5]-f[1]-f[4])*pt[2] 
         + (f[0]+f[6]-f[2]-f[4])*pt[1] 
         + (f[7]+f[4]+f[2]+f[1]-f[6]-f[5]-f[3]-f[0])*pt[1]*pt[2];
  g[1] = (f[2]-f[0]) + (f[0]+f[3]-f[1]-f[2])*pt[2] 
         + (f[0]+f[6]-f[2]-f[4])*pt[1] 
         + (f[7]+f[4]+f[2]+f[1]-f[6]-f[5]-f[3]-f[0])*pt[0]*pt[2];
  g[2] = (f[1]-f[0]) + (f[0]+f[3]-f[1]-f[2])*pt[1] 
         + (f[0]+f[5]-f[1]-f[4])*pt[0] 
         + (f[7]+f[4]+f[2]+f[1]-f[6]-f[5]-f[3]-f[0])*pt[0]*pt[1];

  return g;
}

template<class TInputImage, class TOutput, class TCoordRep>
typename Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>::RealType
Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>
::HFunction(const RealType &x) const
{
  RealType val = 0.0;
  if (x <= 0.0)
  {
    val = 1.0;
  }
  else if (x < 1)
  {
    val = exp(3.0*exp(-1.0/x)/(x-1.0));
  }
  return val;
}

template<class TInputImage, class TOutput, class TCoordRep>
typename Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>::RealType
Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>
::GFunction(const RealType &x) const
{
  RealType val = 0.0;
  if (x <= (0.5-this->m_d))
  {
    val = 1.0;
  }
  else if (x < (0.5+this->m_d))
  {
    RealType d1 = 0.5 - this->m_d;
    RealType d2 = (x - d1)/(2 * this->m_d);
    RealType h1 = this->HFunction(d2);
    RealType h2 = this->HFunction(1.0 - d2);
    val = exp(3.0*exp(-1.0/x)/(x - 1.0));
  }
  return val;
}

template<class TInputImage, class TOutput, class TCoordRep>
typename Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>::RealType
Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>
::WeightFunction(const PointType &pt) const
{
  RealType g = 1.0;

  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    if (pt[i] < -this->m_d || pt[i] > (1.0 + this->m_d))
    {
      return 0.0;
    }
    else
    {
      g *= this->GFunction(vnl_math_abs(pt[i]-0.5));
    }  
  }
  return g;
}

template<class TInputImage, class TOutput, class TCoordRep>
typename Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>::RealType
Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>
::HFunctionDerivative(const RealType &x) const
{
  RealType val = 0.0;
  if ((x > 0.0) && (x < 1.0))
  {
    RealType tmp = (3.0*exp(-1.0/x)/(1.0-x))-(1.0/x);
    val = -3.0*exp(tmp)*((x*tmp)+1.0)/(tmp*tmp*x*x);
  }
  return val;
}

template<class TInputImage, class TOutput, class TCoordRep>
typename Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>::RealType
Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>
::GFunctionDerivative(const RealType &x) const
{
  RealType val = 0.0;
  if ((x > 0.0) && (x < (0.5 + this->m_d)))
  {
    RealType a1 = 0.5 - this->m_d;
    RealType a2 = (x - a1)/(2 * this->m_d);
    RealType h1 =  this->HFunction(a2);
    RealType h2 =  this->HFunction(1.0 - a2);
    RealType d1 =  this->HFunctionDerivative(a2)/(2 * this->m_d);
    RealType d2 = -this->HFunctionDerivative(1.0 - a2)/(2 * this->m_d);
    
    val = d1*(h1 + h2) - h1*(d1 + d2);
    val /= ((h1+h2)*(h1+h2));
  }
  return val;
}

template<class TInputImage, class TOutput, class TCoordRep>
typename Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>::PointType
Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>
::WeightFunctionDerivative(const PointType &pt) const
{
  PointType d(0.0);
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    if (pt[i] < -this->m_d || pt[i] > (1.0+this->m_d))
    {
      return d;
    }
  }
 
  pt -= 0.5;
   
  PointType s(ImageDimension);
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    s[i] = (pt[i] < 0) ? -1.0 : 1.0;
  }

  d[0] = s[0]*this->GFunctionDerivative(s[0]*pt[0]) 
         *this->GFunction(s[1]*pt[1])*this->GFunction(s[2]*pt[2]);
  d[1] = s[1]*this->GFunctionDerivative(s[1]*pt[1]) 
         *this->GFunction(s[2]*pt[2])*this->GFunction(s[0]*pt[0]);
  d[2] = s[2]*this->GFunctionDerivative(s[2]*pt[2]) 
         *this->GFunction(s[0]*pt[0])*this->GFunction(s[1]*pt[1]);

  return d;
}

template <class TInputImage, class TOutput, class TCoordRep>
void
Smooth3DSurfaceRepresentationFunction<TInputImage, TOutput, TCoordRep>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);     
}

} // end namespace itk

#endif

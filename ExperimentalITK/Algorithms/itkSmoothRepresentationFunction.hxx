#ifndef __itkSmoothRepresentationFunction_hxx
#define __itkSmoothRepresentationFunction_hxx

#include "itkSmoothRepresentationFunction.h"
#include "itkStatisticsImageFilter.h"
#include "itkAbsImageFilter.h"
#include "itkCentralDifferenceImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkCastImageFilter.h"

#include "vnl/vnl_math.h"

namespace itk
{ 

template<class TInputImage, class TOutput, class TCoordRep>
SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>
::SmoothRepresentationFunction()
{
}

template<class TInputImage, class TOutput, class TCoordRep>
SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>
::~SmoothRepresentationFunction()
{
}

template<class TInputImage, class TOutput, class TCoordRep>
void
SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>
::SetInputImage( const ImageType *ptr )
{
  Superclass::SetInputImage( ptr );

  itkDebugMacro( << "Finding the minimum of the absolute value of the input image." );
  typedef AbsImageFilter<ImageType, RealImageType> AbsFilterType;
  typename AbsFilterType::Pointer abs = AbsFilterType::New();
  abs->SetInput( this->m_Image ); 
  abs->Update();

  typedef StatisticsImageFilter<RealImageType> StatisticsFilterType;
  typename StatisticsFilterType::Pointer stats = StatisticsFilterType::New();
  stats->SetInput( abs->GetOutput() );  
  stats->Update();
  RealType minF = stats->GetMinimum();

  typedef CastImageFilter<ImageType, RealImageType> CastFilterType;
  typename CastFilterType::Pointer caster = CastFilterType::New();
  caster->SetInput( this->m_Image );  
  caster->Update();

  // Find the maximum gradient component.  
  itkDebugMacro( << "Finding the maximum gradient component of the distance image." );
  typedef CentralDifferenceImageFunction<RealImageType, double> DerivativeFunction;
  typename DerivativeFunction::Pointer deriv = DerivativeFunction::New();
  deriv->SetInputImage( caster->GetOutput() );

  typedef ConstNeighborhoodIterator<RealImageType> NeighborhoodIteratorType;
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );

  RealType maxG = 0.0;
  
  NeighborhoodIteratorType It( radius, caster->GetOutput(), 
        caster->GetOutput()->GetRequestedRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    typename RealImageType::IndexType idx = It.GetIndex();
    typename DerivativeFunction::ContinuousIndexType cidx;
 
    for (unsigned int i = 0; i < It.GetNeighborhood().Size(); i++)
      {
      typename NeighborhoodIteratorType::OffsetType offset;
      if ( this->m_Image->GetRequestedRegion().IsInside( It.GetIndex( i ) ) )
        {
        offset = It.GetOffset( i );      
        bool isCorner = true;
        for ( unsigned int j = 0; j < ImageDimension; j++ )
          {
          if ( offset[j] == 0 )
            {
            isCorner = false;
            break;
            }
          }
        if ( isCorner )
          {
          for ( unsigned int j = 0; j < ImageDimension; j++ )
            {
            cidx[j] = static_cast<RealType>( It.GetIndex()[j] ) 
                    + static_cast<RealType>( offset[j] )*0.4;
            }
          if (this->m_Image->GetRequestedRegion().IsInside(cidx))
            {
            for (unsigned int j = 0; j < ImageDimension; j++)
              {      
              maxG = vnl_math_max(maxG, deriv->EvaluateAtContinuousIndex(cidx)[j]);  
              }  
            }
          }
        }
      }
    }              
  this->m_d = vnl_math_min( 0.2, minF/( 2.0*sqrt( 3.0 )*maxG ) );
}

template<class TInputImage, class TOutput, class TCoordRep>
TOutput
SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>
::EvaluateAtIndex( const IndexType &idx ) const
{
  PointType pt;
  this->m_Image->TransformIndexToPhysicalPoint( idx, pt );
  return this->Evaluate( pt );
}

template<class TInputImage, class TOutput, class TCoordRep>
TOutput
SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>
::EvaluateAtContinuousIndex( const ContinuousIndexType &cidx ) const
{
  PointType pt;
  this->m_Image->TransformContinuousIndexToPhysicalPoint( cidx, pt );
  return this->Evaluate( pt );
}

template<class TInputImage, class TOutput, class TCoordRep>
TOutput
SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>
::Evaluate( const PointType &point ) const
{
  itkDebugMacro( "Evaluating function at " << point );

  IndexType idx, tmp;
  bool isInside = this->m_Image->TransformPhysicalPointToIndex( point, idx );
  if ( !isInside )
    {
    itkExceptionMacro( << "The specified point is outside the function domain." ); 
    }

  typedef ConstNeighborhoodIterator<ImageType> NeighborhoodIteratorType;
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );

  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    idx[i] = static_cast<unsigned int>( 
      ( point[i] - this->m_Image->GetOrigin()[i] )
      / this->m_Image->GetSpacing()[i] );
    if ( idx[i] == this->m_Image->GetRequestedRegion().GetSize()[i] - 1 )
      {
      idx[i]--; 
      }
    }

  NeighborhoodIteratorType It( radius, this->m_Image, 
               this->m_Image->GetRequestedRegion() );
  It.SetLocation( idx );

  RealType num = 0.0;
  RealType den = 0.0;
  for ( unsigned int i = 0; i < pow( 3, ImageDimension ); i++ )
    {
    PointType localPoint = this->FromGlobalToLocal( point, It.GetIndex( i ) );
    RealType W = this->WeightFunction( localPoint );
    if ( W > 0.0 )
      {
      num += ( this->LinearInterpolationFunction( localPoint, It.GetIndex( i ) )*W );
      den += W;
      }
    }
  if ( den != NumericTraits<RealType>::Zero )
    {
    return num/den;
    }
  return 0.0;         
} 

template<class TInputImage, class TOutput, class TCoordRep>
typename SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>::PointType
SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>
::EvaluateGradientAtIndex( const IndexType &idx ) const
{
  PointType pt;
  this->m_Image->TransformIndexToPhysicalPoint( idx, pt );
  return this->EvaluateGradient( pt );
}

template<class TInputImage, class TOutput, class TCoordRep>
typename SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>::PointType
SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>
::EvaluateGradientAtContinuousIndex( const ContinuousIndexType &cidx ) const
{
  PointType pt;
  this->m_Image->TransformContinuousIndexToPhysicalPoint( cidx, pt );
  return this->EvaluateGradient( pt );
}

template<class TInputImage, class TOutput, class TCoordRep>
typename SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>::PointType
SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>
::EvaluateGradient( const PointType &point ) const
{
  itkDebugMacro( "Evaluating function gradient at " << point );

  IndexType idx, tmp;
  bool isInside = this->m_Image->TransformPhysicalPointToIndex( point, idx );
  if ( !isInside )
    {
    itkExceptionMacro( << "The specified point is outside the function domain."); 
    }

  typedef ConstNeighborhoodIterator<ImageType> NeighborhoodIteratorType;
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );

  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    idx[i] = static_cast<unsigned int>( 
      ( point[i] - this->m_Image->GetOrigin()[i] )
      / this->m_Image->GetSpacing()[i] );
    if ( idx[i] == this->m_Image->GetRequestedRegion().GetSize()[i] - 1 )
      {
      idx[i]--; 
      }
    }

  NeighborhoodIteratorType It( radius, this->m_Image, 
               this->m_Image->GetRequestedRegion() );
  It.SetLocation( idx );
    
  PointType grad; 
  grad.Fill( 0.0 );
  for ( unsigned int i = 0; i < pow( 3, ImageDimension ); i++ )
    {
    PointType localPoint = this->FromGlobalToLocal( point, It.GetIndex( i ) );
    RealType W = this->WeightFunction( localPoint );
    if ( W > 0.0 )
      { 
      RealType F = this->LinearInterpolationFunction( localPoint, It.GetIndex( i ) );
      PointType dW = this->WeightFunctionDerivative( localPoint );
      PointType dF = this->LinearInterpolationGradient( localPoint, It.GetIndex( i ) );
      for ( unsigned int j = 0; j < ImageDimension; j++ )
        {
        grad[j] += ( dW[j]*F + dF[j]*W );
        }
      }
    }     
  for (unsigned int i = 0; i < ImageDimension; i++)
    {
    grad[i] /= this->m_Image->GetSpacing()[i];  // Assuming isotropic voxels
    }    

  return grad;
} 

template<class TInputImage, class TOutput, class TCoordRep>
typename SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>::RealType
SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>
::LinearInterpolationFunction(const PointType& pt, const IndexType& idx) const
{
  IndexType K;
  K[0] = 1;
  for (unsigned int i = 1; i < ImageDimension; i++)
  {
    K[i] = K[i-1]*2;
  }
  
  RealType sum = 0.0;

  for ( unsigned int i = 0; i < pow( 2.0, static_cast<RealType>( ImageDimension ) ); i++ )
    {
    unsigned int index = i;  
    typename ImageType::OffsetType offset;
    for ( unsigned int j = 0; j < ImageDimension; j++ )
      {
      offset[ImageDimension-j-1] = static_cast<unsigned int>( index / K[ImageDimension-j-1] );
      index %= K[ImageDimension-j-1];
      }

    if ( this->m_Image->GetRequestedRegion().IsInside( idx+offset ) )
      { 
      RealType val = 1.0;
      for ( unsigned int j = 0; j < ImageDimension; j++ ) 
        {
        if ( offset[j] = 0 )
          {
          val *= ( 1.0 - pt[j] );
          }
        else
          {
          val *= pt[j];
          }
        }
      sum += ( val*this->m_Image->GetPixel( idx+offset ) );
    }
  }  
  return sum;
}

template<class TInputImage, class TOutput, class TCoordRep>
typename SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>::PointType
SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>
::LinearInterpolationGradient(const PointType& pt, const IndexType& idx) const
{
  IndexType K;
  K[0] = 1;
  for ( unsigned int i = 1; i < ImageDimension; i++ )
    {
    K[i] = K[i-1]*2;
    }
  
  PointType grad;
  grad.Fill(0.0);

  for ( unsigned int d = 0; d < ImageDimension; d++ )
    {
    RealType sum = 0.0;
    for ( unsigned int i = 0; i < pow( 2.0, static_cast<RealType>( ImageDimension ) ); i++ )
      {
      unsigned int index = i;  
      typename ImageType::OffsetType offset;
      for ( unsigned int j = 0; j < ImageDimension; j++ )
        {
        offset[ImageDimension-j-1] = static_cast<unsigned int>(index/K[ImageDimension-j-1] );
        index %= K[ImageDimension-j-1];
        }
      if ( this->m_Image->GetRequestedRegion().IsInside( idx+offset ) )
        { 
        RealType val = 1.0;
        for ( unsigned int j = 0; j < ImageDimension; j++ ) 
          {
          if ( j != d )
            {
            if ( offset[j] == 0 )
              {
              val *= ( 1.0 - pt[j] );
              }
            else
              {
              val *= pt[j];
              }
            }
          else
            {
            if ( offset[j] == 0 )
              {
              val *= -1.0;
              }
            else
              {
              val *= 1.0;
              }
            }
          }
        sum += ( val*this->m_Image->GetPixel( idx+offset ) );
        }
      }
    grad[d] = sum;  
    }
  return grad;
}

template<class TInputImage, class TOutput, class TCoordRep>
typename SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>::PointType
SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>
::FromGlobalToLocal( const PointType& pt, const IndexType& idx ) const
{
  PointType tmp;
  typename ImageType::SpacingType spacing = this->m_Image->GetSpacing();
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    tmp[i] = ( pt[i] - ( this->m_Image->GetOrigin()[i] + 
              static_cast<RealType>( idx[i] )*spacing[i] ) )
              / spacing[i];
    }
  return tmp;
}

template<class TInputImage, class TOutput, class TCoordRep>
typename SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>::RealType
SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>
::HFunction( const RealType &x ) const
{
  RealType H = 0.0;
  if ( x <= 0.0 )
    {
    H = 1.0;
    }
  else if ( x < 1 )
    {
    H = exp( 3.0*exp( -1.0/x )/( x-1.0 ) );
    }
  return H;
}

template<class TInputImage, class TOutput, class TCoordRep>
typename SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>::RealType
SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>
::GFunction( const RealType &x ) const
{
  RealType G;
  if ( x <= 0.5-this->m_d )
    {
    G = 1.0;
    }
  else if ( x < 0.5+this->m_d )
    {
    RealType d1 = 0.5 - this->m_d;
    RealType d2 = ( x - d1 )/( 2 * this->m_d );
    RealType H1 = this->HFunction( d2 );
    RealType H2 = this->HFunction( 1.0 - d2 );
    G = H1 / ( H1 + H2 );
    }
  return G;
}

template<class TInputImage, class TOutput, class TCoordRep>
typename SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>::RealType
SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>
::WeightFunction( const PointType &pt ) const
{
  RealType W = 1.0;
  for (unsigned int i = 0; i < ImageDimension; i++)
    {
    if ( pt[i] < -this->m_d || pt[i] > 1.0 + this->m_d )
      {
      return 0.0;
      }
    else
      {
      W *= this->GFunction( vnl_math_abs( pt[i]-0.5 ) );
      }  
    }
  return W;
}

template<class TInputImage, class TOutput, class TCoordRep>
typename SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>::RealType
SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>
::HFunctionDerivative(const RealType &x) const
{
  RealType dH = 0.0;
  if ( x > 0.0 && x < 1.0 )
    {
    RealType tmp = ( 3.0*exp( -1.0/x )/( 1.0-x ) )-( 1.0/x );
    dH = -3.0*exp( tmp )*( x*( x-1.0 )+1.0 )/( ( x-1.0 )*( x-1.0 )*tmp*x*x );
    }
  return dH;
}

template<class TInputImage, class TOutput, class TCoordRep>
typename SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>::RealType
SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>
::GFunctionDerivative( const RealType &x ) const
{
  RealType dG = 0.0;
  if ( x > 0.5 - this->m_d && x < 0.5 + this->m_d )
    {
    RealType a1 = 0.5 - this->m_d;
    RealType a2 = ( x - a1 ) / ( 2 * this->m_d );
    RealType H1 = this->HFunction( a2 );
    RealType H2 = this->HFunction( 1.0 - a2 );
    RealType dH1 =  this->HFunctionDerivative( a2 )/( 2 * this->m_d );
    RealType dH2 = -this->HFunctionDerivative( 1.0 - a2 )/( 2 * this->m_d );
    
    dG = dH1*(H1 + H2) - H1*(dH1 + dH2);
    dG /= ( ( H1+H2 )*( H1+H2 ) );
    }
  return dG;
}

template<class TInputImage, class TOutput, class TCoordRep>
typename SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>::PointType
SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>
::WeightFunctionDerivative( const PointType &pt ) const
{
  PointType dW;
  dW.Fill( 0.0 );
 
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    if ( pt[i] < -this->m_d || pt[i] > ( 1.0+this->m_d ) )
      {
      return dW;
      }
    }
  
  PointType s;
  PointType G;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {     
    if ( pt[i] < 0.5 )
      {
      s[i] = -1.0;
      }
    else
      {
      s[i] = 1.0;
      }
    G[i] = this->GFunction( s[i]*( pt[i]-0.5 ) );
    dW[i] = s[i]*this->GFunctionDerivative( s[i]*( pt[i]-0.5 ) );
    } 
  
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    for ( unsigned int j = 0; j < ImageDimension; j++ )
      {
      if ( i != j )
        {
        dW[i] *= G[j];
        }
      }
    }

  return dW;
}

template <class TInputImage, class TOutput, class TCoordRep>
void
SmoothRepresentationFunction<TInputImage, TOutput, TCoordRep>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);     
}

} // end namespace itk

#endif

#ifndef _itkWebSplineWeightFunction_hxx_
#define _itkWebSplineWeightFunction_hxx_

// disable debug warnings in MS compiler
#ifdef _MSC_VER
#pragma warning(disable: 4786)
#endif
 
#include "itkWebSplineWeightFunction.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkImageRegionIterator.h"

namespace itk {

template<class TInputImage, class TOutput, class TCoordRep>
WebSplineWeightFunction<TInputImage, TOutput, TCoordRep>
::WebSplineWeightFunction()
{
  this->m_Gamma = 1.0;
  this->m_Delta = 1.0;

  this->m_Interpolator = InterpolateFunctionType::New();
}

template<class TInputImage, class TOutput, class TCoordRep>
WebSplineWeightFunction<TInputImage, TOutput, TCoordRep>
::~WebSplineWeightFunction()
{
}

template<class TInputImage, class TOutput, class TCoordRep>
void
WebSplineWeightFunction<TInputImage, TOutput, TCoordRep>
::SetInputImage(const ImageType *ptr)
{
  Superclass::SetInputImage(ptr);
  
  itkDebugMacro( << "Calculating the signed distance transform." );

  typedef SignedMaurerDistanceMapImageFilter
          <ImageType, RealImageType> DistanceFilterType;
  typename DistanceFilterType::Pointer 
           distancer = DistanceFilterType::New();
  distancer->SetInput( this->m_Image );
  distancer->SetSquaredDistance( false );
  distancer->SetBackgroundValue( this->m_BackgroundValue );
  distancer->SetInsideIsPositive( true );
  distancer->Update();

  ImageRegionIterator<RealImageType> It( distancer->GetOutput(),
    distancer->GetOutput()->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( It.Get() < 0.0 )
      {
      It.Set( 0.0 ); 
      }
    }
  this->m_Interpolator->SetInputImage( distancer->GetOutput() ); 
}

template<class TInputImage, class TOutput, class TCoordRep>
TOutput
WebSplineWeightFunction<TInputImage, TOutput, TCoordRep>
::Evaluate( const PointType &point ) const
{
  RealType E = 1.0 - this->m_Interpolator->Evaluate( point )/this->m_Delta;
  if ( E > 0.0 )
    {
    return static_cast<OutputType>( 1.0 - pow( E, this->m_Gamma ) );
    }
  else
    {
    return NumericTraits<OutputType>::One;
    }  
} 

template<class TInputImage, class TOutput, class TCoordRep>
typename WebSplineWeightFunction <TInputImage, TOutput, TCoordRep>::PointType
WebSplineWeightFunction<TInputImage, TOutput, TCoordRep>
::EvaluateGradient( const PointType &point ) const
{
  itkDebugMacro( "Evaluating function gradient at " << point );

  PointType dF;
  dF.Fill( 0.0 );

  RealType E = 1.0 - this->m_Interpolator->Evaluate( point )/this->m_Delta;
  if ( E > 0.0 )
    {
    dF = this->m_Interpolator->EvaluateGradient( point );
    RealType t = this->m_Gamma * pow( E, this->m_Gamma-1 ) / this->m_Delta;
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      dF[i] *= t;  
      }
    }
  return dF;
}  

template <class TInputImage, class TOutput, class TCoordRep>
void
WebSplineWeightFunction<TInputImage, TOutput, TCoordRep>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent  );
  os << indent << "Background value: "
     << this->m_BackgroundValue << std::endl;
  os << indent << "Delta: "
     << this->m_Delta << std::endl;
  os << indent << "Gamma: "
     << this->m_Gamma << std::endl;
}

}
#endif

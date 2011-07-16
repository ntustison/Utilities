#ifndef __itkMultipleLabelToDistanceMapImageFilter_hxx
#define __itkMultipleLabelToDistanceMapImageFilter_hxx

#include "itkMultipleLabelToDistanceMapImageFilter.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkStatisticsImageFilter.h"

#include "vnl/vnl_math.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
MultipleLabelToDistanceMapImageFilter<TInputImage, TOutputImage>
::MultipleLabelToDistanceMapImageFilter() : m_UseImageSpacing( true ),
                                            m_SquaredDistance( false ),
                                            m_NormalizeImage( true ),
                                            m_Sigma( 1.0 )
{
}

template <class TInputImage, class TOutputImage>
MultipleLabelToDistanceMapImageFilter<TInputImage, TOutputImage>
::~MultipleLabelToDistanceMapImageFilter()
{  
}

template <class TInputImage, class TOutputImage>
void
MultipleLabelToDistanceMapImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  typename OutputImageType::Pointer output = OutputImageType::New();
  output->SetOrigin( this->GetInput()->GetOrigin() );
  output->SetRegions( this->GetInput()->GetRequestedRegion() );
  output->SetSpacing( this->GetInput()->GetSpacing() );
  output->Allocate();

  typedef StatisticsImageFilter<InputImageType> StatsFilterType;
  typename StatsFilterType::Pointer stats = StatsFilterType::New();
  stats->SetInput( this->GetInput() );
  stats->Update();

  OutputPixelType distances( stats->GetMaximum() );
  distances.fill( 0 );
  output->FillBuffer( distances );
  
  for ( unsigned int label = 1; label <= stats->GetMaximum(); label++ )
    {

    typedef BinaryThresholdImageFilter
      <InputImageType, InputImageType> ThresholderType;
    typename ThresholderType::Pointer thresholder = ThresholderType::New();
    thresholder->SetInput( this->GetInput() );
    thresholder->SetInsideValue( 1 );
    thresholder->SetOutsideValue( 0 );
    thresholder->SetLowerThreshold( label );
    thresholder->SetUpperThreshold( label );
    thresholder->Update();

    typedef SignedMaurerDistanceMapImageFilter
      <InputImageType, RealImageType> DistanceFilterType;
    typename DistanceFilterType::Pointer distancer = DistanceFilterType::New();
    distancer->SetInput( thresholder->GetOutput() );
    distancer->SetBackgroundValue( 0 );
    distancer->SetInsideIsPositive( false );
    distancer->SetUseImageSpacing( this->m_UseImageSpacing );
    distancer->SetSquaredDistance( this->m_SquaredDistance );
    distancer->Update();

    ImageRegionIterator<OutputImageType> ItO( output, 
      output->GetLargestPossibleRegion() );
    ImageRegionIterator<RealImageType> ItD( distancer->GetOutput(), 
      distancer->GetOutput()->GetLargestPossibleRegion() );

    ItO.GoToBegin();
    ItD.GoToBegin();
    while ( !ItO.IsAtEnd() )
      {    
      OutputPixelType distances = ItO.Get();
      if ( this->m_NormalizeImage )
        {
        distances[label-1] = vcl_exp( -vnl_math_sqr( ItD.Get() / this->m_Sigma ) );
        }   
      else
        {
        distances[label-1] = ItD.Get();
        }            
      ItO.Set( distances );
      ++ItO;
      ++ItD;
      }
    } 


  this->GraftOutput( output );
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TOutputImage>
void
MultipleLabelToDistanceMapImageFilter<TInputImage, TOutputImage>
::PrintSelf(
  std::ostream& os, 
  Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Use image spacing: "
     << this->m_UseImageSpacing << std::endl;
  os << indent << "Squared distance: "
     << this->m_SquaredDistance << std::endl;
  os << indent << "Normalize: "
     << this->m_NormalizeImage << std::endl;
  os << indent << "Sigma: "
     << this->m_Sigma << std::endl;
}



}  //end namespace itk

#endif

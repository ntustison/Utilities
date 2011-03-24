#ifndef __itkBiasFieldEstimationBSplineFilter_txx
#define __itkBiasFieldEstimationBSplineFilter_txx

#include "itkBiasFieldEstimationBSplineFilter.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkPointSet.h"

namespace itk
{


template<class TInputImage, class TOutputImage>
BiasFieldEstimationBSplineFilter<TInputImage, TOutputImage>
::BiasFieldEstimationBSplineFilter() 
{
  this->m_IgnorePixelValue = NumericTraits<InputPixelType>::max();
  this->m_SplineOrder = 3;
  this->m_MinimumLevel = 0;
  this->m_MaximumLevel = 0;  
  
  this->m_ConfidenceImage = NULL;
}

template<class TInputImage, class TOutputImage>
BiasFieldEstimationBSplineFilter<TInputImage, TOutputImage>
::~BiasFieldEstimationBSplineFilter()
{
}

template<class TInputImage, class TOutputImage>
void
BiasFieldEstimationBSplineFilter<TInputImage, TOutputImage>
::GenerateData()
{
  typedef PointSet<VectorType, ImageDimension> 
      DeformationFieldPointSetType;
  typename DeformationFieldPointSetType::Pointer fieldPoints = 
             DeformationFieldPointSetType::New();    

  typedef Image<VectorType, ImageDimension> VectorImageType;
  typedef BSplineScatteredDataPointSetToImageFilter
    <DeformationFieldPointSetType, VectorImageType> BSplineFilterType;
  typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
  
  typename BSplineFilterType::WeightsContainerType::Pointer confidenceValues = 
    BSplineFilterType::WeightsContainerType::New();    

  ImageRegionConstIteratorWithIndex<InputImageType> 
    It( this->GetInput(), this->GetInput()->GetRequestedRegion() );

  itkDebugMacro( << "Extracting points from input deformation field. " )
  
  VectorType data;

  unsigned int N = 0;   
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    data[0] = It.Get();
            
    if ( data[0] != this->m_IgnorePixelValue )
      { 
      typename DeformationFieldPointSetType::PointType point;
      this->GetInput()->TransformIndexToPhysicalPoint( It.GetIndex(), point );
      
      fieldPoints->SetPointData( N, data );
      fieldPoints->SetPoint( N, point );

      if ( this->m_ConfidenceImage )
        {
        confidenceValues->InsertElement
          ( N, this->m_ConfidenceImage->GetPixel( It.GetIndex() ) );
        }
      N++;  
      }  
    }

  itkDebugMacro( << "Calculating the B-spline deformation field. " );
  
  typename OutputImageType::PointType origin;
  typename OutputImageType::SpacingType spacing;
  typename OutputImageType::SizeType size;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    origin[i] = this->GetInput( 0 )->GetOrigin()[i];
    spacing[i] = this->GetInput( 0 )->GetSpacing()[i];
    size[i] = this->GetInput( 0 )->GetLargestPossibleRegion().GetSize()[i];   
    }

  typename BSplineFilterType::ArrayType ncps;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    unsigned int spans = 1;
    for ( unsigned int j = 0; j < this->m_MinimumLevel; j++ )
      {
      spans *= 2;
      }
    ncps[i] = spans + this->m_SplineOrder;  
    }
  typename BSplineFilterType::ArrayType close;
  close.Fill( false );

  bspliner->SetOrigin( origin );
  bspliner->SetSpacing( spacing );
  bspliner->SetSize( size );
  bspliner->SetGenerateOutputImage( true );
  bspliner->SetNumberOfLevels( this->m_MaximumLevel-this->m_MinimumLevel+1 );
  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetNumberOfControlPoints( ncps );
  bspliner->SetCloseDimension( close );
  bspliner->SetInput( fieldPoints );
  if ( this->m_ConfidenceImage )
    {
    bspliner->SetPointWeights( confidenceValues );
    }
  bspliner->Update();
  
  this->GetOutput()->SetSpacing( bspliner->GetOutput()->GetSpacing() );
  this->GetOutput()->SetOrigin( bspliner->GetOutput()->GetOrigin() );
  this->GetOutput()->SetRegions( bspliner->GetOutput()->GetLargestPossibleRegion() );
  this->GetOutput()->Allocate();

  ImageRegionIterator<OutputImageType> 
    ItO( this->GetOutput(), this->GetOutput()->GetRequestedRegion() );
  ImageRegionConstIterator<VectorImageType> 
    ItB( bspliner->GetOutput(), bspliner->GetOutput()->GetRequestedRegion() );

  for ( ItB.GoToBegin(), ItO.GoToBegin(); !ItB.IsAtEnd(); ++ItB, ++ItO )
    {
    ItO.Set( static_cast<OutputPixelType>( ItB.Get()[0] ) ); 
    }
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TOutputImage>
void
BiasFieldEstimationBSplineFilter<TInputImage, TOutputImage>
::PrintSelf( std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Ignore pixel value: "
     << this->m_IgnorePixelValue << std::endl;
  os << indent << "Minimum level: "
     << this->m_MinimumLevel << std::endl;
  os << indent << "Maximum level: "
     << this->m_MaximumLevel << std::endl;
  os << indent << "Spline order: "
     << this->m_SplineOrder << std::endl;
}

} // end namespace itk

#endif

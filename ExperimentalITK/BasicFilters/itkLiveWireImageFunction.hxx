/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLiveWireImageFunction.hxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:52 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkLiveWireImageFunction_hxx
#define __itkLiveWireImageFunction_hxx

#include "itkLiveWireImageFunction.h"

#include "itkCastImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhood.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkZeroCrossingBasedEdgeDetectionImageFilter.h"

#include "vnl/vnl_math.h"

namespace itk
{

// Constructor
template<class TInputImage>
LiveWireImageFunction<TInputImage>
::LiveWireImageFunction()
{
  this->m_GradientMagnitudeWeight = 0.43;
  this->m_ZeroCrossingWeight = 0.43;
  this->m_GradientDirectionWeight = 0.14; 
  this->m_UseFaceConnectedness = true;
  this->m_UseImageSpacing = true;
  this->m_FindStepEdges = true;

  this->m_ZeroCrossingImage = NULL;
  this->m_MaskImage = NULL;

  this->m_InsidePixelValue = NumericTraits<MaskPixelType>::One;

  this->m_AnchorSeed.Fill( 0 );
}
 
// Destructor
template<class TInputImage>
LiveWireImageFunction<TInputImage>
::~LiveWireImageFunction()
{
}

/**
 * Initialize by setting the input image
 */
template <class TInputImage>
void
LiveWireImageFunction<TInputImage>
::SetInputImage( const InputImageType * ptr )
{
  Superclass::SetInputImage( ptr );

  /**
   * Generate images for generating the weights from the input image
   */ 
  typedef CastImageFilter<InputImageType, RealImageType> CasterType;
  typename CasterType::Pointer caster = CasterType::New();
  caster->SetInput( this->GetInputImage() );
  caster->Update();

  typename GradientFilterType::Pointer gradient = GradientFilterType::New();
  gradient->SetInput( this->GetInputImage() );
  gradient->SetUseImageSpacing( this->m_UseImageSpacing );
  gradient->Update();
  this->m_GradientImage = gradient->GetOutput();

  this->m_GradientMagnitudeImage = RealImageType::New();
  this->m_GradientMagnitudeImage->SetOrigin( 
    this->GetInputImage()->GetOrigin() );
  this->m_GradientMagnitudeImage->SetSpacing( 
    this->GetInputImage()->GetSpacing() );
  this->m_GradientMagnitudeImage->SetRegions( 
    this->GetInputImage()->GetRequestedRegion() );
  this->m_GradientMagnitudeImage->Allocate();

  ImageRegionIterator<GradientImageType> ItG( this->m_GradientImage, 
    this->m_GradientImage->GetRequestedRegion() );
  ImageRegionIterator<RealImageType> ItM( this->m_GradientMagnitudeImage, 
    this->m_GradientMagnitudeImage->GetRequestedRegion() );
  for ( ItG.GoToBegin(), ItM.GoToBegin(); !ItG.IsAtEnd(); ++ItG, ++ItM )
    {
    ItM.Set( ( ItG.Get() ).GetNorm() );
    } 

  typedef RescaleIntensityImageFilter<RealImageType, RealImageType> RescalerType;
  typename RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetInput( this->m_GradientMagnitudeImage );
  rescaler->SetOutputMinimum( NumericTraits<RealType>::Zero );
  rescaler->SetOutputMaximum( NumericTraits<RealType>::One );
  rescaler->Update();  
  this->m_RescaledGradientMagnitudeImage = rescaler->GetOutput();

  if ( !this->m_ZeroCrossingImage )
    {
    typedef ZeroCrossingBasedEdgeDetectionImageFilter<RealImageType, RealImageType> 
      ZeroCrossingFilterType;
    typename ZeroCrossingFilterType::Pointer 
      zeroCrossing = ZeroCrossingFilterType::New();
    zeroCrossing->SetInput( caster->GetOutput() );
    zeroCrossing->SetForegroundValue( NumericTraits<RealType>::One ); 
    zeroCrossing->SetBackgroundValue( NumericTraits<RealType>::Zero ); 
    zeroCrossing->Update();
    this->m_ZeroCrossingImage = zeroCrossing->GetOutput();
    }

  this->GeneratePathDirectionImage();
}

template <class TInputImage>
void
LiveWireImageFunction<TInputImage>
::GeneratePathDirectionImage()
{
  if ( !this->IsInsideBuffer( this->m_AnchorSeed ) )
    {
    itkExceptionMacro( "m_AnchorSeed is not inside buffer." );
    }

  /**
   * Initialize data structures
   */
  typedef Image<bool, ImageDimension> BooleanImageType;
  typename BooleanImageType::Pointer expanded = BooleanImageType::New();
  expanded->SetOrigin( this->GetInputImage()->GetOrigin() );
  expanded->SetSpacing( this->GetInputImage()->GetSpacing() );
  expanded->SetRegions( this->GetInputImage()->GetRequestedRegion() );
  expanded->Allocate();
  expanded->FillBuffer( false );
  expanded->SetPixel( this->m_AnchorSeed, true );

  this->m_PathDirectionImage = OffsetImageType::New();
  this->m_PathDirectionImage->SetOrigin( this->GetInputImage()->GetOrigin() );
  this->m_PathDirectionImage->SetSpacing( this->GetInputImage()->GetSpacing() );
  this->m_PathDirectionImage->SetRegions( this->GetInputImage()->GetRequestedRegion() );
  this->m_PathDirectionImage->Allocate();

  this->m_PathDirectionCostImage = RealImageType::New();
  this->m_PathDirectionCostImage->SetOrigin( this->GetInputImage()->GetOrigin() );
  this->m_PathDirectionCostImage->SetSpacing( this->GetInputImage()->GetSpacing() );
  this->m_PathDirectionCostImage->SetRegions( this->GetInputImage()->GetRequestedRegion() );
  this->m_PathDirectionCostImage->Allocate();
  this->m_PathDirectionCostImage->FillBuffer( 0 );
 
  typename ConstNeighborhoodIterator<InputImageType>::RadiusType radius;
  radius.Fill( 1 );
  ConstNeighborhoodIterator<BooleanImageType> It( radius, expanded, 
    expanded->GetRequestedRegion() );
  unsigned int numberOfNeighbors = 1;
  for ( unsigned int d = 0; d < ImageDimension; d++ )
    {
    numberOfNeighbors *= ( 2*radius[d] + 1 ); 
    }
  Neighborhood<RealType, ImageDimension> scaleFactors;
  scaleFactors.SetRadius( radius );
  for ( unsigned int n = 0; n < numberOfNeighbors; n++ )
    {
    scaleFactors[n] = 0;
    if ( n == static_cast<unsigned int>( 0.5 * numberOfNeighbors ) )
      {
      continue; 
      } 
    typename Neighborhood<RealType, ImageDimension>::OffsetType offset 
      = scaleFactors.GetOffset( n );
   
    bool isFaceConnected = true;
    unsigned int sumOffset = 0; 
    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      sumOffset += vnl_math_abs( offset[d] );
      if ( this->m_UseImageSpacing )
        {
        scaleFactors[n] += ( static_cast<RealType>( offset[d] * offset[d] ) 
          * this->GetInputImage()->GetSpacing()[d] 
          * this->GetInputImage()->GetSpacing()[d] );  
        }
      else
        {
        scaleFactors[n] += static_cast<RealType>( offset[d] * offset[d] );  
        }
      if ( sumOffset > 1 && this->m_UseFaceConnectedness )
        {
        isFaceConnected = false;
        break;
        }
      }      
    if ( !isFaceConnected )
      {
      scaleFactors[n] = 0;
      }
    if ( scaleFactors[n] > 0 )
      {
      scaleFactors[n] = vcl_sqrt( scaleFactors[n] );
      } 
    }   

  /**
   * Generate the PathDirectionImage with Dijkstra's algorithm
   */
  
  NodeType anchorNode;
  anchorNode.cost = 0.0;
  anchorNode.index = this->m_AnchorSeed;

  if ( this->m_MaskImage && 
       this->m_MaskImage->GetPixel( this->m_AnchorSeed ) 
         != this->m_InsidePixelValue )
    {
    itkWarningMacro( "The anchor seed is outside the user-defined mask region." ); 
    return;
    }    

  PriorityQueueType Q;
  Q.push( anchorNode );

  while ( !Q.empty() )
    {
    NodeType centerNode = Q.top(); 
    Q.pop();

    expanded->SetPixel( centerNode.index, true );
    It.SetLocation( centerNode.index ); 
    PointType centerPoint;
    this->GetInputImage()->TransformIndexToPhysicalPoint( 
      centerNode.index, centerPoint );

    typename GradientImageType::PixelType centerGradient 
      = this->m_GradientImage->GetPixel( centerNode.index );
    RealType centerNorm 
      = this->m_GradientMagnitudeImage->GetPixel( centerNode.index );

    for ( unsigned int n = 0; n < numberOfNeighbors; n++ )
      {
      if ( scaleFactors[n] == 0 )
        {
        continue;
        } 
      bool inBounds;
      bool isExpanded = It.GetPixel( n, inBounds );
      if ( isExpanded || !inBounds )
        {
        continue;
        }  

      IndexType index = It.GetIndex( n );
      NodeType neighborNode;
      neighborNode.index = index;

      if ( this->m_MaskImage && 
           this->m_MaskImage->GetPixel( index ) != this->m_InsidePixelValue )
        {
        continue;
        }    

      typename GradientImageType::PixelType neighborGradient 
        = this->m_GradientImage->GetPixel( index );
      RealType neighborNorm = this->m_GradientMagnitudeImage->GetPixel( index );
      
      RealType fz = 0.0;
      RealType fg = 0.0;
      RealType fd = 0.0;

      if ( this->m_ZeroCrossingWeight > 0 && this->m_ZeroCrossingImage )
        {
        fz = 1.0 - this->m_ZeroCrossingImage->GetPixel( index );
        } 
      if ( this->m_FindStepEdges )
        {
        if ( this->m_GradientMagnitudeWeight > 0.0 )
          {
          fg = 1.0 - this->m_RescaledGradientMagnitudeImage->GetPixel( index );
          }
        if ( this->m_GradientDirectionWeight > 0.0 
             && neighborNorm > 0 && centerNorm > 0 )
          {
          PointType neighborPoint;
          this->GetInputImage()->TransformIndexToPhysicalPoint( 
            neighborNode.index, neighborPoint );
          typename PointType::VectorType vector 
            = neighborPoint - centerPoint;
          RealType vectorNorm = vector.GetNorm();            

          RealType centerMin = vnl_math_min( centerGradient * vector, 
            centerGradient * -vector );            
          RealType neighborMin = vnl_math_min( neighborGradient * vector, 
            neighborGradient * -vector );            

          fd = 1.0 - ( vcl_acos( centerMin / ( centerNorm * vectorNorm ) ) +  
            vcl_acos( neighborMin / ( neighborNorm * vectorNorm ) ) )
            / vnl_math::pi;
          }
        }
      else
        {
        if ( this->m_GradientMagnitudeWeight > 0.0 )
          {
          fg = this->m_RescaledGradientMagnitudeImage->GetPixel( index );
          } 
        if ( this->m_GradientDirectionWeight > 0.0 
             && neighborNorm > 0 && centerNorm > 0 )
          {
          PointType neighborPoint;
          this->GetInputImage()->TransformIndexToPhysicalPoint( 
            neighborNode.index, neighborPoint );
          typename PointType::VectorType vector 
            = neighborPoint - centerPoint;
          RealType vectorNorm = vector.GetNorm();            

          RealType centerMin = vnl_math_min( centerGradient * vector, 
            centerGradient * -vector );            
          RealType neighborMin = vnl_math_min( neighborGradient * vector, 
            neighborGradient * -vector );            

          fd = 1.0 - ( vcl_acos( centerMin / ( centerNorm * vectorNorm ) ) + 
            vcl_acos( neighborMin / ( neighborNorm * vectorNorm ) ) )
            / vnl_math::pi;
          } 
        }

      neighborNode.cost = centerNode.cost + scaleFactors[n] *
        ( fg * this->m_GradientMagnitudeWeight
        + fz * this->m_ZeroCrossingWeight
        + fd * this->m_GradientDirectionWeight );

      Q.push( neighborNode ); 

      this->m_PathDirectionImage->SetPixel( neighborNode.index, 
        centerNode.index - neighborNode.index ); 

      this->m_PathDirectionCostImage->SetPixel( neighborNode.index, 
        neighborNode.cost ); 

      NodeType node = Q.top();
      while ( expanded->GetPixel( node.index ) )
        {
        Q.pop();
        node = Q.top();
        } 
      }  
    }  
}

template <class TInputImage>
typename LiveWireImageFunction<TInputImage>::OutputType::Pointer 
LiveWireImageFunction<TInputImage>
::EvaluateAtIndex( const IndexType &index ) const
{ 
  if ( !this->IsInsideBuffer( index ) )
    {
    itkWarningMacro( "Requested index is not inside buffer." );
    return NULL;
    }
  if ( this->m_MaskImage && 
       this->m_MaskImage->GetPixel( index ) 
         != this->m_InsidePixelValue )
    {
    itkWarningMacro( "The index is outside the user-defined mask region." ); 
    return NULL;
    }    

  typename OutputType::Pointer output = OutputType::New();
  output->Initialize();

  IndexType currentIndex = index;
  while ( currentIndex != this->m_AnchorSeed )
    {
    output->AddVertex( VertexType( currentIndex ) );
    currentIndex += this->m_PathDirectionImage->GetPixel( currentIndex );
    } 
  output->AddVertex( VertexType( currentIndex ) );

  return output;
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputImage>
void
LiveWireImageFunction<TInputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "AnchorSeed: " 
     << this->m_AnchorSeed << std::endl;
  os << indent << "GradientMagnitudeWeight: " 
     << this->m_GradientMagnitudeWeight << std::endl;
  os << indent << "GradientDirectionWeight: " 
     << this->m_GradientDirectionWeight << std::endl;
  os << indent << "ZeroCrossingWeight: " 
     << this->m_ZeroCrossingWeight << std::endl;
  os << indent << "UseImageSpacing: " 
     << this->m_UseImageSpacing << std::endl;
  os << indent << "UseFaceConnectedness: " 
     << this->m_UseFaceConnectedness << std::endl;

  os << indent << "MaskImage"
     << this->m_MaskImage << std::endl;
  os << indent << "InsidePixelValue"
     << this->m_InsidePixelValue << std::endl;
  
}
  
} // end namespace itk

#endif

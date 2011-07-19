/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBSplineVectorFlowImageFunction.hxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:51 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkBSplineVectorFlowImageFunction_hxx
#define _itkBSplineVectorFlowImageFunction_hxx

#include "itkBSplineVectorFlowImageFunction.h"

#include "itkAddImageFilter.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace itk
{
template <typename TVectorFlowField>
BSplineVectorFlowImageFunction<TVectorFlowField>
::BSplineVectorFlowImageFunction()
{
  this->m_ControlPointLattice = NULL;
  this->m_MaskImage = NULL;
  this->m_ForegroundValue = NumericTraits<typename MaskImageType::PixelType>::One;

  this->m_SplineOrder.Fill( 3 );
  this->m_NumberOfControlPoints.Fill( 4 );
  this->m_NumberOfLevels.Fill( 1 );

  this->m_CurrentTimePoint = 0;

  this->m_Size.Fill( 100 );
  this->m_Origin.Fill( 0 );
  this->m_Spacing.Fill( 1 );

}

template <typename TVectorFlowField>
typename BSplineVectorFlowImageFunction<TVectorFlowField>
::VectorType
BSplineVectorFlowImageFunction<TVectorFlowField>
::Evaluate( const PointType &point ) const
{

/*
  if ( this->m_ControlPointLattice == NULL )
    {
    itkExceptionMacro( "Error: The vector flow field has not yet been instantiated." );
    } 
*/

  typedef BSplineControlPointImageFilter<ControlPointLatticeType, VectorFlowFieldType> BSplinerType;
  typename BSplinerType::Pointer bspliner = BSplinerType::New();
  ArrayType close;
  close.Fill( false );

  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetCloseDimension( close );
  bspliner->SetInput( this->m_ControlPointLattice );
  bspliner->SetOrigin( this->m_Origin );
  bspliner->SetSize( this->m_Size );
  bspliner->SetSpacing( this->m_Spacing );
  
  typename VectorFlowFieldType::PointType pt;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    pt[i] = point[i];
    } 
  VectorType vector;
  bspliner->EvaluateAtPoint( pt, vector );

  return vector;
}

template <typename TVectorFlowField>
void BSplineVectorFlowImageFunction<TVectorFlowField>
::UpdateWithDeformationFieldAtTimePoint( 
    typename DeformationFieldType::Pointer deformationField, RealType timePoint )
{

  ImageRegionIteratorWithIndex<DeformationFieldType> It( deformationField,
    deformationField->GetLargestPossibleRegion() );
  
  typename PointSetType::Pointer points = PointSetType::New();
  points->Initialize();

  unsigned long count = 0;
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if ( !this->m_MaskImage || 
         ( this->m_MaskImage->GetPixel( It.GetIndex() ) == this->m_ForegroundValue ) )
      {
      typename DeformationFieldType::PointType dPoint;
      deformationField->TransformIndexToPhysicalPoint( It.GetIndex(), dPoint );

      typename PointSetType::PointType bPoint;  
      for ( unsigned int i = 0; i < ImageDimension-1; i++ )
        {
        bPoint[i] = dPoint[i];
        }       
      bPoint[ImageDimension-1] = timePoint;
      
      points->SetPoint( count, bPoint );
      points->SetPointData( count, It.Get() );
      count++;
      } 
    }

  typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
  ArrayType close;
  close.Fill( false );

  bspliner->SetInput( points );
  bspliner->SetOrigin( this->m_Origin );
  bspliner->SetSize( this->m_Size );
  bspliner->SetSpacing( this->m_Spacing );
  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetCloseDimension( close );
  bspliner->SetNumberOfLevels( this->m_NumberOfLevels );              
  bspliner->SetNumberOfControlPoints( this->m_NumberOfControlPoints );       
  bspliner->SetGenerateOutputImage( false );
  bspliner->Update();

  if ( !this->m_ControlPointLattice )
    {
    this->m_ControlPointLattice = ControlPointLatticeType::New();
    this->m_ControlPointLattice = bspliner->GetPhiLattice();
    }
  else
    { 
    typedef AddImageFilter<ControlPointLatticeType, ControlPointLatticeType, 
      ControlPointLatticeType> AdderType;
    typename AdderType::Pointer adder = AdderType::New();
    adder->SetInput1( bspliner->GetPhiLattice() );
    adder->SetInput2( this->m_ControlPointLattice );
    adder->Update(); 
    
    this->m_ControlPointLattice = adder->GetOutput(); 
    }  
 
}

template <typename TVectorFlowField>
typename BSplineVectorFlowImageFunction<TVectorFlowField>
::DeformationFieldType::Pointer
BSplineVectorFlowImageFunction<TVectorFlowField>
::GenerateDeformationFieldAtTimePoint( RealType t )
{
  typedef BSplineControlPointImageFilter<ControlPointLatticeType, VectorFlowFieldType> BSplinerType;
  typename BSplinerType::Pointer bspliner = BSplinerType::New();
  ArrayType close;
  close.Fill( false );

  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetCloseDimension( close );
  bspliner->SetInput( this->m_ControlPointLattice );
  bspliner->SetOrigin( this->m_Origin );
  bspliner->SetSize( this->m_Size );
  bspliner->SetSpacing( this->m_Spacing );

  typename ControlPointLatticeType::PointType P;
  P.Fill( 1e10 );
  P[ImageDimension-1] = t;

  typename VectorFlowFieldType::Pointer vectorFlowField = 
    VectorFlowFieldType::New();
  vectorFlowField = bspliner->GenerateOutputImageAt( P );

  typedef ExtractImageFilter<VectorFlowFieldType, DeformationFieldType> ExtracterType;
  typename ExtracterType::Pointer extracter = ExtracterType::New();
  extracter->SetInput( vectorFlowField );

  typename ExtracterType::InputImageRegionType region;
  typename ExtracterType::InputImageSizeType regionSize;
  typename ExtracterType::InputImageIndexType regionIndex;
   
  for ( unsigned int i = 0; i < ImageDimension-1; i++ )
    {
    regionSize[i] = vectorFlowField->GetLargestPossibleRegion().GetSize()[i];
    regionIndex[i] = vectorFlowField->GetLargestPossibleRegion().GetIndex()[i];
    }
  regionSize[ImageDimension-1] = 0;
  regionIndex[ImageDimension-1] = 0;  
  region.SetSize( regionSize );
  region.SetIndex( regionIndex );

  extracter->SetExtractionRegion( region );
  extracter->Update();
  return extracter->GetOutput();
};

template <typename TVectorFlowField>
void
BSplineVectorFlowImageFunction<TVectorFlowField>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
}
  
} // end namespace itk

#endif

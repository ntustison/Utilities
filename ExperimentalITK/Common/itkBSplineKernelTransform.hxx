/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBSplineKernelTransform.hxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:20:03 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkBSplineKernelTransform_hxx
#define _itkBSplineKernelTransform_hxx
#include "itkBSplineKernelTransform.h"

namespace itk
{

/**
 *
 */
template <class TScalarType, unsigned int NDimensions>
BSplineKernelTransform<TScalarType, NDimensions>::
BSplineKernelTransform()
{
  this->m_SourceLandmarks = PointSetType::New();
  this->m_TargetLandmarks = PointSetType::New();
  this->m_Displacements   = VectorSetType::New();
  
  this->m_BSplineTransform = BSplineTransformType::New();
  this->m_BSplineTransformIsCalculated = false;
  
  this->m_SplineOrder = 3;
  this->m_NumberOfControlPoints.Fill( this->m_SplineOrder + 1 );
  this->m_NumberOfLevels = 3;
  
  this->m_DomainExpansionFactor = 0.5;
  this->m_Origin.Fill( 0.0 );
  this->m_Spacing.Fill( 0.0 );
  this->m_Size.Fill( 0 );
}

/**
 *
 */
template <class TScalarType, unsigned int NDimensions>
BSplineKernelTransform<TScalarType, NDimensions>::
~BSplineKernelTransform()
{
}

/**
 *
 */
template <class TScalarType, unsigned int NDimensions>
void
BSplineKernelTransform<TScalarType, NDimensions>
::GenerateBSplineTransformation()
{

  typename BSplinePointSetType::Pointer bsplinePoints
    = BSplinePointSetType::New();
  bsplinePoints->Initialize();
  
  typename BSplinePointSetType::PointsContainerIterator ItB 
    = bsplinePoints->GetPoints()->Begin();
  typename BSplinePointSetType::PointDataContainerIterator ItV 
    = bsplinePoints->GetPointData()->Begin();
  
  PointsConstIterator ItS 
    = this->m_SourceLandmarks->GetPoints()->Begin();
  PointsConstIterator ItT 
    = this->m_TargetLandmarks->GetPoints()->Begin();
    
  while ( ItB != this->m_SourceLandmarks->GetPoints()->End() )
    {
//    ItB->Value() = ItS->Value();
//    ItV->Value() = ItT->Value() - ItS->Value(); 
    
    ++ItB;
    ++ItV;
    ++ItS;
    ++ItT;
    }
    
  /**
   * Determine the transformation domain
   */
 bool useExpansionFactor = false;  
 for ( unsigned int d = 0; d < SpaceDimension; d++ )
   {  
   if ( this->m_Size[d] == 0 || this->m_Spacing[d] == 0.0 )
     {
     useExpansionFactor = true;
     }  
   } 
    
  if ( useExpansionFactor )
    {
    InputPointType sourceMinPoint 
      = this->m_SourceLandmarks->GetBoundingBox()->GetMinimum();  
    InputPointType sourceMaxPoint 
      = this->m_SourceLandmarks->GetBoundingBox()->GetMaximum();  
    InputPointType targetMinPoint 
      = this->m_TargetLandmarks->GetBoundingBox()->GetMinimum();  
    InputPointType targetMaxPoint 
      = this->m_TargetLandmarks->GetBoundingBox()->GetMaximum();  
     
    InputPointType minPoint;
    minPoint.Fill( NumericTraits<ScalarType>::max() );
    InputPointType maxPoint;
    maxPoint.Fill( NumericTraits<ScalarType>::min() );
    for ( unsigned int d = 0; d < SpaceDimension; d++ )
      {
      minPoint[d] = vnl_math_min( sourceMinPoint[d], targetMinPoint[d] ); 
      maxPoint[d] = vnl_math_min( sourceMaxPoint[d], targetMaxPoint[d] ); 
      } 

    /**
      * The size and spacing are irrelevant other than 
      * the fact that they define the transformation domain.
      */
    this->m_Size.Fill( 100 );
  
    for ( unsigned int d = 0; d < SpaceDimension; d++ )
      {
      /**
        * Expand the bounding box to ensure coverage
        */
      ScalarType expandedMinimum = minPoint[d] 
        - this->m_DomainExpansionFactor * ( maxPoint[d] - minPoint[d] );
      ScalarType expandedMaximum = maxPoint[d] 
        + this->m_DomainExpansionFactor * ( maxPoint[d] - minPoint[d] );
      this->m_Origin[d] = expandedMinimum;
      this->m_Spacing[d] = ( expandedMaximum - expandedMinimum ) 
        / static_cast<ScalarType>( this->m_Size[d] - 1 );
      }
    }  
    
  this->m_BSplineTransform->SetInput( bsplinePoints );
  this->m_BSplineTransform->SetOrigin( this->m_Origin );
  this->m_BSplineTransform->SetSpacing( this->m_Spacing );
  this->m_BSplineTransform->SetSize( this->m_Size );
  this->m_BSplineTransform->SetNumberOfLevels( this->m_NumberOfLevels );
  this->m_BSplineTransform->SetSplineOrder( this->m_SplineOrder );
  this->m_BSplineTransform->SetNumberOfControlPoints( 
    this->m_NumberOfControlPoints );
  this->m_BSplineTransform->SetGenerateOutputImage( false );
  this->m_BSplineTransform->Update();
 
  this->m_BSplineTransformIsCalculated = true; 
} 

/**
 *
 */
template <class TScalarType, unsigned int NDimensions>
typename BSplineKernelTransform<TScalarType, NDimensions>::OutputPointType
BSplineKernelTransform<TScalarType, NDimensions>
::TransformPoint( const InputPointType& thisPoint ) const
{
  if ( !this->m_BSplineTransformIsCalculated )
    {
    this->GenerateBSplineTransformation(); 
    } 

  OutputPointType result;

  InputVectorType vector;
  typename BSplinePointSetType::PointType bsplinePoint;
  bsplinePoint.CastFrom( thisPoint );
  this->m_BSplineTransform->EvaluateAtPoint( bsplinePoint, vector );

  result = thisPoint + vector;

  return result;
}

// Compute the Jacobian in one position 
template <class TScalarType, unsigned int NDimensions>
const typename BSplineKernelTransform<TScalarType,NDimensions>::JacobianType & 
BSplineKernelTransform< TScalarType,NDimensions>::
GetJacobian( const InputPointType &thisPoint ) const
{
  if ( !this->m_BSplineTransformIsCalculated )
    {
    this->GenerateBSplineTransformation(); 
    } 

  typename BSplinePointSetType::PointType bsplinePoint;
  bsplinePoint.CastFrom( thisPoint );

  typename BSplineTransformType::GradientType J;
  this->m_BSplineTransform->EvaluateGradientAtPoint( bsplinePoint, J );

  this->m_Jacobian.SetSize( J.Rows(), J.Cols() );

  for ( unsigned int i = 0; i < J.Rows(); i++ )
    {
    for ( unsigned int j = 0; j < J.Cols(); j++ )
      {
      this->m_Jacobian(i, j) = J(i, j);
      if ( i == j )
        {
        this->m_Jacobian(i, j) += 1.0; 
        } 
      } 
    }

  return this->m_Jacobian;
}

// Return the Displacements 
template <class TScalarType, unsigned int NDimensions>
typename BSplineKernelTransform<TScalarType,NDimensions>::VectorSetType * 
BSplineKernelTransform< TScalarType,NDimensions>::
GetDisplacements()
{
  if ( !this->m_BSplineTransformIsCalculated )
    {
    this->GenerateBSplineTransformation(); 
    } 

  if ( this->m_Displacements->Size() == 0 )
    {
    typename VectorSetType::Iterator ItD 
      = this->m_Displacements->Begin();
    typename VectorSetType::Iterator ItS 
      = this->m_SourceLandmarks->GetPoints()->Begin();

    while ( ItD != this->m_Displacements->End() )
      {
      InputVectorType vector;
      typename BSplinePointSetType::PointType point;
      for ( unsigned int d = 0; d < SpaceDimension; d++ )
        {
        point[d] = ItS.Value()[d]; 
        }
      this->m_BSplineTransform->EvaluateAtPoint( point, vector );
    
      ItD.Value() = vector;
      
      ++ItD;
      ++ItS;
      }
    }
     
  return this->m_Displacements.GetPointer();  
}


template <class TScalarType, unsigned int NDimensions>
void
BSplineKernelTransform<TScalarType, NDimensions>::
PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  if( this->m_SourceLandmarks )
    {
    os << indent << "SourceLandmarks: " << std::endl;
    this->m_SourceLandmarks->Print(os,indent.GetNextIndent());
    }
  if( this->m_TargetLandmarks )
    {
    os << indent << "TargetLandmarks: " << std::endl;
    this->m_TargetLandmarks->Print(os,indent.GetNextIndent());
    }
  if( this->m_Displacements )
    {
    os << indent << "Displacements: " << std::endl;
    this->m_Displacements->Print(os,indent.GetNextIndent());
    }
  os << indent << "Stiffness: " << this->m_Stiffness << std::endl;
}
} // namespace itk

#endif

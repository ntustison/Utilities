/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBSplineVectorFlowImageFunction.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:51 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBSplineVectorFlowImageFunction_h
#define __itkBSplineVectorFlowImageFunction_h

#include "itkImageFunction.h"

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkContinuousIndex.h"
#include "itkPointSet.h"

namespace itk
{
/** \class BSplineVectorFlowImageFunction
 *
 */
template <typename TVectorFlowField>
class ITK_EXPORT BSplineVectorFlowImageFunction :
    public ImageFunction<TVectorFlowField, 
             typename TVectorFlowField::PixelType,
             typename TVectorFlowField::PixelType::ValueType>
{
public:
  /** Standard class typedefs. */
  typedef BSplineVectorFlowImageFunction                   Self;
  typedef ImageFunction<TVectorFlowField, 
           typename TVectorFlowField::PixelType,
           typename TVectorFlowField::PixelType::
                     ValueType>                            Superclass;
  typedef SmartPointer<Self>                               Pointer;
  typedef SmartPointer<const Self>                         ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods) */
  itkTypeMacro( BSplineVectorFlowImageFunction, ImageFunction );

  /** The dimensionality of the input deformation field. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TVectorFlowField::ImageDimension );

  /** 
   * Extract some information from the image types.  
   */
  typedef TVectorFlowField                                 VectorFlowFieldType;
  typedef typename TVectorFlowField::PixelType             VectorType;
  typedef typename VectorType::ValueType                   RealType;

  /**
   * B-spline related typedefs 
   */ 
  typedef PointSet<VectorType, 
    itkGetStaticConstMacro( ImageDimension )>              PointSetType;
  typedef Image<VectorType, 
    itkGetStaticConstMacro( ImageDimension ) - 1>          DeformationFieldType;
  typedef BSplineScatteredDataPointSetToImageFilter
    <PointSetType, VectorFlowFieldType>                    BSplineFilterType;
  typedef typename BSplineFilterType::PointDataImageType   ControlPointLatticeType;
  typedef typename BSplineFilterType::ArrayType            ArrayType;  

  /**
   * typedefs to define the domain of the vector flow field
   */
  typedef typename VectorFlowFieldType::SizeType           SizeType;
  typedef typename VectorFlowFieldType::PointType          OriginType;
  typedef typename VectorFlowFieldType::SpacingType        SpacingType; 

  /**
   * other typedefs
   */
  typedef typename Superclass::PointType                   PointType;
  typedef typename Superclass::IndexType                   IndexType;
  typedef typename Superclass::ContinuousIndexType         ContinuousIndexType;
  typedef Image<unsigned int, 
    itkGetStaticConstMacro( ImageDimension ) - 1>          MaskImageType;
  

  itkSetMacro( Size, SizeType );
  itkGetConstMacro( Size, SizeType );

  itkSetMacro( Origin, OriginType );
  itkGetConstMacro( Origin, OriginType );

  itkSetMacro( Spacing, SpacingType );
  itkGetConstMacro( Spacing, SpacingType );

  itkSetMacro( CurrentTimePoint, RealType );
  itkGetConstMacro( CurrentTimePoint, RealType );

  itkSetMacro( SplineOrder, ArrayType );
  itkGetConstMacro( SplineOrder, ArrayType ); 

  itkSetMacro( NumberOfLevels, ArrayType );
  itkGetConstMacro( NumberOfLevels, ArrayType ); 

  itkSetMacro( NumberOfControlPoints, ArrayType );
  itkGetConstMacro( NumberOfControlPoints, ArrayType ); 

  itkSetMacro( MaskImage, typename MaskImageType::Pointer );
  itkGetConstMacro( MaskImage, typename MaskImageType::Pointer );

  itkSetMacro( ForegroundValue, typename MaskImageType::PixelType );
  itkGetConstMacro( ForegroundValue, typename MaskImageType::PixelType );

   

  void UpdateWithDeformationFieldAtTimePoint( typename DeformationFieldType::Pointer, RealType );

  typename DeformationFieldType::Pointer
    GenerateDeformationFieldAtTimePoint( RealType );

  /**
   * Evaluate at a point in the vector flow field
   */
  virtual VectorType Evaluate( const PointType& ) const;

  virtual VectorType EvaluateAtIndex( const IndexType & index ) const
    {
    PointType point;
    this->GetInputImage()->TransformIndexToPhysicalPoint( index, point );
    return this->Evaluate( point );     
    } 
  virtual VectorType EvaluateAtContinuousIndex( const ContinuousIndexType & index ) const
    {
    PointType point;
    this->GetInputImage()->TransformContinuousIndexToPhysicalPoint( index, point );
    return this->Evaluate( point );     
    }  


protected:
  BSplineVectorFlowImageFunction();
  virtual ~BSplineVectorFlowImageFunction() {}

  void PrintSelf ( std::ostream& os, Indent indent ) const;

private:
  BSplineVectorFlowImageFunction(const Self&); //purposely not implemented
  VectorType operator=(const Self&); //purposely not implemented


  SizeType                                                 m_Size;                                 
  OriginType                                               m_Origin;                                 
  SpacingType                                              m_Spacing;                                  

  ArrayType                                                m_SplineOrder;
  ArrayType                                                m_NumberOfLevels;
  ArrayType                                                m_NumberOfControlPoints;

  RealType                                                 m_CurrentTimePoint;

  typename MaskImageType::Pointer                          m_MaskImage;
  typename MaskImageType::PixelType                        m_ForegroundValue;

  typename ControlPointLatticeType::Pointer                m_ControlPointLattice;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBSplineVectorFlowImageFunction.txx"
#endif

#endif

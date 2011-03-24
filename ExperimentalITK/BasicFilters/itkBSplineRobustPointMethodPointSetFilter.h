/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBSplineRobustPointMethodPointSetFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:51 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBSplineRobustPointMethodPointSetFilter_h
#define __itkBSplineRobustPointMethodPointSetFilter_h

#include "itkPointSetToImageFilter.h"

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkPointSet.h"
#include "itkVariableLengthVector.h"
#include "itkVariableSizeMatrix.h"
#include "itkVector.h"

namespace itk
{

/** \class BSplineRobustPointMethodPointSetFilter.h
 * \brief point set filter.
 */

template <class TPointSet, 
          class TOutputImage = Image< 
                Vector<float, GetPointSetDimension<TPointSet>::PointDimension>, 
                GetPointSetDimension<TPointSet>::PointDimension> >
class ITK_EXPORT BSplineRobustPointMethodPointSetFilter 
: public PointSetToImageFilter<TPointSet, TOutputImage>
{
public:
  typedef BSplineRobustPointMethodPointSetFilter              Self;
  typedef PointSetToImageFilter<TPointSet, TOutputImage>      Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Extract dimension from output image. */
  itkStaticConstMacro( Dimension, unsigned int,
                       TOutputImage::ImageDimension );
        
  /** Image typedef support. */
  typedef TPointSet                                           InputPointSetType;
  typedef TOutputImage                                        VectorFieldType;
  typedef typename VectorFieldType::PixelType                 VectorType;
  typedef typename VectorType::ValueType                      RealType;

  /** Other typedef */
  typedef VariableSizeMatrix<RealType>                        MatrixType;
  typedef VariableLengthVector<RealType>                      OutlierVectorType;

  /** B-spline typedefs */
  typedef PointSet<VectorType, 
    itkGetStaticConstMacro( Dimension )>                      PointSetType;
  typedef typename PointSetType::PointType                    PointType;
  typedef BSplineScatteredDataPointSetToImageFilter
    <PointSetType, VectorFieldType>                           BSplineFitterType;
  typedef typename BSplineFitterType::ArrayType               ArrayType;
  typedef VectorFieldType                                     ControlPointLatticeType;
  typedef typename BSplineFitterType::WeightsContainerType    WeightsContainerType;

  /** Helper functions */

  itkSetMacro( SplineOrder, unsigned int );
  itkGetConstMacro( SplineOrder, unsigned int );
  
  itkSetMacro( NumberOfLevels, unsigned int );
  itkGetConstMacro( NumberOfLevels, unsigned int );
  
  itkSetMacro( NumberOfControlPoints, ArrayType );
  itkGetConstMacro( NumberOfControlPoints, ArrayType );

  itkSetClampMacro( AnnealingRate, RealType, 0.9, 0.99 );
  itkGetConstMacro( AnnealingRate, RealType );

  itkSetClampMacro( InitialTemperature, RealType, 0, NumericTraits<RealType>::max() );
  itkGetConstMacro( InitialTemperature, RealType );

  itkSetClampMacro( FinalTemperature, RealType, 0, NumericTraits<RealType>::max() );
  itkGetConstMacro( FinalTemperature, RealType );

  itkSetMacro( NumberOfIterationsPerTemperature, unsigned int );
  itkGetConstMacro( NumberOfIterationsPerTemperature, unsigned int );

  itkBooleanMacro( SolveSimplerLeastSquaresProblem );
  itkSetMacro( SolveSimplerLeastSquaresProblem, bool );
  itkGetConstMacro( SolveSimplerLeastSquaresProblem, bool );

  itkBooleanMacro( UseBoundingBox );
  itkSetMacro( UseBoundingBox, bool );
  itkGetConstMacro( UseBoundingBox, bool );

  itkGetConstMacro( VPoints, typename InputPointSetType::Pointer );
  itkGetConstMacro( ControlPointLattice, typename ControlPointLatticeType::Pointer );


protected:
  BSplineRobustPointMethodPointSetFilter();
  virtual ~BSplineRobustPointMethodPointSetFilter();
  void PrintSelf( std::ostream& os, Indent indent ) const;

  void GenerateData();

private:
  void Initialize();
  void CalculateInitialAndFinalTemperatures();
  void UpdateCorrespondenceMatrix();
  void NormalizeCorrespondenceMatrix();
  void UpdateTransformation();
 

  void VisualizeCurrentState();

private:
  BSplineRobustPointMethodPointSetFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  typename InputPointSetType::Pointer                        m_VPoints;

  MatrixType                                                 m_CorrespondenceMatrix;
  OutlierVectorType                                          m_OutlierRow;
  OutlierVectorType                                          m_OutlierColumn;
  typename InputPointSetType::PointType                      m_OutlierPointX;
  typename InputPointSetType::PointType                      m_OutlierPointV;

  RealType                                                   m_AnnealingRate;
  RealType                                                   m_InitialTemperature;
  RealType                                                   m_FinalTemperature;
  RealType                                                   m_CurrentTemperature;
  unsigned int                                               m_NumberOfIterationsPerTemperature;
  bool                                                       m_SolveSimplerLeastSquaresProblem;
  bool                                                       m_UseBoundingBox;
    
  typename ControlPointLatticeType::Pointer                  m_ControlPointLattice;

  ArrayType                                                  m_NumberOfControlPoints;
  unsigned int                                               m_SplineOrder;
  unsigned int                                               m_NumberOfLevels;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBSplineRobustPointMethodPointSetFilter.txx"
#endif

#endif


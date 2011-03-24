/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkThinPlateSplineRobustPointMethodPointSetFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:53 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkThinPlateSplineRobustPointMethodPointSetFilter_h
#define __itkThinPlateSplineRobustPointMethodPointSetFilter_h

#include "itkPointSetToImageFilter.h"

#include "itkPointSet.h"
#include "itkKernelTransform.h"
#include "itkVariableLengthVector.h"
#include "itkVariableSizeMatrix.h"
#include "itkVector.h"

namespace itk
{

/** \class ThinPlateSplineRobustPointMethodPointSetFilter.h
 * \brief point set filter.
 */

template <class TPointSet, 
          class TOutputImage = Image< 
                Vector<float, GetPointSetDimension<TPointSet>::PointDimension>, 
                GetPointSetDimension<TPointSet>::PointDimension> >
class ITK_EXPORT ThinPlateSplineRobustPointMethodPointSetFilter 
: public PointSetToImageFilter<TPointSet, TOutputImage>
{
public:
  typedef ThinPlateSplineRobustPointMethodPointSetFilter      Self;
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
  typedef TOutputImage                                        OutputImageType;
  typedef TOutputImage                                        VectorFieldType;
  typedef typename VectorFieldType::PixelType                 VectorType;
  typedef typename VectorType::ValueType                      RealType;

  /** Other typedef */
  typedef VariableSizeMatrix<RealType>                        MatrixType;
  typedef VariableLengthVector<RealType>                      OutlierVectorType;

  /** thin-plate spline typedefs */
  typedef KernelTransform<RealType,                      
    itkGetStaticConstMacro( Dimension )>                      TransformType;
  typedef typename TransformType::PointSetType                PointSetType;
  typedef typename PointSetType::PointType                    PointType;

  /** Helper functions */

  itkSetClampMacro( AnnealingRate, RealType, 0.9, 0.99 );
  itkGetConstMacro( AnnealingRate, RealType );

  itkSetClampMacro( InitialLambda, RealType, 0, NumericTraits<RealType>::max() );
  itkGetConstMacro( InitialLambda, RealType );

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

protected:
  ThinPlateSplineRobustPointMethodPointSetFilter();
  virtual ~ThinPlateSplineRobustPointMethodPointSetFilter();
  void PrintSelf( std::ostream& os, Indent indent ) const;

  void GenerateData();

private:
  void Initialize();
  void CalculateInitialAndFinalTemperatures();
  void UpdateCorrespondenceMatrix();
  void UpdateTransformation();
 
  void VisualizeCurrentState();

private:
  ThinPlateSplineRobustPointMethodPointSetFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  typename InputPointSetType::Pointer                        m_VPoints;

  MatrixType                                                 m_CorrespondenceMatrix;
  OutlierVectorType                                          m_OutlierRow;
  OutlierVectorType                                          m_OutlierColumn;
  typename InputPointSetType::PointType                      m_OutlierPointX;
  typename InputPointSetType::PointType                      m_OutlierPointV;

  RealType                                                   m_InitialLambda;
  RealType                                                   m_CurrentLambda;
  RealType                                                   m_AnnealingRate;
  RealType                                                   m_InitialTemperature;
  RealType                                                   m_FinalTemperature;
  RealType                                                   m_CurrentTemperature;
  unsigned int                                               m_NumberOfIterationsPerTemperature;
  bool                                                       m_SolveSimplerLeastSquaresProblem;
  bool                                                       m_UseBoundingBox;
    
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkThinPlateSplineRobustPointMethodPointSetFilter.txx"
#endif

#endif


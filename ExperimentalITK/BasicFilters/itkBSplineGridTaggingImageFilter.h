/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBSplineGridTaggingImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:51 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBSplineGridTaggingImageFilter_h
#define __itkBSplineGridTaggingImageFilter_h

#include "itkImageToImageFilter.h"

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkPointSet.h"
#include "itkVector.h"

#include "vnl/vnl_sparse_matrix.h"
#include "vnl/vnl_vector.h"

namespace itk
{

/** \class BSplineGridTaggingImageFilter.h
 * \brief Image filter.
 */

template <class TControlPointLattice, 
          class TCandidatePointImage, 
          class TOutputImage = TControlPointLattice>
class ITK_EXPORT BSplineGridTaggingImageFilter 
: public ImageToImageFilter<TControlPointLattice, TOutputImage>
{
public:
  typedef BSplineGridTaggingImageFilter                       Self;
  typedef ImageToImageFilter<TControlPointLattice, 
                             TOutputImage>                    Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Extract dimension from input image. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TControlPointLattice::ImageDimension );
        
  /** Image typedef support. */
  typedef TControlPointLattice                                ControlPointLatticeType;
  typedef TCandidatePointImage                                CandidatePointImageType;
  typedef TCandidatePointImage                                LabelImageType;
  typedef TOutputImage                                        OutputImageType;

  /** Other typedef */
  typedef float                                               RealType;
  typedef Image<RealType, 
    itkGetStaticConstMacro( ImageDimension )>                 RealImageType; 
  typedef FixedArray<unsigned int, 
    itkGetStaticConstMacro( ImageDimension )>                 ArrayType; 
  typedef FixedArray<unsigned int, 1>                         LabelType;
  typedef PointSet<LabelType, 
    itkGetStaticConstMacro( ImageDimension )>                 LabelPointSetType;
  typedef vnl_sparse_matrix<RealType>                         SparseMatrixType;
  typedef vnl_vector<RealType>                                OutlierArrayType;

  /** B-spline typedefs */
  typedef Vector<RealType, 1>                                 ScalarType; 
  typedef Image<ScalarType, 
    itkGetStaticConstMacro( ImageDimension )>                 ScalarFieldType; 
  typedef PointSet<ScalarType, 
    itkGetStaticConstMacro( ImageDimension )>                 PointSetType;
  typedef typename PointSetType::PointType                    PointType;
  typedef BSplineScatteredDataPointSetToImageFilter
    <PointSetType, ScalarFieldType>                           BSplineFitterType;
  typedef ScalarFieldType                                     ScalarControlPointLatticeType;
  typedef typename BSplineFitterType::WeightsContainerType    WeightsContainerType;

  /** Helper functions */

  void SetInput1( const TControlPointLattice *);
  void SetInput2( const TCandidatePointImage *);

  itkSetClampMacro( WhichParametricDimension, unsigned int, 
    0, itkGetStaticConstMacro( ImageDimension ) - 1 );
  itkGetConstMacro( WhichParametricDimension, unsigned int );

  itkSetMacro( NumberOfTagLabels, unsigned int );
  itkGetConstMacro( NumberOfTagLabels, unsigned int );

  itkSetMacro( BoundingBoxMinimum, typename CandidatePointImageType::PointType );
  itkGetConstMacro( BoundingBoxMinimum, typename CandidatePointImageType::PointType );

  itkSetMacro( BoundingBoxMaximum, typename CandidatePointImageType::PointType );
  itkGetConstMacro( BoundingBoxMaximum, typename CandidatePointImageType::PointType );

  itkSetMacro( FirstTagLocation, RealType );
  itkGetConstMacro( FirstTagLocation, RealType );

  itkSetMacro( InitialTagSpacing, RealType );
  itkGetConstMacro( InitialTagSpacing, RealType );

  itkSetMacro( SplineOrder, unsigned int );
  itkGetConstMacro( SplineOrder, unsigned int );
  
  itkSetMacro( NumberOfLevels, unsigned int );
  itkGetConstMacro( NumberOfLevels, unsigned int );
  
  itkSetMacro( NumberOfControlPoints, ArrayType );
  itkGetConstMacro( NumberOfControlPoints, ArrayType );

  itkSetClampMacro( AnnealingRate, RealType, 0.9, 0.99 );
  itkGetConstMacro( AnnealingRate, RealType );

  itkSetMacro( NumberOfIterationsPerTemperature, unsigned int );
  itkGetConstMacro( NumberOfIterationsPerTemperature, unsigned int );

  itkBooleanMacro( SolveSimplerLeastSquaresProblem );
  itkSetMacro( SolveSimplerLeastSquaresProblem, bool );
  itkGetConstMacro( SolveSimplerLeastSquaresProblem, bool );

protected:
  BSplineGridTaggingImageFilter();
  virtual ~BSplineGridTaggingImageFilter();
  void PrintSelf( std::ostream& os, Indent indent ) const;

  void GenerateData();

private:
  void InitializeLabelPointSets();
  void UpdateCorrespondenceMatrix();
  void NormalizeCorrespondenceMatrix();
  void UpdateTransformation();

private:
  BSplineGridTaggingImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  typename LabelPointSetType::Pointer                        m_XPoints;
  typename LabelPointSetType::Pointer                        m_VPoints;

  SparseMatrixType                                           m_CorrespondenceMatrix;
  OutlierArrayType                                           m_OutlierRow;
  OutlierArrayType                                           m_OutlierColumn;

  RealType                                                   m_AnnealingRate;
  RealType                                                   m_InitialTemperature;
  RealType                                                   m_FinalTemperature;
  RealType                                                   m_CurrentTemperature;
  unsigned int                                               m_NumberOfIterationsPerTemperature;
  bool                                                       m_SolveSimplerLeastSquaresProblem;
    
  typename ScalarControlPointLatticeType::Pointer            m_ScalarControlPointLattice;

  unsigned int                                               m_NumberOfTagLabels;
  unsigned int                                               m_WhichParametricDimension;
  typename CandidatePointImageType::PointType                m_BoundingBoxMinimum;
  typename CandidatePointImageType::PointType                m_BoundingBoxMaximum;
  RealType                                                   m_InitialTagSpacing;
  RealType                                                   m_FirstTagLocation;

  ArrayType                                                  m_NumberOfControlPoints;
  unsigned int                                               m_SplineOrder;
  unsigned int                                               m_NumberOfLevels;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBSplineGridTaggingImageFilter.txx"
#endif

#endif


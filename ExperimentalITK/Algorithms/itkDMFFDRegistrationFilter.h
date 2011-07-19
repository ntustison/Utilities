/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDMFFDRegistrationFilter.h,v $
  Language:  C++
  Date:      $Date: 2009/04/10 16:12:07 $
  Version:   $Revision: 1.2 $

  Copyright ( c ) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkDMFFDRegistrationFilter_h_
#define _itkDMFFDRegistrationFilter_h_

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkInterpolateImageFunction.h"
#include "itkNeighborhoodIterator.h"
#include "itkAvantsPDEDeformableRegistrationFunction.h"
#include "itkPointSet.h"
#include "itkRecursiveMultiResolutionPyramidImageFilter.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkVector.h"
#include "itkVectorContainer.h"

#include "vnl/vnl_vector.h"

#include <iostream>
#include <string>

#define itkSetElementObjectMacro( name,type ) \
  virtual void Set##name ( type* _arg, const unsigned int _i ) \
  { \
    itkDebugMacro( "setting element " << _i << " of " #name " to " << _arg ); \
    if (this->m_##name[_i] != _arg) \
      { \
      this->m_##name[_i] = _arg; \
      this->Modified(); \
      } \
  }
#define itkGetElementObjectMacro( name,type ) \
  virtual type * Get##name ( const unsigned int _i = 0 ) \
  { \
    itkDebugMacro( "returning element " << _i << " of " << #name " of " << this->m_##name  ); \
    return this->m_##name[_i].GetPointer(); \
  }

namespace itk {
/** \class DMFFDRegistrationFilter
    \brief DMFFD Image registration filter.
     \par
*/

template<class TMovingImage, class TFixedImage = TMovingImage,
  class TWarpedImage = TMovingImage>
class ITK_EXPORT DMFFDRegistrationFilter
: public ImageToImageFilter<TFixedImage, TWarpedImage>
{
public:
  typedef DMFFDRegistrationFilter                              Self;
  typedef ImageToImageFilter<TMovingImage, TWarpedImage>       Superclass;
  typedef SmartPointer<Self>                                   Pointer;
  typedef SmartPointer<const Self>                             ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ) */
  itkTypeMacro( DMFFDRegistrationFilter, ImageToImageFilter );

  typedef TMovingImage                                         MovingImageType;
  typedef TFixedImage                                          FixedImageType;
  typedef TWarpedImage                                         WarpedImageType;

  typedef typename FixedImageType::PixelType                   PixelType;
  typedef typename FixedImageType::SizeType                    SizeType;
  
  /** Dimensionality of input and output data is assumed to be the same. */
  itkStaticConstMacro( ImageDimension, unsigned int,
    FixedImageType::ImageDimension );

  typedef float                                                RealType;
  typedef Image<RealType,
          itkGetStaticConstMacro( ImageDimension )>            RealImageType;
  typedef Image<RealType,
          itkGetStaticConstMacro( ImageDimension )>            WeightImageType;
  typedef FixedArray<unsigned int,
          itkGetStaticConstMacro( ImageDimension )>            ArrayType;
  typedef VectorContainer<unsigned, ArrayType>                 ArrayContainerType;
  typedef Array<unsigned int>                                  ResizableUIntArrayType;
  typedef Array<RealType>                                      ResizableRealArrayType;

  /** Typedefs for the deformation field */
  typedef Vector<RealType,
    itkGetStaticConstMacro( ImageDimension )>                  VectorType;
  typedef Image<VectorType,
    itkGetStaticConstMacro( ImageDimension )>                  DeformationFieldType;

  typedef NeighborhoodIterator<DeformationFieldType>           NeighborhoodIteratorType;
  typedef typename NeighborhoodIteratorType::RadiusType        RadiusType;

  typedef RecursiveMultiResolutionPyramidImageFilter
          <FixedImageType, FixedImageType>                     FixedImagePyramidType;
  typedef typename FixedImagePyramidType::ScheduleType         FixedScheduleType;
  typedef RecursiveMultiResolutionPyramidImageFilter
          <MovingImageType, MovingImageType>                   MovingImagePyramidType;
  typedef typename FixedImagePyramidType::ScheduleType         MovingScheduleType;
          
  typedef RecursiveMultiResolutionPyramidImageFilter
          <WeightImageType, WeightImageType>                   WeightImagePyramidType;

  /** Typedef support for the interpolation function */
  typedef InterpolateImageFunction<FixedImageType, double>     ImageInterpolatorType;
  typedef typename ImageInterpolatorType::Pointer              ImageInterpolatorPointer;

  /** Typedefs for image metrics */
  typedef AvantsPDEDeformableRegistrationFunction
     <FixedImageType, MovingImageType, DeformationFieldType>   PDEDeformableMetricType;
  typedef typename PDEDeformableMetricType::Pointer            PDEDeformableMetricPointer;
  typedef typename PDEDeformableMetricType::RadiusType         MetricRadiusType;

  /** Typedefs for B-spline filter */
  typedef PointSet<VectorType,
    itkGetStaticConstMacro( ImageDimension )>                  PointSetType;
  typedef BSplineScatteredDataPointSetToImageFilter
    <PointSetType, DeformationFieldType>                       BSplineFilterType;
  typedef typename BSplineFilterType::PointDataImageType       ControlPointLatticeType;
  typedef typename ControlPointLatticeType::Pointer            ControlPointLatticePointer;

  /** Main functions */
  void RunRegistration()
    { this->Update(); }

  /** Set/Get functions ( variables are explained below ). */

  itkSetMacro( WeightImage, typename WeightImageType::Pointer );
  itkGetConstMacro( WeightImage, typename WeightImageType::Pointer );

  itkSetMacro( InitialDeformationFieldControlPoints,
    typename ControlPointLatticeType::Pointer );
  itkGetConstMacro( InitialDeformationFieldControlPoints,
    typename ControlPointLatticeType::Pointer );

  itkGetConstMacro( TotalDeformationFieldControlPoints,
    typename ControlPointLatticeType::Pointer );

  itkSetMacro( SplineOrder, unsigned int );
  itkGetConstMacro( SplineOrder, unsigned int );

  itkSetMacro( MaximumNumberOfIterations, ResizableUIntArrayType );
  itkGetConstMacro( MaximumNumberOfIterations, ResizableUIntArrayType );

  itkSetMacro( LineSearchMaximumIterations, unsigned int );
  itkGetConstMacro( LineSearchMaximumIterations, unsigned int );

  itkSetMacro( GradientScalingFactor, ResizableRealArrayType );
  itkGetConstMacro( GradientScalingFactor, ResizableRealArrayType );

  itkSetMacro( EnforceDiffeomorphism, bool );
  itkGetConstMacro( EnforceDiffeomorphism, bool );
  itkBooleanMacro( EnforceDiffeomorphism );

  itkSetMacro( EmploySteepestDescent, bool );
  itkGetConstMacro( EmploySteepestDescent, bool );
  itkBooleanMacro( EmploySteepestDescent );

  itkSetMacro( LineSearchMaximumStepSize, RealType );
  itkGetConstMacro( LineSearchMaximumStepSize, RealType );

  itkSetMacro( WhichGradient, unsigned short );
  itkGetConstMacro( WhichGradient, unsigned short );

  itkSetMacro( MinimumJacobian, RealType );
  itkGetConstMacro( MinimumJacobian, RealType );

  itkSetMacro( MaximumJacobian, RealType );
  itkGetConstMacro( MaximumJacobian, RealType );

  itkSetMacro( InteriorPenaltyParameter, RealType );
  itkGetConstMacro( InteriorPenaltyParameter, RealType );

  itkSetMacro( InitialMeshResolution, ArrayType );
  itkGetConstMacro( InitialMeshResolution, ArrayType );

  itkSetMacro( DoubleMeshResolutionAtEachLevel, bool );
  itkGetConstMacro( DoubleMeshResolutionAtEachLevel, bool );
  itkBooleanMacro( DoubleMeshResolutionAtEachLevel );

  itkSetMacro( FixedImageShrinkFactors, ArrayType );
  itkGetConstMacro( FixedImageShrinkFactors, ArrayType );

  itkSetMacro( MovingImageShrinkFactors, ArrayType );
  itkGetConstMacro( MovingImageShrinkFactors, ArrayType );

  /** Set/Get the image metric.  ( PDEDeformableRegistrationFunction ) */
  itkSetElementObjectMacro( PDEDeformableMetric, PDEDeformableMetricType );
  itkGetElementObjectMacro( PDEDeformableMetric, PDEDeformableMetricType );

  /**
   * Print to screen the course of optimization.
   */
  itkSetMacro( Verbose, bool );
  itkGetConstMacro( Verbose, bool );
  itkBooleanMacro( Verbose );

protected :
  /** de/constructor */
  DMFFDRegistrationFilter();
  ~DMFFDRegistrationFilter();

  void PrintSelf( std::ostream& os, Indent indent ) const;
  void GenerateData();

private :
  DMFFDRegistrationFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  void InitializeImages();

  RealType EvaluateGradientFieldOverImageRegion();
  RealType EvaluateMetricOverImageRegion( RealType );

  void IterativeSolve();  // Conjugate gradient descent

  RealType EvaluateEnergyForLineSearch( RealType );
  RealType FindBracketingTriplet( RealType*, RealType*, RealType* );
  void LineMinimization( RealType*, RealType* );

  void BrentSearch( RealType, RealType, RealType, RealType, RealType*,
    RealType* );

private :

  /**
   * Image variables
   */
  ArrayType                             m_FixedImageShrinkFactors;
  ArrayType                             m_MovingImageShrinkFactors;
  typename WeightImageType::Pointer     m_WeightImage;

  typename MovingImageType::Pointer     m_CurrentMovingImage[2];
  typename FixedImageType::Pointer      m_CurrentFixedImage[2];
  typename WeightImageType::Pointer     m_CurrentWeightImage;
  
  FixedScheduleType                     m_FixedPyramidSchedule;
  MovingScheduleType                    m_MovingPyramidSchedule;
  std::vector<SizeType>                 m_PyramidLevelImageSizes;

  /**
   * B-spline variables
   */
  unsigned int                          m_SplineOrder;
  ArrayType                             m_InitialMeshResolution;
  bool                                  m_DoubleMeshResolutionAtEachLevel;

  /**
   * Optimization variables
   */
  ResizableUIntArrayType                m_MaximumNumberOfIterations;
  ResizableRealArrayType                m_GradientScalingFactor;
  unsigned int                          m_LineSearchMaximumIterations;
  RealType                              m_LineSearchMaximumStepSize;
  unsigned int                          m_CurrentLevel;
  int                                   m_CurrentIteration;

  unsigned int                          m_NumberOfLevels;
  bool                                  m_EmploySteepestDescent;

  RealType                              m_MinimumJacobian;
  RealType                              m_MaximumJacobian;
  bool                                  m_EnforceDiffeomorphism;
  RealType                              m_InteriorPenaltyParameter;

  unsigned short                        m_WhichGradient;

  /**
   * Image metric variables
   */
  PDEDeformableMetricPointer            m_PDEDeformableMetric[2];
  bool                                  m_MaximizeMetric[2];
  MetricRadiusType                      m_MetricRadius[2];
  ImageInterpolatorPointer              m_ImageInterpolator;

  /**
   * Control point lattices
   */
  ControlPointLatticePointer            m_InitialDeformationFieldControlPoints;
  ControlPointLatticePointer            m_TotalDeformationFieldControlPoints;
  ControlPointLatticePointer            m_CurrentDeformationFieldControlPoints;
  ControlPointLatticePointer            m_GradientFieldControlPoints;
  
  /**
   * Other variables
   */
  bool                                  m_Verbose; 
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDMFFDRegistrationFilter.hxx"
#endif

#endif

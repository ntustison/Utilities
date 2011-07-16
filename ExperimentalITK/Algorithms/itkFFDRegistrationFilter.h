/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFFDRegistrationFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:13:41 $
  Version:   $Revision: 1.1.1.1 $

  Copyright ( c ) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkFFDRegistrationFilter_h_
#define _itkFFDRegistrationFilter_h_

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkInterpolateImageFunction.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
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

/** Set/Get built-in type to handle elements of arrays; */
#define itkSetElementMacro( name,type ) \
  void Set##name ( const type _arg, const unsigned int _i ) \
  { \
    itkDebugMacro( "setting element " << _i << " of " #name " to " << _arg ); \
    if ( this->m_##name[_i] != _arg ) \
      { \
      this->m_##name[_i] = _arg; \
      this->Modified(  ); \
      } \
  }
#define itkSetAllElementsMacro( name,type ) \
  void Set##name ( const type _arg ) \
  { \
    itkDebugMacro( "setting the entire vector to " << _arg ); \
    this->m_##name.fill( _arg ); \
    this->Modified(  ); \
  }
#define itkGetElementConstMacro( name,type ) \
  virtual type Get##name ( const unsigned int _i = 0 ) const \
  { \
    itkDebugMacro( "returning element " << _i << " of " << #name " of " << this->m_##name  ); \
    return this->m_##name[_i]; \
  }
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

#define itkGetConstSTDVectorMacro(name,type) \
  virtual type Get##name () const \
  { \
    return this->m_##name; \
  }


namespace itk {
/** \class FFDRegistrationFilter
    \brief FFD Image registration filter.
     \par
*/

template<class TMovingImage, class TFixedImage, class TWarpedImage = TFixedImage>
class ITK_EXPORT FFDRegistrationFilter
: public ImageToImageFilter<TFixedImage, TWarpedImage>
{
public:
  typedef FFDRegistrationFilter                                Self;
  typedef ImageToImageFilter<TMovingImage, TWarpedImage>       Superclass;
  typedef SmartPointer<Self>                                   Pointer;
  typedef SmartPointer<const Self>                             ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ) */
  itkTypeMacro( FFDRegistrationFilter, ImageToImageFilter );

  typedef TMovingImage                                         MovingImageType;
  typedef TFixedImage                                          FixedImageType;
  typedef TWarpedImage                                         WarpedImageType;

  typedef typename FixedImageType::PixelType                   PixelType;
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
  typedef Vector<RealType,
          itkGetStaticConstMacro( ImageDimension ) + 1>        LandmarkType;
  typedef VectorContainer<unsigned, LandmarkType>              LandmarkContainerType;

  /** Typedefs for the deformation field */
  typedef Vector<RealType,
    itkGetStaticConstMacro( ImageDimension )>                  VectorType;
  typedef Image<VectorType,
    itkGetStaticConstMacro( ImageDimension )>                  DeformationFieldType;

  typedef NeighborhoodIterator<DeformationFieldType>           NeighborhoodIteratorType;
  typedef ShapedNeighborhoodIterator<DeformationFieldType>     ShapedNeighborhoodIteratorType;
  typedef typename NeighborhoodIteratorType::RadiusType        RadiusType;

  typedef RecursiveMultiResolutionPyramidImageFilter
                          <FixedImageType, FixedImageType>     FixedImagePyramidType;
  typedef RecursiveMultiResolutionPyramidImageFilter
                          <MovingImageType, MovingImageType>   MovingImagePyramidType;
  typedef RecursiveMultiResolutionPyramidImageFilter
                  <WeightImageType, WeightImageType>           WeightImagePyramidType;

  /** Typedef support for the interpolation function */
  typedef InterpolateImageFunction<FixedImageType,
                                   double>                     ImageInterpolatorType;

  /** Typedefs for image metrics */
  typedef AvantsPDEDeformableRegistrationFunction
     <FixedImageType, MovingImageType, DeformationFieldType>   PDEDeformableMetricType;
  typedef typename PDEDeformableMetricType::RadiusType         MetricRadiusType;

  /** Typedefs for B-spline filter */
  typedef PointSet<VectorType,
    itkGetStaticConstMacro( ImageDimension )>                  PointSetType;
  typedef BSplineScatteredDataPointSetToImageFilter
    <PointSetType, DeformationFieldType>                       BSplineFilterType;
  typedef typename BSplineFilterType::PointDataImageType       ControlPointLatticeType;
  typedef typename Statistics
     ::MersenneTwisterRandomVariateGenerator                   RandomizerType;

  /** Main functions */
  void RunRegistration()
    { this->Update(); }

  /** Set/Get functions ( variables are explained below ). */

  itkSetMacro( WeightImage, typename WeightImageType::Pointer );
  itkGetConstMacro( WeightImage, typename WeightImageType::Pointer );

  itkSetMacro( InitialControlPointLattice, typename ControlPointLatticeType::Pointer );
  itkGetConstMacro( InitialControlPointLattice, typename ControlPointLatticeType::Pointer );

  itkGetConstMacro( TotalDeformationFieldControlPoints,
    typename ControlPointLatticeType::Pointer );

  itkSetMacro( FixedLandmarkContainer, typename LandmarkContainerType::Pointer );
  itkGetConstMacro( FixedLandmarkContainer, typename LandmarkContainerType::Pointer );
  itkSetMacro( MovingLandmarkContainer, typename LandmarkContainerType::Pointer );
  itkGetConstMacro( MovingLandmarkContainer, typename LandmarkContainerType::Pointer );

  itkSetMacro( SplineOrder, unsigned int );
  itkGetConstMacro( SplineOrder, unsigned int );

  itkSetMacro( Directionality, ArrayType );
  itkGetConstMacro( Directionality, ArrayType );

  itkSetElementMacro( MaximumNumberOfIterations, unsigned int );
  itkSetAllElementsMacro( MaximumNumberOfIterations, unsigned int );
  itkGetElementConstMacro( MaximumNumberOfIterations, unsigned int );
  itkSetMacro( MaximumNumberOfIterations, vnl_vector<unsigned int> );
  itkGetConstMacro( MaximumNumberOfIterations, vnl_vector<unsigned int> );

  void SetNumberOfLevels( unsigned int );
  itkGetConstMacro( NumberOfLevels, unsigned int );

  itkSetMacro( LineSearchMaximumIterations, unsigned int );
  itkGetConstMacro( LineSearchMaximumIterations, unsigned int );

  itkSetMacro( InitializeWithLandmarks, bool );
  itkGetConstMacro( InitializeWithLandmarks, bool );
  itkBooleanMacro( InitializeWithLandmarks );

  itkSetMacro( EnforceDiffeomorphism, bool );
  itkGetConstMacro( EnforceDiffeomorphism, bool );
  itkBooleanMacro( EnforceDiffeomorphism );

  itkSetMacro( EmploySteepestDescent, bool );
  itkGetConstMacro( EmploySteepestDescent, bool );
  itkBooleanMacro( EmploySteepestDescent );

  itkSetMacro( WhichGradient, unsigned short );
  itkGetConstMacro( WhichGradient, unsigned short );

  itkSetMacro( NumberOfAffineNeighborhoodSamplesPerIteration, unsigned int );
  itkGetConstMacro( NumberOfAffineNeighborhoodSamplesPerIteration, unsigned int );

  itkSetMacro( AffineNeighborhoodRadius, unsigned int );
  itkGetConstMacro( AffineNeighborhoodRadius, unsigned int );

  itkSetClampMacro( AffineNeighborhoodSamplingDensity, RealType, 0, 1 );
  itkGetConstMacro( AffineNeighborhoodSamplingDensity, RealType );

  itkSetMacro( LandmarkWeighting, RealType );
  itkGetConstMacro( LandmarkWeighting, RealType );

  itkSetMacro( MinimumJacobian, RealType );
  itkGetConstMacro( MinimumJacobian, RealType );

  itkSetMacro( MaximumJacobian, RealType );
  itkGetConstMacro( MaximumJacobian, RealType );

  itkSetMacro( InteriorPenaltyParameter, RealType );
  itkGetConstMacro( InteriorPenaltyParameter, RealType );

  itkSetMacro( MeshResolution, ArrayType );
  itkGetConstMacro( MeshResolution, ArrayType );

  itkSetMacro( DoubleMeshResolutionAtEachLevel, bool );
  itkGetConstMacro( DoubleMeshResolutionAtEachLevel, bool );
  itkBooleanMacro( DoubleMeshResolutionAtEachLevel );

  itkSetMacro( FixedImageShrinkFactors, ArrayType );
  itkGetConstMacro( FixedImageShrinkFactors, ArrayType );

  itkSetMacro( MovingImageShrinkFactors, ArrayType );
  itkGetConstMacro( MovingImageShrinkFactors, ArrayType );

  /** Set/Get the image interpolator. */
  itkSetObjectMacro( ImageInterpolator, ImageInterpolatorType );
  itkGetObjectMacro( ImageInterpolator, ImageInterpolatorType );

  /** Set/Get the image metric.  ( PDEDeformableRegistrationFunction ) */
  itkSetElementObjectMacro( PDEDeformableMetric, PDEDeformableMetricType );
  itkGetElementObjectMacro( PDEDeformableMetric, PDEDeformableMetricType );

  itkSetElementMacro( MetricRadius, MetricRadiusType );
  itkGetElementConstMacro( MetricRadius, MetricRadiusType );

  itkSetElementMacro( MaximizeMetric, bool );
  itkGetElementConstMacro( MaximizeMetric, bool );

  itkSetStringMacro( FilePrefix );
  itkGetStringMacro( FilePrefix );

  itkSetMacro( Prolificacy, bool );
  itkGetConstMacro( Prolificacy, bool );
  itkBooleanMacro( Prolificacy );

  itkGetConstSTDVectorMacro( GradientComputationTimes, std::vector<RealType> );
  itkGetConstSTDVectorMacro( NumberOfGradientPoints, std::vector<unsigned int> );
  itkGetConstSTDVectorMacro( EnergyValues, std::vector<RealType> );
  itkGetConstSTDVectorMacro( LevelNumbers, std::vector<unsigned int> );

  RealType CalculateLandmarkError();

protected :
  /** de/constructor */
  FFDRegistrationFilter();
  ~FFDRegistrationFilter();

  void PrintSelf( std::ostream& os, Indent indent ) const;
  void GenerateData();

private :
  FFDRegistrationFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  void InitializeImages();

  RealType EvaluateGradientFieldOverImageRegion();
  void EvaluateGradientFieldOverImageRegionForAffineTransform();
  RealType EvaluateMetricOverImageRegion( RealType );
  void InitializeDeformationFieldWithLandmarks();


  void IterativeSolve();  // Conjugate gradient descent

  RealType EvaluateEnergyForLineSearch( RealType );
  void FindBracketingTriplet( RealType*, RealType*, RealType* );
  void LineMinimization( RealType*, RealType* );

  void BrentSearch( RealType, RealType, RealType, RealType*, RealType* );
  void BruteForceSearch( RealType, RealType, RealType, RealType*, RealType* );
  void GoldenSectionSearch( RealType, RealType, RealType, RealType*, RealType* );

private :

  /** Filter Outputs */
  typename ControlPointLatticeType::Pointer         m_TotalDeformationFieldControlPoints;
  typename ControlPointLatticeType::Pointer         m_CurrentDeformationFieldControlPoints;
  typename ControlPointLatticeType::Pointer         m_GradientFieldControlPoints;
  typename ControlPointLatticeType::Pointer         m_InitialControlPointLattice;

  typename WeightImageType::Pointer                 m_WeightImage;

  typename MovingImageType::Pointer                 m_CurrentMovingImage[2];
  typename FixedImageType::Pointer                  m_CurrentFixedImage[2];
  typename WeightImageType::Pointer                 m_CurrentWeightImage;

  typename FixedImagePyramidType::ScheduleType      m_FixedPyramidSchedule;
  typename MovingImagePyramidType::ScheduleType     m_MovingPyramidSchedule;
  std::vector<typename FixedImageType::SizeType>    m_PyramidLevelImageSizes;

  /** Set/get variables */
  unsigned int                                      m_NumberOfLevels;
  ArrayType                                         m_FixedImageShrinkFactors;
  ArrayType                                         m_MovingImageShrinkFactors;
  unsigned int                                      m_LineSearchMaximumIterations;
  unsigned int                                      m_SplineOrder;
  ArrayType                                         m_Directionality;
  unsigned short                                    m_WhichGradient;

  ArrayType                                         m_MeshResolution;
  bool                                              m_DoubleMeshResolutionAtEachLevel;
  vnl_vector<unsigned int>                          m_MaximumNumberOfIterations;

  unsigned int                                      m_CurrentLevel;
  typename LandmarkContainerType::Pointer           m_FixedLandmarkContainer;
  typename LandmarkContainerType::Pointer           m_MovingLandmarkContainer;
  RealType                                          m_LandmarkWeighting;
  RealType                                          m_MinimumJacobian;
  RealType                                          m_MaximumJacobian;
  RealType                                          m_InteriorPenaltyParameter;

  bool                                              m_EmploySteepestDescent;
  bool                                              m_EnforceDiffeomorphism;
  bool                                              m_InitializeWithLandmarks;
  typename ImageInterpolatorType::Pointer           m_ImageInterpolator;

  typename PDEDeformableMetricType::Pointer         m_PDEDeformableMetric[2];
  bool                                              m_MaximizeMetric[2];
  MetricRadiusType                                  m_MetricRadius[2];

  std::vector<RealType>                             m_GradientComputationTimes;
  std::vector<unsigned int>                         m_NumberOfGradientPoints;
  std::vector<RealType>                             m_EnergyValues;
  std::vector<unsigned int>                         m_LevelNumbers;

  std::string                                       m_FilePrefix;
  bool                                              m_Prolificacy;

  unsigned int                                      m_NumberOfAffineNeighborhoodSamplesPerIteration;
  RealType                                          m_AffineNeighborhoodSamplingDensity;
  unsigned int                                      m_AffineNeighborhoodRadius;
  typename RandomizerType::Pointer                  m_Randomizer;


};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFFDRegistrationFilter.hxx"
#endif

#endif

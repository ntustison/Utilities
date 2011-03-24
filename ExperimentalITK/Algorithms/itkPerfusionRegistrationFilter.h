/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPerfusionRegistrationFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:13:41 $
  Version:   $Revision: 1.1.1.1 $

  Copyright ( c ) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef _itkPerfusionRegistrationFilter_h_
#define _itkPerfusionRegistrationFilter_h_

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkInterpolateImageFunction.h"
#include "itkNeighborhoodIterator.h"
#include "itkPDEDeformableRegistrationFunction.h"
#include "itkPointSet.h"
#include "itkRecursiveMultiResolutionPyramidImageFilter.h"
#include "itkVector.h"
#include "itkVectorContainer.h"

#include "vnl/vnl_vector.h"

#include <iostream>
#include <string>
#include <vector>

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

namespace itk {
/** \class PerfusionRegistrationFilter
    \brief FFD Image registration filter.
     \par
*/

template<class TImage, class TWarpedImage = TImage>
class ITK_EXPORT PerfusionRegistrationFilter
: public ImageToImageFilter<TImage, TWarpedImage>
{
public:
  typedef PerfusionRegistrationFilter                              Self;
  typedef ImageToImageFilter<TImage, TWarpedImage>             Superclass;
  typedef SmartPointer<Self>                                   Pointer;
  typedef SmartPointer<const Self>                             ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ) */
  itkTypeMacro( PerfusionRegistrationFilter, ImageToImageFilter );

  typedef TImage                                               ImageType;
  typedef TWarpedImage                                         WarpedImageType;

  typedef typename ImageType::PixelType                        PixelType;
  /** Dimensionality of input and output data is assumed to be the same. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       ImageType::ImageDimension );

  typedef float                                                RealType;
  typedef Image<RealType,
          itkGetStaticConstMacro( ImageDimension )>            RealImageType;
  typedef Image<RealType,
          itkGetStaticConstMacro( ImageDimension ) + 1>        TimeDependentRealImageType;
  typedef FixedArray<unsigned int,
          itkGetStaticConstMacro( ImageDimension )>            ArrayType;
  typedef VectorContainer<unsigned, ArrayType>                 ArrayContainerType;

  /** Typedefs for the deformation field */
  typedef Vector<RealType,
          itkGetStaticConstMacro( ImageDimension )>            VectorType;
  typedef Image<VectorType,
          itkGetStaticConstMacro( ImageDimension )>            DeformationFieldType;
  typedef NeighborhoodIterator<DeformationFieldType>           NeighborhoodIteratorType;
  typedef typename NeighborhoodIteratorType::RadiusType        RadiusType;

  typedef RecursiveMultiResolutionPyramidImageFilter
                          <ImageType, ImageType>               ImagePyramidType;

  /** Typedef support for the interpolation function */
  typedef InterpolateImageFunction<ImageType, double>          ImageInterpolatorType;

  /** Typedefs for image metrics */
  typedef PDEDeformableRegistrationFunction
    <ImageType, ImageType, DeformationFieldType>               PDEDeformableMetricType;
  typedef typename PDEDeformableMetricType::RadiusType         MetricRadiusType;

  /** Typedefs for B-spline filter */
  typedef PointSet<VectorType,
    itkGetStaticConstMacro( ImageDimension ) + 1>              PointSetType;
  typedef Image<VectorType,
          itkGetStaticConstMacro( ImageDimension ) + 1>        TimeDependentFieldType;
  typedef BSplineScatteredDataPointSetToImageFilter
    <PointSetType, TimeDependentFieldType>                     BSplineFilterType;
  typedef typename BSplineFilterType::PointDataImageType       ControlPointLatticeType;
  typedef typename BSplineFilterType::WeightsContainerType     WeightsContainerType;

  /** Main functions */
  void RunRegistration()
    { this->Update(); }

  /** Set/Get functions ( variables are explained below ). */

  itkSetMacro( SplineOrder, unsigned int );
  itkGetConstMacro( SplineOrder, unsigned int );

  itkSetElementMacro( MaximumNumberOfIterations, unsigned int );
  itkSetAllElementsMacro( MaximumNumberOfIterations, unsigned int );
  itkGetElementConstMacro( MaximumNumberOfIterations, unsigned int );
  itkSetMacro( MaximumNumberOfIterations, vnl_vector<unsigned int> );
  itkGetConstMacro( MaximumNumberOfIterations, vnl_vector<unsigned int> );

  itkSetMacro( LineSearchMaximumIterations, unsigned int );
  itkGetConstMacro( LineSearchMaximumIterations, unsigned int );

  itkSetMacro( MeshResolution, ArrayType );
  itkGetConstMacro( MeshResolution, ArrayType );

  itkSetMacro( TemporalMeshResolution, unsigned int );
  itkGetConstMacro( TemporalMeshResolution, unsigned int );

  itkSetMacro( ImageShrinkFactors, ArrayType );
  itkGetConstMacro( ImageShrinkFactors, ArrayType );

  /** Set/Get the image interpolator. */
  itkSetObjectMacro( ImageInterpolator, ImageInterpolatorType );
  itkGetObjectMacro( ImageInterpolator, ImageInterpolatorType );

  /** Set/Get the image metric.  ( PDEDeformableRegistrationFunction ) */
  itkSetObjectMacro( PDEDeformableMetric, PDEDeformableMetricType );
  itkGetObjectMacro( PDEDeformableMetric, PDEDeformableMetricType );

  itkSetMacro( MetricRadius, MetricRadiusType );
  itkGetConstMacro( MetricRadius, MetricRadiusType );

  itkSetMacro( MaximizeMetric, bool );
  itkGetConstMacro( MaximizeMetric, bool );

  itkSetMacro( TemporalOrigin, RealType );
  itkGetConstMacro( TemporalOrigin, RealType );

  itkSetMacro( TemporalEnd, RealType );
  itkGetConstMacro( TemporalEnd, RealType );

  itkSetMacro( GradientScalingFactor, RealType );
  itkGetConstMacro( GradientScalingFactor, RealType );

  itkSetMacro( WrapTime, bool );
  itkGetConstMacro( WrapTime, bool );
  itkBooleanMacro( WrapTime );

  itkSetMacro( TemporalSplineOrder, unsigned int );
  itkGetConstMacro( TemporalSplineOrder, unsigned int );

  void SetNthTimePoint( unsigned int i, RealType t )
    {
    itkDebugMacro( "setting element " << i << " of m_TimePoints to " << t );
    if ( this->m_TimePoints.size() <= i || this->m_TimePoints.empty() )
      {
      this->m_TimePoints.resize( i + 1 );
      }
    this->m_TimePoints[i] = t;
    this->Modified();
    }
  RealType GetNthTimePoint( unsigned int i )
    {
    if ( this->m_TimePoints.size() > i )
      {
      return this->m_TimePoints[i];
      }
    }
  void SetNthImageFileName( unsigned int i, std::string s )
    {
    itkDebugMacro( "setting element " << i << " of m_ImageFileNames to " << s.c_str() );
    if ( this->m_ImageFileNames.size() <= i || this->m_ImageFileNames.empty() )
      {
      this->m_ImageFileNames.resize( i + 1 );
      }
    this->m_ImageFileNames[i] = s;
    this->Modified();
    }
  void SetNthImageFileName( unsigned int i, const char *c )
    {
    itkDebugMacro( "setting element " << i << " of m_ImageFileNames to " << c );
    if ( this->m_ImageFileNames.size() <= i || this->m_ImageFileNames.empty() )
      {
      this->m_ImageFileNames.resize( i + 1 );
      }
    this->m_ImageFileNames[i] = std::string( c );
    this->Modified();
    }
  const char* GetNthImageFileName( unsigned int i )
    {
    if ( this->m_ImageFileNames.size() > i )
      {
      return this->m_ImageFileNames[i].c_str();
      }
    }

  void GenerateOutputAtTimePoint( unsigned int );

protected :
  /** de/constructor */
  PerfusionRegistrationFilter();
  ~PerfusionRegistrationFilter();

  void PrintSelf( std::ostream& os, Indent indent ) const;
  void GenerateData();

private :

  PerfusionRegistrationFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  void InitializeImages();

  void IterativeSolve();

  void EvaluateGradientFieldOverImageRegion();
  RealType EvaluateMetricOverImageRegion( RealType );
  typename ImageType::Pointer EvaluateImageAtPyramidLevel( unsigned int );
  typename DeformationFieldType::Pointer EvaluateFieldFromControlPointsAtTimePoint
    ( typename ControlPointLatticeType::Pointer, unsigned int, bool );
  void InitializeDeformationFieldWithLandmarks();
//  void InitializeDeformationFieldWithDeformationFields();

  RealType EvaluateEnergyForLineSearch( RealType );
  void FindBracketingTriplet( RealType*, RealType*, RealType* );
  void LineMinimization( RealType*, RealType* );

  void BrentSearch( RealType, RealType, RealType, RealType*, RealType* );
  void BruteForceSearch( RealType, RealType, RealType, RealType*, RealType* );
  void GoldenSectionSearch( RealType, RealType, RealType, RealType*, RealType* );

private:

  typename ControlPointLatticeType::Pointer            m_TotalDeformationFieldControlPoints;
  typename ControlPointLatticeType::Pointer            m_CurrentDeformationFieldControlPoints;
  typename ControlPointLatticeType::Pointer            m_GradientFieldControlPoints;

  std::vector<typename ImageType::SizeType>            m_PyramidLevelImageSizes;

  unsigned int                                         m_LineSearchMaximumIterations;

  ArrayType                                            m_MeshResolution;
  unsigned int                                         m_TemporalMeshResolution;
  ArrayType                                            m_ImageShrinkFactors;
  vnl_vector<unsigned int>                             m_MaximumNumberOfIterations;
  MetricRadiusType                                     m_MetricRadius;
  unsigned int                                         m_SplineOrder;

  unsigned int                                         m_TemporalSplineOrder;
  bool                                                 m_WrapTime;
  RealType                                             m_TemporalOrigin;
  RealType                                             m_TemporalEnd;
  std::vector<RealType>                                m_TimePoints;

  unsigned int                                         m_CurrentLevel;
  RealType                                             m_GradientScalingFactor;
  std::vector<std::string>                             m_ImageFileNames;

  typename ImageInterpolatorType::Pointer              m_ImageInterpolator;

  typename PDEDeformableMetricType::Pointer            m_PDEDeformableMetric;
  bool                                                 m_MaximizeMetric;

};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPerfusionRegistrationFilter.txx"
#endif

#endif


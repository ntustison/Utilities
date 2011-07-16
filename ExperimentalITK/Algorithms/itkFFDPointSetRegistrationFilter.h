/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFFDPointSetRegistrationFilter.h,v $
  Language:  C++
  Date:      $Date: 2009/01/14 02:56:30 $
  Version:   $Revision: 1.3 $

  Copyright ( c ) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef _itkFFDPointSetRegistrationFilter_h_
#define _itkFFDPointSetRegistrationFilter_h_

#include "itkPointSetToPointSetFilter.h"

#include "itkBSplineControlPointImageFilter.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkImage.h"
#include "itkMacro.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkPointSet.h"
#include "itkVector.h"

#include "vnl/vnl_vector.h"

#include <map>

namespace itk {
/** \class FFDPointSetRegistrationFilter
    \brief FFD Image registration filter.
     \par
*/
template<class TFixedPointSet,
          class TMovingPointSet = TFixedPointSet,
          class TWarpedPointSet = TFixedPointSet>
class ITK_EXPORT FFDPointSetRegistrationFilter
: public PointSetToPointSetFilter<TMovingPointSet, TWarpedPointSet>
{
public:
  typedef FFDPointSetRegistrationFilter                        Self;
  typedef PointSetToPointSetFilter<TMovingPointSet, TWarpedPointSet>
                                                               Superclass;
  typedef SmartPointer<Self>                                   Pointer;
  typedef SmartPointer<const Self>                             ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ) */
  itkTypeMacro( FFDPointSetRegistrationFilter, PointSetToPointSetFilter );

  /**
   * Point set typedefs
   */
  typedef TFixedPointSet                                       FixedPointSetType;
  typedef typename FixedPointSetType::PointType                FixedPointType;
  typedef typename FixedPointSetType::PixelType                FixedPixelType;
  typedef TMovingPointSet                                      MovingPointSetType;
  typedef typename MovingPointSetType::PointType               MovingPointType;
  typedef typename MovingPointSetType::PixelType               MovingPixelType;
  typedef TWarpedPointSet                                      WarpedPointSetType;
  typedef typename WarpedPointSetType::PointType               WarpedPointType;
  typedef typename WarpedPointSetType::PixelType               WarpedPixelType;

  /**
   * Dimensionality of input and output data is assumed to be the same.
   */
  itkStaticConstMacro( Dimension, unsigned int,
                       FixedPointSetType::PointDimension );

  typedef float                                                RealType;
  typedef Image<RealType,
          itkGetStaticConstMacro( Dimension )>                 RealImageType;
  typedef FixedArray<unsigned int,
          itkGetStaticConstMacro( Dimension )>                 ArrayType;
  typedef Array<unsigned int>                                  ResizableArrayType;

  typedef std::map<WarpedPixelType, RealType>                  LabelWeightsMapType;

  /** Typedefs for B-spline filter */
  typedef Vector<RealType,
          itkGetStaticConstMacro( Dimension )>                 VectorType;
  typedef Image<VectorType,
          itkGetStaticConstMacro( Dimension )>                 DeformationFieldType;
  typedef PointSet<VectorType,
          itkGetStaticConstMacro( Dimension )>                 BSplinePointSetType;
  typedef BSplineScatteredDataPointSetToImageFilter
          <BSplinePointSetType, DeformationFieldType>          BSplineFilterType;
  typedef typename BSplineFilterType::PointDataImageType       ControlPointLatticeType;
  typedef typename BSplineFilterType::WeightsContainerType     WeightsContainerType;
  typedef BSplineControlPointImageFilter
    <ControlPointLatticeType, DeformationFieldType>            BSplineControlPointFilterType;
  typedef typename BSplineControlPointFilterType::OriginType   OriginType;
  typedef typename BSplineControlPointFilterType::SpacingType  SpacingType;
  typedef typename BSplineControlPointFilterType::SizeType     SizeType;

  typedef typename Statistics
     ::MersenneTwisterRandomVariateGenerator                   RandomizerType;

  /** Main functions */
  void RunRegistration()
    { this->Update(); }

  /** Set/Get functions ( variables are explained below ). */

  itkSetMacro( SplineOrder, unsigned int );
  itkGetConstMacro( SplineOrder, unsigned int );

  itkSetMacro( Directionality, ArrayType );
  itkGetConstMacro( Directionality, ArrayType );

  itkSetMacro( InitialMeshResolution, ArrayType );
  itkGetConstMacro( InitialMeshResolution, ArrayType );

  itkSetMacro( MaximumNumberOfIterations, ResizableArrayType );
  itkGetConstMacro( MaximumNumberOfIterations, ResizableArrayType );

  itkSetMacro( NumberOfFixedSamples, unsigned long );
  itkGetConstMacro( NumberOfFixedSamples, unsigned long );

  itkSetMacro( NumberOfMovingSamples, unsigned long );
  itkGetConstMacro( NumberOfMovingSamples, unsigned long );

  void SetNumberOfLevels( unsigned int  );
  itkGetConstMacro( NumberOfLevels, unsigned int );

  itkSetMacro( LineSearchMaximumIterations, unsigned int );
  itkGetConstMacro( LineSearchMaximumIterations, unsigned int );

  itkSetMacro( WhichGradient, unsigned short );
  itkGetConstMacro( WhichGradient, unsigned short );

  itkSetMacro( EmploySteepestDescent, bool );
  itkGetConstMacro( EmploySteepestDescent, bool );
  itkBooleanMacro( EmploySteepestDescent );

  itkSetMacro( UseInputAsSamples, bool );
  itkGetConstMacro( UseInputAsSamples, bool );
  itkBooleanMacro( UseInputAsSamples );

  itkSetMacro( UseAnisotropicCovariances, bool );
  itkGetConstMacro( UseAnisotropicCovariances, bool );
  itkBooleanMacro( UseAnisotropicCovariances );

  itkSetMacro( EmployTerm2, bool );
  itkGetConstMacro( EmployTerm2, bool );
  itkBooleanMacro( EmployTerm2 );

  itkSetMacro( Prolificacy, bool );
  itkGetConstMacro( Prolificacy, bool );
  itkBooleanMacro( Prolificacy );

  itkSetMacro( Alpha, RealType );
  itkGetConstMacro( Alpha, RealType );

  itkSetMacro( RegularizationSigma, RealType );
  itkGetConstMacro( RegularizationSigma, RealType );

  itkSetMacro( KernelSigma, RealType );
  itkGetConstMacro( KernelSigma, RealType );

  itkSetMacro( CovarianceKNeighborhood, unsigned int );
  itkGetConstMacro( CovarianceKNeighborhood, unsigned int );

  itkSetMacro( EvaluationKNeighborhood, unsigned int );
  itkGetConstMacro( EvaluationKNeighborhood, unsigned int );

  itkSetStringMacro( FilePrefix );
  itkGetStringMacro( FilePrefix );

  itkSetClampMacro( AnnealingRate, RealType, 0.00001, 1.0 );
  itkGetMacro( AnnealingRate, RealType );

  itkSetClampMacro( ExpansionFactor, RealType, 0.0, 100.0 );
  itkGetMacro( ExpansionFactor, RealType );

  itkSetMacro( Origin, OriginType );
  itkGetConstMacro( Origin, OriginType );

  itkSetMacro( Spacing, SpacingType );
  itkGetConstMacro( Spacing, SpacingType );

  itkSetMacro( Size, SizeType );
  itkGetConstMacro( Size, SizeType );
  
  itkSetMacro( BindBoundary, bool );
  itkGetConstMacro( BindBoundary, bool );
  itkBooleanMacro( BindBoundary )
  
  itkSetMacro( LabelWeights, LabelWeightsMapType );
  itkGetConstMacro( LabelWeights, LabelWeightsMapType );
  
  itkSetMacro( InitialDeformationFieldControlPoints, 
    typename ControlPointLatticeType::Pointer );
  itkGetConstMacro( InitialDeformationFieldControlPoints, 
    typename ControlPointLatticeType::Pointer );    

  typename ControlPointLatticeType::Pointer
    GetTotalDeformationFieldControlPoints()
    {
    return this->m_TotalDeformationFieldControlPoints;
    }


protected :
  /** de/constructor */
  FFDPointSetRegistrationFilter();
  ~FFDPointSetRegistrationFilter();

  void PrintSelf( std::ostream& os, Indent indent ) const;
  void GenerateData();

private :
  FFDPointSetRegistrationFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  void Initialize();

  RealType EvaluateMetricAndGradient();
  RealType EvaluateMetric( RealType );

  void IterativeSolve();  // Conjugate gradient descent

  RealType EvaluateEnergyForLineSearch( RealType );
  RealType FindBracketingTriplet( RealType*, RealType*, RealType* );
  void LineMinimization( RealType*, RealType* );

  void BrentSearch( RealType, RealType, RealType, RealType, RealType*, RealType* );
  void BruteForceSearch( RealType, RealType, RealType, RealType*, RealType* );
  void GoldenSectionSearch( RealType, RealType, RealType, RealType*, RealType* );

private :

  /** Filter Outputs */
  unsigned int                                      m_NumberOfLevels;
  ArrayType                                         m_FixedImageShrinkFactors;
  ArrayType                                         m_MovingImageShrinkFactors;
  unsigned int                                      m_LineSearchMaximumIterations;
  unsigned int                                      m_SplineOrder;
  ArrayType                                         m_InitialMeshResolution;
  ArrayType                                         m_Directionality;
  unsigned short                                    m_WhichGradient;

  ResizableArrayType                                m_MaximumNumberOfIterations;
  unsigned int                                      m_CurrentLevel;
  int                                               m_CurrentIteration;
  bool                                              m_EmploySteepestDescent;
  bool                                              m_GenerateMeanPointSet;
  bool                                              m_UseInputAsSamples;
  bool                                              m_EmployTerm2;
  bool                                              m_UseAnisotropicCovariances;

  unsigned long                                     m_NumberOfFixedSamples;
  unsigned long                                     m_NumberOfMovingSamples;
  RealType                                          m_RegularizationSigma;
  RealType                                          m_KernelSigma;
  RealType                                          m_AnnealingRate;
  RealType                                          m_CurrentAnnealing;
  unsigned int                                      m_CovarianceKNeighborhood;
  unsigned int                                      m_EvaluationKNeighborhood;
  typename RandomizerType::Pointer                  m_Randomizer;

  typename FixedPointSetType::Pointer               m_FixedSamplePoints;
  typename MovingPointSetType::Pointer              m_MovingSamplePoints;
  typename WarpedPointSetType::Pointer              m_WarpedSamplePoints;

  FixedPointType                                    m_MaxBoundary;
  FixedPointType                                    m_MinBoundary;
  bool                                              m_BindBoundary;

  typename ControlPointLatticeType::Pointer         m_InitialDeformationFieldControlPoints;
  typename ControlPointLatticeType::Pointer         m_TotalDeformationFieldControlPoints;
  typename ControlPointLatticeType::Pointer         m_CurrentDeformationFieldControlPoints;
  typename ControlPointLatticeType::Pointer         m_GradientFieldControlPoints;

  RealType                                          m_Alpha;
  RealType                                          m_ExpansionFactor;
  std::string                                       m_FilePrefix;
  bool                                              m_Prolificacy;
  
  LabelWeightsMapType                               m_LabelWeights;

  OriginType m_Origin;
  SpacingType m_Spacing;
  SizeType m_Size;

};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFFDPointSetRegistrationFilter.hxx"
#endif

#endif

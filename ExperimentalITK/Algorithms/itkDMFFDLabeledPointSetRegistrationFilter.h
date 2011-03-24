/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDMFFDLabeledPointSetRegistrationFilter.h,v $
  Language:  C++
  Date:      $Date: 2009/06/02 17:32:46 $
  Version:   $Revision: 1.4 $

  Copyright ( c ) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDMFFDLabeledPointSetRegistrationFilter_h
#define __itkDMFFDLabeledPointSetRegistrationFilter_h

#include "itkPointSetToPointSetFilter.h"

#include "itkBSplineControlPointImageFilter.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkImage.h"
#include "itkPointSet.h"
#include "itkVector.h"

#include <map>

namespace itk {
/** \class DMFFDLabeledPointSetRegistrationFilter
    \brief Point-set registration filter.
    \par Information-theoretic labeled point-set registration filter which
    uses a directly manipulated free-form deformation model.  Input consists
    of two labeled point-sets (fixed and moving).  Output is the warped
    point-set.  The total deformation field, described by control points,
    can subsequently be used to generate the resulting sampled deformation
    field.
*/
template<class TFixedPointSet, class TMovingPointSet = TFixedPointSet,
  class TWarpedPointSet = TFixedPointSet>
class ITK_EXPORT DMFFDLabeledPointSetRegistrationFilter
: public PointSetToPointSetFilter<TMovingPointSet, TWarpedPointSet>
{
public:
  typedef DMFFDLabeledPointSetRegistrationFilter         Self;
  typedef PointSetToPointSetFilter
    <TMovingPointSet, TWarpedPointSet>                   Superclass;
  typedef SmartPointer<Self>                             Pointer;
  typedef SmartPointer<const Self>                       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ) */
  itkTypeMacro( DMFFDLabeledPointSetRegistrationFilter,
    PointSetToPointSetFilter );

  /**
   * Point set typedefs
   */
  typedef TFixedPointSet                            FixedPointSetType;
  typedef typename FixedPointSetType::PointType     FixedPointType;
  typedef typename FixedPointSetType::PixelType     FixedPixelType;
  typedef TMovingPointSet                           MovingPointSetType;
  typedef typename MovingPointSetType::PointType    MovingPointType;
  typedef typename MovingPointSetType::PixelType    MovingPixelType;
  typedef TWarpedPointSet                           WarpedPointSetType;
  typedef typename WarpedPointSetType::PointType    WarpedPointType;
  typedef typename WarpedPointSetType::PixelType    WarpedPixelType;

  /**
   * Dimensionality of input and output data is assumed to be the same.
   */
  itkStaticConstMacro( Dimension, unsigned int,
    FixedPointSetType::PointDimension );

  typedef float                                     RealType;
  typedef Image<RealType,
    itkGetStaticConstMacro( Dimension )>            RealImageType;
  typedef FixedArray<unsigned int,
    itkGetStaticConstMacro( Dimension )>            ArrayType;
  typedef Array<unsigned int>                       ResizableUIntArrayType;
  typedef Array<RealType>                           ResizableRealArrayType;
  typedef std::map<WarpedPixelType, RealType>       LabelWeightsMapType;

  /** Typedefs for B-spline filter */
  typedef Vector<RealType,
    itkGetStaticConstMacro( Dimension )>            VectorType;
  typedef Image<VectorType,
    itkGetStaticConstMacro( Dimension )>            DeformationFieldType;
  typedef PointSet<VectorType,
    itkGetStaticConstMacro( Dimension )>            BSplinePointSetType;
  typedef BSplineScatteredDataPointSetToImageFilter
    <BSplinePointSetType, DeformationFieldType>     BSplineFilterType;
  typedef typename
    BSplineFilterType::PointDataImageType           ControlPointLatticeType;
  typedef typename ControlPointLatticeType::Pointer ControlPointLatticePointer;
  typedef typename
    BSplineFilterType::WeightsContainerType         WeightsContainerType;
  typedef BSplineControlPointImageFilter
    <ControlPointLatticeType, DeformationFieldType> ControlPointFilterType;
  typedef typename
    ControlPointFilterType::OriginType              OriginType;
  typedef typename
    ControlPointFilterType::SpacingType             SpacingType;
  typedef typename ControlPointFilterType::SizeType SizeType;
  typedef typename
    ControlPointFilterType::DirectionType           DirectionType;

  /** Main functions */
  void RunRegistration()
    { this->Update(); }

  /** Set/Get functions */

  /**
   * Variables governing the B-spline transformation model
   */

  /**
   * B-spline order
   */
  itkSetMacro( SplineOrder, unsigned int );
  itkGetConstMacro( SplineOrder, unsigned int );

  /**
   * B-spline mesh size at the lowest resolution.
   */
  itkSetMacro( InitialMeshResolution, ArrayType );
  itkGetConstMacro( InitialMeshResolution, ArrayType );

  /**
   * special ivar used for taggine registration.
   */
  itkSetMacro( Directionality, ArrayType );
  itkGetConstMacro( Directionality, ArrayType );

  /**
   * turn on/off the initial similarity transform calculation
   */
  itkSetMacro( CalculateInitialSimilarityTransform, bool );
  itkGetConstMacro( CalculateInitialSimilarityTransform, bool );
  itkBooleanMacro( CalculateInitialSimilarityTransform );

  /**
   * Weights for each label to emphasize different labels in the gradient.
   */
  void SetLabelWeights( LabelWeightsMapType labelWeights )
    {
    this->m_LabelWeights = labelWeights;
    this->Modified();
    }

  /**
   * The origin, spacing, and size define the transformation domain.
   */
  itkSetMacro( Origin, OriginType );
  itkGetConstMacro( Origin, OriginType );

  itkSetMacro( Spacing, SpacingType );
  itkGetConstMacro( Spacing, SpacingType );

  itkSetMacro( Size, SizeType );
  itkGetConstMacro( Size, SizeType );

  itkSetMacro( Direction, DirectionType );
  itkGetConstMacro( Direction, DirectionType );

  /**
   * Variables governing the gradient descent iterative optimization.
   */

  /**
   * Maximum number of iterations for each resolution level.
   */
  itkSetMacro( MaximumNumberOfIterations, ResizableUIntArrayType );
  itkGetConstMacro( MaximumNumberOfIterations, ResizableUIntArrayType );

  /**
   * Gradient scaling factor.
   */
  itkSetMacro( GradientScalingFactor, ResizableRealArrayType );
  itkGetConstMacro( GradientScalingFactor, ResizableRealArrayType );

  /**
   * Maximum number of iterations for each line iteration.
   */
  itkSetMacro( LineSearchMaximumIterations, unsigned int );
  itkGetConstMacro( LineSearchMaximumIterations, unsigned int );

  /**
   * Maximum step size for each line iteration.
   */
  itkSetMacro( LineSearchMaximumStepSize, RealType );
  itkGetConstMacro( LineSearchMaximumStepSize, RealType );

  /**
   * Variables governing the Jensen-Havrda-Charvat-Tsallis point-set
   * similarity measure.
   */

  /**
   * To calculate the gradient, one can either set the number of fixed
   * and moving samples to be generated randomly from the input point-sets or
   * one can use the input point-sets as samples.
   */
  itkSetMacro( NumberOfFixedSamples, unsigned long );
  itkGetConstMacro( NumberOfFixedSamples, unsigned long );

  itkSetMacro( NumberOfMovingSamples, unsigned long );
  itkGetConstMacro( NumberOfMovingSamples, unsigned long );

  itkSetMacro( UseInputAsSamples, bool );
  itkGetConstMacro( UseInputAsSamples, bool );
  itkBooleanMacro( UseInputAsSamples );

  /**
   * Anisotropic covariances are calculated from the local point-set
   * structure (using m_KernelSigma and m_Covariance KNeighborhood as a
   * weighting parameter) and the covariance generated from
   * m_PointSetSigma. and the current annealing value  Otherwise,
   * isotropic covariances are used employing m_PointSetSigma and the
   * current annealing value.
   */
  itkSetMacro( UseAnisotropicCovariances, bool );
  itkGetConstMacro( UseAnisotropicCovariances, bool );
  itkBooleanMacro( UseAnisotropicCovariances );

  itkSetMacro( KernelSigma, RealType );
  itkGetConstMacro( KernelSigma, RealType );

  itkSetMacro( CovarianceKNeighborhood, unsigned int );
  itkGetConstMacro( CovarianceKNeighborhood, unsigned int );

  itkSetMacro( PointSetSigmas, ResizableRealArrayType );
  itkGetConstMacro( PointSetSigmas, ResizableRealArrayType );

  /**
   * Describes the deterministic annealing schedule.
   */
  itkSetClampMacro( AnnealingRate, RealType, 0.00001, 1.0 );
  itkGetMacro( AnnealingRate, RealType );

  /**
   * Use the second term of the JHCT divergence.
   */
  itkSetMacro( UseRegularizationTerm, bool );
  itkGetConstMacro( UseRegularizationTerm, bool );
  itkBooleanMacro( UseRegularizationTerm );

  /**
   * The alpha parameter of the JHCT divergence controls the robustness
   * sensitivity trade-off with alpha = 1 equivalent to the ML solution and
   * alpha = 2 the robust L2 solution.
   */
  itkSetMacro( Alpha, RealType );
  itkGetConstMacro( Alpha, RealType );

  /**
   * Specifies the number of nearest neighbors (k-d tree) to be used in
   * calculating the metric and gradient values.
   */
  itkSetMacro( EvaluationKNeighborhood, unsigned int );
  itkGetConstMacro( EvaluationKNeighborhood, unsigned int );

  /**
   * Print to screen the course of optimization.
   */
  itkSetMacro( Verbose, bool );
  itkGetConstMacro( Verbose, bool );
  itkBooleanMacro( Verbose );

  /**
   * Allow one to specify and initial deformation field defined by a
   * control point lattic.e
   */
  itkSetObjectMacro( InitialDeformationFieldControlPoints,
    ControlPointLatticeType );
  itkGetObjectMacro( InitialDeformationFieldControlPoints,
    ControlPointLatticeType );

  /**
   * Get the resulting transformation model defined by a grid of control point
   * grids.
   */
  itkGetObjectMacro( TotalDeformationFieldControlPoints,
    ControlPointLatticeType );

  /**
   * Get the resulting similarity transformation model parameters.
   */
  itkGetConstMacro( FixedPointSetCentroid, FixedPointType );
  itkGetConstMacro( MovingPointSetCentroid, MovingPointType );
  itkGetConstMacro( ScaleFactor, RealType );

protected :
  /** de/constructor */
  DMFFDLabeledPointSetRegistrationFilter();
  ~DMFFDLabeledPointSetRegistrationFilter();

  void PrintSelf( std::ostream& os, Indent indent ) const;
  void GenerateData();

private :
  //purposely not implemented
  DMFFDLabeledPointSetRegistrationFilter( const Self& );
  void operator=( const Self& ); //purposely not implemented

  void Initialize();
  void CalculateInitialCentroidsAndNorms();

  RealType EvaluateMetricAndGradient();
  RealType EvaluateMetric( RealType );

  /**
   * Conjugate gradient descent and ancillary functions
   */
  void IterativeSolve();
  RealType EvaluateEnergyForLineSearch( RealType );
  RealType FindBracketingTriplet( RealType*, RealType*, RealType* );
  void LineMinimization( RealType*, RealType* );
  void BrentSearch( RealType, RealType, RealType, RealType, RealType*,
    RealType* );

private :

  /**
   * Variables governing the B-spline transformation model
   */
  unsigned int                          m_SplineOrder;
  ArrayType                             m_InitialMeshResolution;
  ArrayType                             m_Directionality;
  LabelWeightsMapType                   m_LabelWeights;
  OriginType                            m_Origin;
  SpacingType                           m_Spacing;
  SizeType                              m_Size;
  DirectionType                         m_Direction;

  FixedPointType                        m_FixedPointSetCentroid;
  MovingPointType                       m_MovingPointSetCentroid;
  RealType                              m_ScaleFactor;
  typename MovingPointSetType::Pointer  m_SimilarityTransformMovingPointSet;
  bool                                  m_CalculateInitialSimilarityTransform;

  /**
   * Variables governing the conjugate gradient descent iterative optimization
   */
  ResizableUIntArrayType                m_MaximumNumberOfIterations;
  ResizableRealArrayType                m_GradientScalingFactor;
  unsigned int                          m_LineSearchMaximumIterations;
  RealType                              m_LineSearchMaximumStepSize;
  unsigned int                          m_CurrentLevel;
  int                                   m_CurrentIteration;

  /**
   * Variables governing the Jensen-Havrda-Charvat-Tsallis point-set
   * similarity measure
   */
  bool                                  m_UseInputAsSamples;
  unsigned long                         m_NumberOfFixedSamples;
  unsigned long                         m_NumberOfMovingSamples;
  bool                                  m_UseAnisotropicCovariances;
  RealType                              m_KernelSigma;
  unsigned int                          m_CovarianceKNeighborhood;
  RealType                              m_CurrentPointSetSigma;
  ResizableRealArrayType                m_PointSetSigmas;
  RealType                              m_AnnealingRate;
  RealType                              m_CurrentAnnealing;
  bool                                  m_UseRegularizationTerm;
  RealType                              m_Alpha;
  unsigned int                          m_EvaluationKNeighborhood;

  /**
   * Control point lattices
   */
  ControlPointLatticePointer            m_InitialDeformationFieldControlPoints;
  ControlPointLatticePointer            m_TotalDeformationFieldControlPoints;
  ControlPointLatticePointer            m_CurrentDeformationFieldControlPoints;
  ControlPointLatticePointer            m_GradientFieldControlPoints;

  /**
   * Other variables.
   */
  bool                                  m_Verbose;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDMFFDLabeledPointSetRegistrationFilter.txx"
#endif

#endif

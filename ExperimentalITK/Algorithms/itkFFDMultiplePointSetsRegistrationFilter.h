/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFFDMultiplePointSetsRegistrationFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:13:41 $
  Version:   $Revision: 1.1.1.1 $

  Copyright ( c ) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef _itkFFDMultiplePointSetsRegistrationFilter_h_
#define _itkFFDMultiplePointSetsRegistrationFilter_h_

#include "itkPointSetToPointSetFilter.h"

#include "itkBSplineControlPointImageFilter.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkImage.h"
#include "itkMacro.h"
#include "itkManifoldParzenWindowsPointSetFunction.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkPointSet.h"
#include "itkVector.h"

#include "vnl/vnl_vector.h"

namespace itk {
/** \class FFDMultiplePointSetsRegistrationFilter
    \brief FFD Image registration filter.
     \par
*/
template<class TPointSet>
class ITK_EXPORT FFDMultiplePointSetsRegistrationFilter
: public PointSetToPointSetFilter<TPointSet, TPointSet>
{
public:
  typedef FFDMultiplePointSetsRegistrationFilter                        Self;
  typedef PointSetToPointSetFilter<TPointSet, TPointSet>       Superclass;
  typedef SmartPointer<Self>                                   Pointer;
  typedef SmartPointer<const Self>                             ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ) */
  itkTypeMacro( FFDMultiplePointSetsRegistrationFilter, PointSetToPointSetFilter );

  /**
   * Point set typedefs
   */
  typedef TPointSet                                            PointSetType;
  typedef typename PointSetType::PointType                     PointType;
  typedef typename PointSetType::PixelType                     PixelType;

  /**
   * Dimensionality of input and output data is assumed to be the same.
   */
  itkStaticConstMacro( Dimension, unsigned int,
                       PointSetType::PointDimension );

  typedef float                                                RealType;
  typedef Image<RealType,
          itkGetStaticConstMacro( Dimension )>                 RealImageType;
  typedef FixedArray<unsigned int,
          itkGetStaticConstMacro( Dimension )>                 ArrayType;
  typedef Array<unsigned int>                                  ResizableArrayType;

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

  typedef ManifoldParzenWindowsPointSetFunction
    <PointSetType>                                             DensityFunctionType;
  typedef typename DensityFunctionType::GaussianType           GaussianType;
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

  itkSetMacro( NumberOfSamples, unsigned long );
  itkGetConstMacro( NumberOfSamples, unsigned long );

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

  itkSetMacro( UseAnisotropicCovariances, bool );
  itkGetConstMacro( UseAnisotropicCovariances, bool );
  itkBooleanMacro( UseAnisotropicCovariances );

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

  typename ControlPointLatticeType::Pointer
    GetTotalDeformationFieldControlPoints( unsigned int i )
    {
    if ( i < this->m_TotalDeformationFieldControlPoints.size() )
      {
      return this->m_TotalDeformationFieldControlPoints[i];
      }
    else
      {
      return NULL;
      }
    }


protected :
  /** de/constructor */
  FFDMultiplePointSetsRegistrationFilter();
  ~FFDMultiplePointSetsRegistrationFilter();

  void PrintSelf( std::ostream& os, Indent indent ) const;
  void GenerateData();

private :
  FFDMultiplePointSetsRegistrationFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  void Initialize();

  void EvaluateGradient();
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
  unsigned int                                      m_LineSearchMaximumIterations;
  unsigned int                                      m_SplineOrder;
  ArrayType                                         m_InitialMeshResolution;
  ArrayType                                         m_Directionality;
  unsigned short                                    m_WhichGradient;

  ResizableArrayType                                m_MaximumNumberOfIterations;
  unsigned int                                      m_CurrentLevel;
  int                                               m_CurrentIteration;
  bool                                              m_EmploySteepestDescent;
  bool                                              m_UseInputAsSamples;
  bool                                              m_EmployTerm2;
  bool                                              m_UseAnisotropicCovariances;

  unsigned long                                     m_NumberOfSamples;
  RealType                                          m_RegularizationSigma;
  RealType                                          m_KernelSigma;
  RealType                                          m_AnnealingRate;
  RealType                                          m_CurrentAnnealing;
  unsigned int                                      m_CovarianceKNeighborhood;
  unsigned int                                      m_EvaluationKNeighborhood;
  typename RandomizerType::Pointer                  m_Randomizer;

  PointType                                         m_MaxBoundary;
  PointType                                         m_MinBoundary;

  std::vector<typename ControlPointLatticeType::Pointer>     m_TotalDeformationFieldControlPoints;
  std::vector<typename ControlPointLatticeType::Pointer>     m_CurrentDeformationFieldControlPoints;
  std::vector<typename ControlPointLatticeType::Pointer>     m_GradientFieldControlPoints;

  RealType                                          m_Alpha;
  RealType                                          m_ExpansionFactor;
  std::string                                       m_FilePrefix;
  bool                                              m_Prolificacy;

  OriginType m_Origin;
  SpacingType m_Spacing;
  SizeType m_Size;



};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFFDMultiplePointSetsRegistrationFilter.txx"
#endif

#endif

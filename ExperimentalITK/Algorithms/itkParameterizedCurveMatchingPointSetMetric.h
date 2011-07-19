/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkParameterizedCurveMatchingPointSetMetric_h
#define __itkParameterizedCurveMatchingPointSetMetric_h

#include "itkPointSetToPointSetMetric.h"

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkImage.h"

namespace itk {

/** \class ParameterizedCurveMatchingPointSetMetric
 *
 * \brief Implementation of the Jensen Havrda Charvat Tsallis Point Set metric.
 *
 *
 * \par REFERENCE
 *
 * \ingroup ITK-Registration
 */

template<class TPointSet>
class ITK_EXPORT ParameterizedCurveMatchingPointSetMetric :
    public PointSetToPointSetMetric<TPointSet, TPointSet>
{
public:
  /** Standard class typedefs. */
  typedef ParameterizedCurveMatchingPointSetMetric             Self;
  typedef PointSetToPointSetMetric<TPointSet, TPointSet>       Superclass;
  typedef SmartPointer<Self>                                   Pointer;
  typedef SmartPointer<const Self>                             ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods) */
  itkTypeMacro( ParameterizedCurveMatchingPointSetMetric,
    PointSetToPointSetMetric );

  typedef TPointSet                                 PointSetType;
  typedef typename PointSetType::PointsContainer    PointsContainer;
  typedef typename PointsContainer::ConstIterator   PointsContainerConstIterator;

  itkStaticConstMacro( PointDimension, unsigned int, TPointSet::PointDimension );

  /** Types transferred from the base class */
  itkSuperclassTraitMacro( MeasureType );
  itkSuperclassTraitMacro( DerivativeType );
  itkSuperclassTraitMacro( LocalDerivativeType );
  itkSuperclassTraitMacro( PointType );
  itkSuperclassTraitMacro( CoordRepType );
  itkSuperclassTraitMacro( PointIdentifier );
  itkSuperclassTraitMacro( NeighborsIdentifierType );

  typedef MeasureType                                   RealType;

  typedef typename PointType::Vector<RealType,
    itkGetStaticConst( PointDimension )>                CurvePointType;
  typedef Image<ParametricPointType,
    itkGetStaticConst( PointDimension )>                CurveImageType;
  typedef typename CurveImageType::Pointer              CurveImagePointer;
  typedef PointSet<CurvePointType, 1>                   SampledCurveType;
  typedef BSplineScatteredDataPointSetToImageFilter
    <SampledCurveType, CurveImageType>                  BSplineFilterType;

  /**
   * Other typedefs
   */
  typedef ManifoldParzenWindowsPointSetFunction
    <PointSetType, RealType>                            DensityFunctionType;
  typedef typename DensityFunctionType::GaussianType    GaussianType;
  typedef typename DensityFunctionType::Pointer         DensityFunctionPointer;

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void Initialize( void ) throw ( ExceptionObject );

  /**
   * This method returns the value of the metric based on the current
   * transformation(s).
   */
  virtual MeasureType GetValue() const;

  /**
   * This method returns the derivative based on the current
   * transformation(s).
   */
  virtual void GetDerivative( DerivativeType & ) const;

  /**
   * This method returns the derivative and value based on the current
   * transformation(s).
   */
  virtual void GetValueAndDerivative( MeasureType &, DerivativeType & ) const;

  /**
   * Set spline order used to construct the fixed and moving B-spline curves.
   * Default = 3.
   */
  itkSetMacro( SplineOrder, unsigned int );

  /**
   * Get spline order used to construct the fixed and moving B-spline curves.
   * Default = 3.
   */
  itkGetConstMacro( SplineOrder, unsigned int );

  /**
   * Set number of fitting levels.  Default = 4.
   */
  itkSetMacro( NumberOfFittingLevels, unsigned int );

  /**
   * Set number of fitting levels.  Default = 4.
   */
  itkGetConstMacro( NumberOfFittingLevels, unsigned int );

  /**
   * Set sampling rate.  The B-spline curve is parameterized to be between
   * 0 and 1.  The number of samples defining the curve is 1 / m_SamplingRate.
   * Default = 0.0001.
   */
  itkSetMacro( SamplingRate, RealType );

  /**
   * Get sampling rate.  Default = 0.0001.
   */
  itkGetConstMacro( SamplingRate, RealType );

  /**
   * Set boolean definining whether or not the curves are closed.
   * Default = false.
   */
  itkSetMacro( UseClosedCurves, bool );

  /**
   * Get boolean macro definining whether or not the curves are closed.
   * Default = false.
   */
  itkGetConstMacro( UsedClosedCurves, bool );

  /**
   * Set/get boolean macro definining whether or not the curves are closed.
   * Default = false.
   */
  itkBooleanMacro( UsedClosedCurves );

  /**
   * Set parametric window.  The B-spline curve is parameterized to be between
   * 0 and 1.  The parametric window defines the search space for the dynamic
   * programming optimization.
   * Default = 0.1.
   */
  itkSetMacro( ParametricWindow, RealType );

  /**
   * Get parametric window.  Default = 0.1.
   */
  itkGetConstMacro( ParametricWindow, RealType );

  /** Pure virtual function from the parent class not needed in this class */
  virtual MeasureType GetLocalNeighborhoodValue( const PointType ) const
  {
    itkExceptionMacro( "This function should not be accessed." );
    return 0;
  }

  /** Pure virtual function from the parent class not needed in this class */
  virtual void GetLocalNeighborhoodValueAndDerivative( const PointType,
    MeasureType &, LocalDerivativeType & ) const
  {
    itkExceptionMacro( "This function should not be accessed." );
  }

protected:
  ParameterizedCurveMatchingPointSetMetric();
  ~ParameterizedCurveMatchingPointSetMetric() {}

  void PrintSelf( std::ostream& os, Indent indent ) const;

private:
  //purposely not implemented
  ParameterizedCurveMatchingPointSetMetric( const Self& );
  void operator=( const Self& );

  CurveImagePointer GenerateBSplineCurveControlPoints( PointSetType * );

  CurveImagePointer ReconstructBSplineCurve( CurveImageType * );

  unsigned int                   m_SplineOrder;
  unsigned int                   m_NumberOfFittingLevels;
  RealType                       m_SamplingRate;
  bool                           m_UseClosedCurves;
  RealType                       m_ParametricWindow;

  CurveImagePointer              m_FixedBSplineCurveControlPointImage;
  CurveImagePointer              m_MovingBSplineCurveControlPointImage;
};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkParameterizedCurveMatchingPointSetMetric.hxx"
#endif

#endif

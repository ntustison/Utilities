/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkEuclideanDistancePointSetMetric.h,v $
  Language:  C++
  Date:      $Date: 2008/12/06 04:52:56 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkEuclideanDistancePointSetMetric_h
#define __itkEuclideanDistancePointSetMetric_h

#include "itkPointSetToPointSetMetric.h"

#include "itkIdentityTransform.h"
#include "itkKdTreeGenerator.h"
#include "itkListSample.h"
#include "itkVector.h"
#include "itkWeightedCentroidKdTreeGenerator.h"

namespace itk {

/** \class EuclideanDistancePointSetMetric
 *
 *
 *
 *
 */
template<class TPointSet>
class ITK_EXPORT EuclideanDistancePointSetMetric :
    public PointSetToPointSetMetric<TPointSet, TPointSet>
{
public:
  /** Standard class typedefs. */
  typedef EuclideanDistancePointSetMetric       Self;
  typedef PointSetToPointSetMetric<TPointSet, TPointSet> Superclass;
  typedef SmartPointer<Self>                             Pointer;
  typedef SmartPointer<const Self>                       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods) */
  itkTypeMacro( EuclideanDistancePointSetMetric,
    PointSetToPointSetMetric );

  itkStaticConstMacro( PointDimension, unsigned int,
                       TPointSet::PointDimension );

  /** Types transferred from the base class */
  typedef typename Superclass::TransformType            TransformType;
  typedef typename Superclass::TransformPointer         TransformPointer;
  typedef typename Superclass::TransformParametersType  TransformParametersType;
  typedef typename Superclass::TransformJacobianType    TransformJacobianType;

  typedef typename Superclass::MeasureType              MeasureType;
  typedef typename Superclass::DerivativeType           DerivativeType;

  typedef TPointSet                                     PointSetType;
  typedef typename PointSetType::Pointer                PointSetPointer;
  typedef typename PointSetType::PointType              PointType;



  /**
   * Other typedefs
   */
  typedef double                                           RealType;
  typedef IdentityTransform<RealType, PointDimension>      DefaultTransformType;
  typedef Vector<typename PointSetType::CoordRepType,
    itkGetStaticConstMacro( PointDimension )>              MeasurementVectorType;
  typedef typename Statistics::ListSample
    <MeasurementVectorType>                                SampleType;
  typedef typename Statistics
    ::WeightedCentroidKdTreeGenerator<SampleType>          TreeGeneratorType;
  typedef typename TreeGeneratorType::KdTreeType           KdTreeType;
  typedef typename KdTreeType
    ::InstanceIdentifierVectorType                         NeighborhoodIdentifierType;


  /**
   * Public function definitions
   */

  itkSetConstObjectMacro( FixedPointSet, PointSetType );
  itkGetConstObjectMacro( FixedPointSet, PointSetType );

  itkSetConstObjectMacro( MovingPointSet, PointSetType );
  itkGetConstObjectMacro( MovingPointSet, PointSetType );

  itkSetObjectMacro( Transform, TransformType );
  itkGetObjectMacro( Transform, TransformType );

  /**
   * The only transform type used is Identity
   */
  unsigned int GetNumberOfParameters( void ) const
    { return m_Transform->GetNumberOfParameters(); }

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void Initialize( void ) throw ( ExceptionObject );

  /** Get the number of values */
  unsigned int GetNumberOfValues() const;

  /** Get the derivatives of the match measure. */
  void GetDerivative( const TransformParametersType & parameters,
                      DerivativeType & Derivative ) const;

  /**  Get the value for single valued optimizers. */
  MeasureType GetValue( const TransformParametersType & parameters ) const;

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative( const TransformParametersType & parameters,
    MeasureType& Value, DerivativeType& Derivative ) const;

  itkSetMacro( UseWithRespectToTheMovingPointSet, bool );
  itkGetConstMacro( UseWithRespectToTheMovingPointSet, bool );
  itkBooleanMacro( UseWithRespectToTheMovingPointSet );

protected:
  EuclideanDistancePointSetMetric();
  ~EuclideanDistancePointSetMetric() {}

  void PrintSelf( std::ostream& os, Indent indent ) const;

private:
  //purposely not implemented
  EuclideanDistancePointSetMetric(const Self&);
  void operator=(const Self&);

  bool                                     m_UseWithRespectToTheMovingPointSet;

  typename TreeGeneratorType::Pointer      m_KdTreeGenerator;
  typename SampleType::Pointer             m_SamplePoints;

  TransformPointer                         m_Transform;
};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkEuclideanDistancePointSetMetric.hxx"
#endif

#endif

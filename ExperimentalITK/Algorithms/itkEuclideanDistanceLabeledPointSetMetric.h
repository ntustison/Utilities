/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkEuclideanDistanceLabeledPointSetMetric.h,v $
  Language:  C++
  Date:      $Date: 2008/12/06 04:52:56 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkEuclideanDistanceLabeledPointSetMetric_h
#define _itkEuclideanDistanceLabeledPointSetMetric_h

#include "itkPointSetToPointSetMetric.h"

#include "itkIdentityTransform.h"

namespace itk {

/** \class EuclideanDistanceLabeledPointSetMetric
 *
 *
 * This class is templated over the fixed point-set type, moving
 * point-set type.
 *
 */
template<class TPointSet>
class ITK_EXPORT EuclideanDistanceLabeledPointSetMetric :
    public PointSetToPointSetMetric<TPointSet, TPointSet>
{
public:
  /** Standard class typedefs. */
  typedef EuclideanDistanceLabeledPointSetMetric     Self;
  typedef PointSetToPointSetMetric<TPointSet, TPointSet>      Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods) */
  itkTypeMacro( EuclideanDistanceLabeledPointSetMetric,
    PointSetToLabelPointSetMetric );

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
  typedef typename PointSetType::PointType              PointType;
  typedef typename PointSetType::PixelType              PixelType;

  /**
   * Other typedefs
   */
  typedef double                                        RealType;
  typedef IdentityTransform<RealType, PointDimension>   DefaultTransformType;

  typedef std::vector<PixelType>                        LabelSetType;

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

  void SetFixedLabelSet( LabelSetType labels )
    {
    typename LabelSetType::const_iterator iter;
    for ( iter = labels.begin(); iter != labels.end(); ++iter )
      {
      this->m_FixedLabelSet.push_back( *iter );
      }
    this->Modified();
    }
  LabelSetType* GetFixedLabelSet()
    {
    return &this->m_FixedLabelSet;
    }
  void SetMovingLabelSet( LabelSetType labels )
    {
    typename LabelSetType::const_iterator iter;
    for ( iter = labels.begin(); iter != labels.end(); ++iter )
      {
      this->m_MovingLabelSet.push_back( *iter );
      }
    this->Modified();
    }
  LabelSetType* GetMovingLabelSet()
    {
    return &this->m_MovingLabelSet;
    }


protected:
  EuclideanDistanceLabeledPointSetMetric();
  ~EuclideanDistanceLabeledPointSetMetric() {}

  void PrintSelf( std::ostream& os, Indent indent ) const;

private:
  EuclideanDistanceLabeledPointSetMetric(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  bool                                     m_UseWithRespectToTheMovingPointSet;

  TransformPointer                         m_Transform;

  LabelSetType                             m_FixedLabelSet;
  LabelSetType                             m_MovingLabelSet;

};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkEuclideanDistanceLabeledPointSetMetric.hxx"
#endif

#endif


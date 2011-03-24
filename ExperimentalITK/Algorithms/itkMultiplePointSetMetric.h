/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMultiplePointSetMetric.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:13:43 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMultiplePointSetMetric_h
#define __itkMultiplePointSetMetric_h

#include "itkProcessObject.h"

#include "itkArray.h"
#include "itkArray2D.h"

namespace itk
{

/** \class MultiplePointSetMetric
 * \brief Computes similarity between multiple point sets.
 *
 * This Class is templated over the type of the point-sets.  This particular
 * class is the base class for a hierarchy of multiple point-set metrics.
 *
 * This class computes a value that measures the similarity between the
 * multiple point-sets.
 *
 * \ingroup RegistrationMetrics
 *
 */

template <class TPointSet>
class ITK_EXPORT MultiplePointSetMetric : public ProcessObject
{
public:

  /** Standard class typedefs. */
  typedef MultiplePointSetMetric          Self;
  typedef ProcessObject                   Superclass;
  typedef SmartPointer<Self>              Pointer;
  typedef SmartPointer<const Self>        ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( MultiplePointSetMetric, ProcessObject );

  /**  Type of the pointset. */
  typedef TPointSet                                     PointSetType;
  typedef typename PointSetType::ConstPointer           PointSetConstPointer;

  /** Constants for the pointset dimensions */
  itkStaticConstMacro( PointSetDimension, unsigned int,
                       TPointSet::PointDimension );

  typedef typename PointSetType
    ::PointsContainer::ConstIterator                    PointIterator;
  typedef typename PointSetType
    ::PointDataContainer::ConstIterator                 PointDataIterator;

  typedef double                                        RealType;

  /**  Type of the measure. */
  typedef Array<RealType>                               MeasureType;

  /**  Type of the derivative. */
  typedef Array2D<RealType>                             DerivativeType;

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void Initialize( void ) throw ( ExceptionObject );

  /** This method returns the value of the cost function corresponding
    * to the specified parameters.
    * This method MUST be overloaded by derived classes   */
  virtual MeasureType GetValue() const = 0;

  /** Return the number of values that are computed by the
   *  multivalued cost function.
   *  This method MUST be overloaded by derived classes */
  virtual unsigned int GetNumberOfValues() const = 0;

  /** This method returns the derivative of the cost function corresponding
    * to the specified parameters
    * This method MUST be overloaded by derived classes   */
  virtual void GetDerivative( DerivativeType & derivative ) const = 0;

  /** Set/Get the image input of this process object.  */
  virtual void SetInput( const TPointSet *pointSet );
  virtual void SetInput( unsigned int, const TPointSet *pointSet );
  const PointSetType* GetInput( void );
  const PointSetType* GetInput( unsigned int idx );


protected:
  MultiplePointSetMetric();
  virtual ~MultiplePointSetMetric() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  MultiplePointSetMetric(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultiplePointSetMetric.txx"
#endif

#endif


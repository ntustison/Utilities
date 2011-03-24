/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkListSampleDAgostinoPearsonNormalityTest.h,v $
  Language:  C++
  Date:      $$
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkListSampleDAgostinoPearsonNormalityTest_h
#define __itkListSampleDAgostinoPearsonNormalityTest_h

#include "itkProcessObject.h"

namespace itk {
namespace Statistics {

/** \class ListSampleDAgostinoPearsonNormalityTest
 * \brief Base class for filters that take a list sample as an input and output
 * another list sample.
 *
 * ListSampleDAgostinoPearsonNormalityTest is the base class for all process objects that output
 * list sample data, and require list sample data as input. Specifically, this class
 * defines the SetInput() method for defining the input to a filter.
 *
 * \ingroup ListSampleFilters
 *
 */
template <class TListSample>
class ITK_EXPORT ListSampleDAgostinoPearsonNormalityTest : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ListSampleDAgostinoPearsonNormalityTest Self;
  typedef ProcessObject Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( ListSampleDAgostinoPearsonNormalityTest, ProcessObject );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Some convenient typedefs. */
  typedef double                                         RealType;

  typedef TListSample                                    ListSampleType;
  typedef typename ListSampleType::MeasurementVectorType MeasurementVectorType;


  /** Set the list sample input of this object.  */
  void SetInput( const ListSampleType *input );

  /** Get the list sample input of this object.  */
  ListSampleType * GetInput();

  MeasurementVectorType CalculateTestStatistic();

protected:
  ListSampleDAgostinoPearsonNormalityTest();
  ~ListSampleDAgostinoPearsonNormalityTest() {};

private:
  ListSampleDAgostinoPearsonNormalityTest( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  typename ListSampleType::ConstPointer              m_ListSample;

};

} // end namespace Statistics
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkListSampleDAgostinoPearsonNormalityTest.txx"
#endif

#endif

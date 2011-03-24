/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkWatsonDistributionListSampleFilter.h,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkWatsonDistributionListSampleFilter_h
#define __itkWatsonDistributionListSampleFilter_h

#include "itkListSampleToListSampleFilter.h"

#include "itkVariableSizeMatrix.h"

namespace itk {
namespace Statistics {

/** \class WatsonDistributionListSampleFilter
 * \brief Base class of filters intended to generate scalar samples from
 * intensity samples.
 *
 */

template<class TVectorListSample, class TScalarListSample>
class ITK_EXPORT WatsonDistributionListSampleFilter
: public ListSampleToListSampleFilter<TVectorListSample, TScalarListSample>
{
public:
  /**
   * Standard class typedefs.
   */
  typedef WatsonDistributionListSampleFilter                  Self;
  typedef ListSampleToListSampleFilter
    <TIntensityListSample, TScalarListSample>                 Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /**
   * Standard macros
   */
  itkTypeMacro( WatsonDistributionListSampleFilter,
    ListSampleToScalarListSampleFilter );

  /**
   * Method for creation through the object factory.
   */
  itkNewMacro( Self );

  /**
   * Conveneient typedefs
   */
  typedef float                                       RealType;
  typedef TVectorListSample                           VectorListSampleType;
  typedef TScalarListSample                           ScalarListSampleType;
  typedef VectorListSampleType::MeasurementVectorType VectorType;
  typedef VariableSizeMatrix<RealType>                TensorType;

protected:
  WatsonDistributionListSampleFilter();
  virtual ~WatsonDistributionListSampleFilter();

  void PrintSelf( std::ostream& os, Indent indent ) const;

  virtual void GenerateData();

private:
  WatsonDistributionListSampleFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  /**
   * Functions to calculate the various parameters of the distribution
   */
  void CalculateMeanVectorAndInertiaTensor();
  RealType CalculateKappa();
  RealType CalculateM();

  VectorType                                           m_MeanVector;
  TensorType                                           m_InertiaTensor;
  RealType                                             m_M;
  RealType                                             m_Kappa;

}; // end of class

} // end of namespace Statistics
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkWatsonDistributionListSampleFilter.txx"
#endif

#endif

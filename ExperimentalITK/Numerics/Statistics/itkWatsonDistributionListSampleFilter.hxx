/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkWatsonDistributionListSampleFilter.hxx,v $
  Language:  C++
  Date:      $Date:  $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkWatsonDistributionListSampleFilter_hxx
#define __itkWatsonDistributionListSampleFilter_hxx

#include "itkWatsonDistributionListSampleFilter.h"

namespace itk {
namespace Statistics {

template<class TIntensityListSample, class TScalarListSample>
WatsonDistributionListSampleFilter<TIntensityListSample, TScalarListSample>
::WatsonDistributionListSampleFilter()
{
  this->AllocateOutput();
  this->GetOutput()->SetMeasurementVectorSize( 1 );
}

template<class TIntensityListSample, class TScalarListSample>
WatsonDistributionListSampleFilter<TIntensityListSample, TScalarListSample>
::~WatsonDistributionListSampleFilter()
{
}

template<class TIntensityListSample, class TScalarListSample>
void
WatsonDistributionListSampleFilter<TIntensityListSample, TScalarListSample>
::GenerateData()
{
  /**
   * Calculate the governing parameters of the distribution.
   */
  TensorType T = this->CalculateMeanVectorAndInertiaTensor();
  RealType kappa = this->CalculateKappa();
  RealType M = this->CalculateM();

  RealType invM = 1.0 / M;

  const unsigned int scalarMeasurementVectorSize =
    this->GetOutput()->GetMeasurementVectorSize();

  /**
   * Convert each measurement vector of the sample to a scalar sample.
   */
  typename VectorListSampleType::ConstIterator It =
    this->GetInput()->Begin();
  while( It != this->GetInput()->End() )
    {
    typename VectorListSampleType::MeasurementVectorType vector
      = It.GetMeasurementVector();
    ScalarMeasurementVectorType measurement;
    measurement.SetSize( scalarMeasurementVectorSize );

    RealType innerProduct = 0.0;
    for( unsigned int d = 0; d < vector.Size(); d++ )
      {
      innerProduct += ( vector[d] * this->m_MeanVector[d] );
      }
    measurement[0] = invM * vcl_exp( this->m_Kappa *
      vnl_math_sqr( innerProduct ) );

    this->GetOutput()->PushBack( measurement );

    ++It;
    }
}

template<class TIntensityListSample, class TScalarListSample>
void
WatsonDistributionListSampleFilter<TIntensityListSample, TScalarListSample>
::CalculateMeanVectorAndInertiaTensor()
{
  this->m_MeanVector.SetSize( this->GetInput()->GetMeasurementVectorSize() );
  this->m_InertiaTensor.SetSize( this->GetInput()->GetMeasurementVectorSize(),
    this->GetMeasurementVectorSize() );

  this->m_MeanVector.Fill( 0.0 );
  this->m_InertiaTensor.Fill( 0.0 );

  typename VectorListSampleType::ConstIterator It =
    this->GetInput()->Begin();
  while( It != this->GetInput()->End() )
    {
    typename VectorListSampleType::MeasurementVectorType measurement =
      It.GetMeasurementVector();
    measurement.SetSize( this->GetMeasurementVectorSize() );
    for( unsigned int i = 0; i < T.Rows(); i++ )
      {
      this->m_MeanVector[i] += measurement[i];
      for( unsigned int j = 0; j < T.Cols(); j++ )
        {
        this->m_InertiaTensor(i, j) += measurement[i] * measurement[j];
        }
      }
    }

  for( unsigned int i = 0; i < this->m_MeanVector.Size(); i++ )
    {
    this->m_MeanVector[i] /= static_cast<RealType>( this->GetInput()->Size() );
    }
  this->m_InertiaTensor /= static_cast<RealType>( this->GetInput()->Size() );
}

template<class TIntensityListSample, class TScalarListSample>
void
WatsonDistributionListSampleFilter<TIntensityListSample, TScalarListSample>
::CalculateKappa()
{

}

template<class TIntensityListSample, class TScalarListSample>
void
WatsonDistributionListSampleFilter<TIntensityListSample, TScalarListSample>
::CalculateM()
{

}

template<class TIntensityListSample, class TScalarListSample>
void
WatsonDistributionListSampleFilter<TIntensityListSample, TScalarListSample>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );

  os << indent << "Parameters: " << std::endl;
  os << indent << "  M = " << this->m_M << std::endl;
  os << indent << "  Kappa = " << this->m_Kappa << std::endl;
  os << indent << "  Mean vector = " << this->m_MeanVector << std::endl;
  os << indent << "  Inertia tensor = " << this->m_InertiaTensor << std::endl;
}


} // end of namespace Statistics
} // end of namespace itk


#endif

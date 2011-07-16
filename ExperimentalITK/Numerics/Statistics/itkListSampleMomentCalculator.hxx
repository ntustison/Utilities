/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkListSampleMomentCalculator.hxx,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision: $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkListSampleMomentCalculator_hxx
#define _itkListSampleMomentCalculator_hxx

#include "itkListSampleMomentCalculator.h"

namespace itk {
namespace Statistics {

/**
 *
 */
template <class TListSample>
ListSampleMomentCalculator<TListSample>
::ListSampleMomentCalculator()
{
  // Modify superclass default values, can be overridden by subclasses
  this->SetNumberOfRequiredInputs( 1 );
}

template <class TListSample>
void
ListSampleMomentCalculator<TListSample>
::SetInput( const TListSample *input )
{
  this->m_ListSample = const_cast<ListSampleType *>( input );
}

template <class TListSample>
typename ListSampleMomentCalculator<TListSample>::ListSampleType *
ListSampleMomentCalculator<TListSample>
::GetInput()
{
  return const_cast<ListSampleType *>( this->m_ListSample.GetPointer() );
}

template <class TListSample>
typename ListSampleMomentCalculator<TListSample>
::MeasurementVectorType
ListSampleMomentCalculator<TListSample>
::GetMean()
{
  MeasurementVectorType mean;
  mean.SetSize( this->GetInput()->GetMeasurementVectorSize() );
  mean.Fill( 0.0 );

  typename ListSampleType::ConstIterator It = this->GetInput()->Begin();
  while( It != this->GetInput()->End() )
    {
    MeasurementVectorType measurement = It.GetMeasurementVector();
    for( unsigned int d = 0; d < mean.Size(); d++ )
      {
      mean[d] += measurement[d];
      }
    ++It;
    }
  for( unsigned int d = 0; d < mean.Size(); d++ )
    {
    mean[d] /= static_cast<double>( this->GetInput()->Size() );
    }

  return mean;
}

template <class TListSample>
typename ListSampleMomentCalculator<TListSample>
::MeasurementVectorType
ListSampleMomentCalculator<TListSample>
::GetStandardizedMoment( unsigned int k )
{
  MeasurementVectorType moment;
  moment.SetSize( this->GetInput()->GetMeasurementVectorSize() );
  moment.Fill( 0.0 );

  if( k == 0 )
    {
    moment.Fill( 1.0 );
    return moment;
    }
  if( k == 1 )
    {
    moment.Fill( 0.0 );
    return moment;
    }

  MeasurementVectorType mean = this->GetMean();

  MeasurementVectorType numerator;
  numerator.SetSize( this->GetInput()->GetMeasurementVectorSize() );
  numerator.Fill( 0.0 );
  MeasurementVectorType denominator;
  denominator.SetSize( this->GetInput()->GetMeasurementVectorSize() );
  denominator.Fill( 0.0 );

  typename ListSampleType::ConstIterator It = this->GetInput()->Begin();
  while( It != this->GetInput()->End() )
    {
    MeasurementVectorType measurement = It.GetMeasurementVector();
    for( unsigned int d = 0; d < mean.Size(); d++ )
      {
      numerator[d] += vcl_pow( static_cast<double>( measurement[d] - mean[d] ),
        static_cast<double>( k ) );
      denominator[d] += vcl_pow( static_cast<double>( measurement[d] - mean[d] ),
        static_cast<double>( 2.0 ) );
      }
    ++It;
    }
  for( unsigned int d = 0; d < moment.Size(); d++ )
    {
    moment[d] = ( numerator[d] / static_cast<double>( this->GetInput()->Size() ) ) /
      vcl_pow( ( denominator[d] / static_cast<double>( this->GetInput()->Size() ) ),
      0.5 * static_cast<double>( k ) );
    }
  return moment;
}

} // end of namespace Statistics
} // end of namespace itk

#endif

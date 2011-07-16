/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkListSampleDAgostinoPearsonNormalityTest.hxx,v $
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
#ifndef _itkListSampleDAgostinoPearsonNormalityTest_hxx
#define _itkListSampleDAgostinoPearsonNormalityTest_hxx

#include "itkListSampleDAgostinoPearsonNormalityTest.h"

#include "itkChiSquareDistribution.h"
#include "itkListSampleMomentCalculator.h"

namespace itk {
namespace Statistics {

/**
 *
 */
template <class TListSample>
ListSampleDAgostinoPearsonNormalityTest<TListSample>
::ListSampleDAgostinoPearsonNormalityTest()
{
  // Modify superclass default values, can be overridden by subclasses
  this->SetNumberOfRequiredInputs( 1 );
}

template <class TListSample>
void
ListSampleDAgostinoPearsonNormalityTest<TListSample>
::SetInput( const TListSample *input )
{
  this->m_ListSample = const_cast<ListSampleType *>( input );
}

template <class TListSample>
typename ListSampleDAgostinoPearsonNormalityTest<TListSample>::ListSampleType *
ListSampleDAgostinoPearsonNormalityTest<TListSample>
::GetInput()
{
  return const_cast<ListSampleType *>( this->m_ListSample.GetPointer() );
}

template <class TListSample>
typename ListSampleDAgostinoPearsonNormalityTest<TListSample>
::MeasurementVectorType
ListSampleDAgostinoPearsonNormalityTest<TListSample>
::CalculateTestStatistic()
{
  RealType n = static_cast<RealType>( this->GetInput()->Size() );

//   RealType u1_g1 = 0.0;
  RealType u2_g1 = ( 6.0 * ( n - 2.0 ) ) / ( ( n + 1.0 ) * ( n + 3.0 ) );
//   RealType m1_g1 = 0.0;
  RealType m2_g1 = ( 36.0 * ( n - 7.0 ) * ( n * n - 2.0 * n - 5.0 ) ) /
    ( ( n - 2.0 ) * ( n + 5.0 ) * ( n + 7.0 ) * ( n + 9.0 ) );

  RealType u1_g2 = -6.0 / ( n + 1.0 );
  RealType u2_g2 = ( 24.0 * n * ( n - 2.0 ) * ( n - 3.0 ) ) /
    ( vnl_math_sqr( n + 1.0 ) * ( n + 3.0 ) * ( n + 5.0 ) );
  RealType m1_g2 = ( 6.0 * ( n * n - 5.0 * n + 2.0 ) ) /
    ( ( n + 7.0 ) * ( n + 9.0 ) ) * vcl_sqrt(
    ( 6.0 * ( n + 3.0 ) * ( n + 5.0 ) ) / ( n * ( n - 2.0 ) * ( n - 3.0 ) ) );
//   RealType m2_g2 = ( 36.0 * ( 15.0 * vcl_pow( n, 6.0 ) - 36.0 *
//     vcl_pow( n, 5.0 ) - 628.0 * vcl_pow( n, 4.0 ) + 982.0 * vcl_pow( n, 3.0 ) +
//     5777.0 * n * n - 6402.0 * n + 900 ) ) / ( n * ( n - 3.0 ) * ( n - 2.0 ) *
//     ( n + 7.0 ) * ( n + 9.0 ) * ( n + 11.0 ) * ( n + 13.0 ) );

  RealType W = vcl_sqrt( vcl_sqrt( 2.0 * m2_g1 + 4.0 ) - 1.0 );
  RealType d = 1 / vcl_sqrt( vcl_log( W ) );
  RealType a = vcl_sqrt( 2.0 / ( W * W - 1.0 ) );

  RealType A = 6.0 + 8.0 / m1_g2 * ( 2.0 / m1_g2 + vcl_sqrt( 1.0 + 4.0 /
    vnl_math_sqr( m1_g2 ) ) );

  typedef ListSampleMomentCalculator<ListSampleType> CalculatorType;
  typename CalculatorType::Pointer calculator = CalculatorType::New();
  calculator->SetInput( this->GetInput() );

  typename ListSampleType::MeasurementVectorType p;
  p.SetSize( this->GetInput()->GetMeasurementVectorSize() );

  for( unsigned int n = 0; n < p.Size(); n++ )
    {
    RealType g1 = calculator->GetStandardizedMoment( 3 )[n];
    RealType g2 = calculator->GetStandardizedMoment( 4 )[n] - 3.0 ;


    RealType Z1_g1 = d * vcl_log( g1 / ( a * vcl_sqrt( u2_g1 ) ) +
      vcl_sqrt( g1 * g1 / ( a * a * u2_g1 ) + 1.0 ) );

    RealType Z2_g2 = vcl_sqrt( 4.5 * A ) * ( 1.0 - 2.0 / ( 9.0 * A ) -
      vcl_pow( ( 1.0 - 2.0 / A ) / ( 1.0 + g2 - u1_g2 / ( vcl_sqrt( u2_g2 ) *
      vcl_sqrt( 2.0 / ( A - 4.0 ) ) ) ), 1.0 / 3.0 ) );

    RealType K2 = Z1_g1 * Z1_g1 + Z2_g2 * Z2_g2;

    typename Statistics::ChiSquareDistribution::Pointer chidistribution =
      Statistics::ChiSquareDistribution::New();
    p[n] = 1.0 - chidistribution->EvaluateCDF( K2, 2 );
    }

  return p;
}

} // end of namespace Statistics
} // end of namespace itk

#endif

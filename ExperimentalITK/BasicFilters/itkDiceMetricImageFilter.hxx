/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDiceMetricImageFilter.hxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:52 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkDiceMetricImageFilter_hxx
#define _itkDiceMetricImageFilter_hxx

#include "itkDiceMetricImageFilter.h"

#include "itkImageRegionConstIterator.h"
#include "itkNumericTraits.h"

namespace itk {


template<class TInputImage1, class TInputImage2>
DiceMetricImageFilter<TInputImage1, TInputImage2>
::DiceMetricImageFilter()
{

  // this filter requires two input images
  this->SetNumberOfRequiredInputs( 2 );

  this->m_DiceMetric = NumericTraits<RealType>::Zero;
  this->m_PixelLabel = NumericTraits<InputImage1PixelType>::One;
}


template<class TInputImage1, class TInputImage2>
void
DiceMetricImageFilter<TInputImage1, TInputImage2>
::SetInput2( const TInputImage2 * image )
{
  this->SetNthInput( 1, const_cast<TInputImage2 *>( image ) );
}


template<class TInputImage1, class TInputImage2>
const typename DiceMetricImageFilter<TInputImage1, TInputImage2>
::InputImage2Type *
DiceMetricImageFilter<TInputImage1, TInputImage2>
::GetInput2()
{
  return static_cast< const TInputImage2 * >
    ( this->ProcessObject::GetInput( 1 ) );
}



template<class TInputImage1, class TInputImage2>
void
DiceMetricImageFilter<TInputImage1, TInputImage2>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();

  // this filter requires:
  // - the largeset possible region of the first image
  // - the corresponding region of the second image
  if ( this->GetInput1() )
    {
    InputImage1Pointer image1 =
      const_cast< InputImage1Type * >( this->GetInput1() );
    image1->SetRequestedRegionToLargestPossibleRegion();

    if ( this->GetInput2() )
      {
      InputImage2Pointer image2 =
        const_cast< InputImage2Type * >( this->GetInput2() );
      image2->SetRequestedRegion(
        this->GetInput1()->GetRequestedRegion() );
      }

    }
}


template<class TInputImage1, class TInputImage2>
void
DiceMetricImageFilter<TInputImage1, TInputImage2>
::EnlargeOutputRequestedRegion(DataObject *data)
{
  Superclass::EnlargeOutputRequestedRegion(data);
  data->SetRequestedRegionToLargestPossibleRegion();
}


template<class TInputImage1, class TInputImage2>
void
DiceMetricImageFilter<TInputImage1, TInputImage2>
::GenerateData()
{
  if ( this->GetInput1()->GetRequestedRegion().GetSize()
         != this->GetInput2()->GetRequestedRegion().GetSize() )
    {
    itkExceptionMacro( "Requested region sizes are not equal." );
    }

  unsigned long count1 = 0;
  unsigned long count2 = 0;
  unsigned long intersectionCount = 0;

  ImageRegionConstIterator<InputImage1Type> It1( this->GetInput1(),
    this->GetInput1()->GetRequestedRegion() );
  ImageRegionConstIterator<InputImage2Type> It2( this->GetInput2(),
    this->GetInput2()->GetRequestedRegion() );

  for ( It1.GoToBegin(), It2.GoToBegin(); !It1.IsAtEnd(); ++It1, ++It2 )
    {
    InputImage1PixelType label1 = It1.Get();
    InputImage1PixelType label2 = It2.Get();
    if ( label1 == this->m_PixelLabel )
      {
      count1++;
      }
    if ( label2 == this->m_PixelLabel )
      {
      count2++;
      }
    if ( label1 == label2 && label1 == this->m_PixelLabel )
      {
      intersectionCount++;
      }
    }

  this->m_DiceMetric = 2.0 * static_cast<RealType>( intersectionCount )
    / static_cast<RealType>( count1 + count2 );
}



template<class TInputImage1, class TInputImage2>
void
DiceMetricImageFilter<TInputImage1, TInputImage2>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "DiceMetric: "
     << m_DiceMetric << std::endl;
  os << indent << "PixelLabel: "
     << m_PixelLabel << std::endl;
}


}// end namespace itk
#endif

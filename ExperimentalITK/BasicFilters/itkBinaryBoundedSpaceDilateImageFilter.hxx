/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBinaryBoundedSpaceDilateImageFilter.hxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:49 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBinaryBoundedSpaceDilateImageFilter_hxx
#define __itkBinaryBoundedSpaceDilateImageFilter_hxx

#include "itkBinaryBoundedSpaceDilateImageFilter.h"

#include "itkCastImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

namespace itk
{

template <class TInputImage, class TOutputImage, class TKernel, class TBoundedSpaceImage>
BinaryBoundedSpaceDilateImageFilter<TInputImage, TOutputImage, TKernel, TBoundedSpaceImage>
::BinaryBoundedSpaceDilateImageFilter()
{
  this->m_Scaling = 1;
  this->m_BoundedSpaceValue = 1;
}

template< class TInputImage, class TOutputImage, class TKernel, class TBoundedSpaceImage>
void
BinaryBoundedSpaceDilateImageFilter< TInputImage, TOutputImage, TKernel, TBoundedSpaceImage>
::GenerateData()
{
  typedef ImageDuplicator<InputImageType> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage( this->GetInput() );
  duplicator->Update();

  typename InputImageType::ConstPointer input = duplicator->GetOutput();
  typename OutputImageType::Pointer output;

  for ( unsigned int i = 0; i < this->GetScaling(); i++ )
    {
    typename Superclass::Pointer dilater = Superclass::New();
    dilater->SetInput( input );
    dilater->SetKernel( this->GetKernel() );
    dilater->SetForegroundValue( this->GetForegroundValue() );
    dilater->SetBackgroundValue( this->GetBackgroundValue() );
    dilater->Update();

    bool isFilled = true;

    ImageRegionIterator<OutputImageType> ItD( dilater->GetOutput(),
      dilater->GetOutput()->GetRequestedRegion() );
    ImageRegionIterator<BoundedSpaceImageType> ItB( this->m_BoundedSpaceImage,
      this->m_BoundedSpaceImage->GetRequestedRegion() );
    ItD.GoToBegin();
    ItB.GoToBegin();
    while ( !ItD.IsAtEnd() )
      {
      if ( ItD.Get() == this->GetForegroundValue()
           && ItB.Get() != this->GetBoundedSpaceValue() )
        {
        ItD.Set( static_cast<typename OutputImageType::PixelType>( this->GetBackgroundValue() ) );
        }
      if ( ItB.Get() == this->m_BoundedSpaceValue &&
           ItD.Get() != this->GetForegroundValue() )
        {
        isFilled = false;
        }
      ++ItD;
      ++ItB;
      }

    output = dilater->GetOutput();

    if ( isFilled )
      {
      break;
      }

    typedef CastImageFilter<OutputImageType, InputImageType> CasterType;
    typename CasterType::Pointer caster = CasterType::New();
    caster->SetInput( dilater->GetOutput() );
    caster->Update();

    input = caster->GetOutput();
    }

  this->GraftOutput( output );
}


/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TOutput, class TKernel, class TBoundedSpaceImage>
void
BinaryBoundedSpaceDilateImageFilter<TInputImage, TOutput, TKernel, TBoundedSpaceImage>
::PrintSelf( std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}

} // end namespace itk

#endif

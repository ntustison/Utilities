/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBinaryBoundedSpaceErodeImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:49 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBinaryBoundedSpaceErodeImageFilter_txx
#define __itkBinaryBoundedSpaceErodeImageFilter_txx

#include "itkBinaryBoundedSpaceErodeImageFilter.h"

#include "itkBinaryBoundedSpaceDilateImageFilter.h"

namespace itk
{
 
template <class TInputImage, class TOutputImage, class TKernel, class TBoundedSpaceImage>
BinaryBoundedSpaceErodeImageFilter<TInputImage, TOutputImage, TKernel, TBoundedSpaceImage>
::BinaryBoundedSpaceErodeImageFilter()
{
  this->m_Scaling = 1;
  this->m_BoundedSpaceValue = 1;
}

template< class TInputImage, class TOutputImage, class TKernel, class TBoundedSpaceImage>
void
BinaryBoundedSpaceErodeImageFilter< TInputImage, TOutputImage, TKernel, TBoundedSpaceImage>
::GenerateData()
{
  typedef ImageDuplicator<InputImageType> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage( this->GetInput() );
  duplicator->Update();

  typename InputImageType::ConstPointer input = duplicator->GetOutput();

  ImageRegionIterator<InputImageType> ItI( input,
    input->GetRequestedRegion() );   
  ImageRegionIterator<BoundedSpaceImageType> ItB( this->m_BoundedSpaceImage,
    this->m_BoundedSpaceImage->GetRequestedRegion() );   

  // Now perform bounded space erosion
  for ( ItI.GoToBegin(), ItB.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItB )
    {
    typename InputImageType::PixelType pixel = ItI.Get();
    if ( pixel == this->GetForegroundValue() )
      {
      ItI.Set( this->GetBackgroundValue() );
      }
    else if ( pixel == this->GetBackgroundValue() 
      && ItB.Get() == this->m_BoundedSpaceValue )
      {
      ItI.Set( this->GetForegroundValue() );
      }
    } 

  typedef BinaryBoundedSpaceDilateImageFilter<InputImageType, 
    OutputImageType, KernelType> DilaterType;
  typename DilaterType::Pointer dilater = DilaterType::New();
  dilater->SetInput( input );
  dilater->SetKernel( this->GetKernel() );
  dilater->SetSpacing( this->GetSpacing() );
  dilater->SetForegroundValue( this->GetForegroundValue() );
  dilater->SetBoundedSpaceImage( this->GetBoundedSpaceImage() );
  dilater->SetBoundedSpaceValue( this->GetBoundedSpaceValue() );
  dilater->Update();

  typename OutputImageType::Pointer output = dilater->GetOutput();

  ImageRegionIterator<OutputImageType> It( output,
    output->GetRequestedRegion() );   
  It.GoToBegin();
  ItB.GoToBegin();
  while ( !It.IsAtEnd() && !ItB.IsAtEnd() ) 
    {
    typename OutputImageType::PixelType pixel = It.Get();
    if ( pixel == this->m_ForegroundValue )
      {
      It.Set( this->m_BackgroundValue );
      }
    else if ( pixel == this->m_BackgroundValue 
      && ItB.Get() == this->m_BoundedSpaceValue )
      {
      It.Set( this->m_ForegroundValue );
      }
    ++It;
    ++ItB;
    } 

  this->GraftOutput( output );
}


/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TOutput, class TKernel, class TBoundedSpaceImage>
void
BinaryBoundedSpaceErodeImageFilter<TInputImage, TOutput, TKernel, TBoundedSpaceImage>
::PrintSelf( std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}

} // end namespace itk

#endif

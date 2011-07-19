/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: itkConvolutionImageFilter.hxx,v $
Language:  C++
Date:      $Date: 2008/12/01 17:51:58 $
Version:   $Revision: 1.2 $

Copyright (c) Insight Software Consortium. All rights reser
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for detail.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef _itkConvolutionImageFilter_hxx_
#define _itkConvolutionImageFilter_hxx_

// disable debug warnings in MS compiler
#ifdef _MSC_VER
#pragma warning(disable: 4786)
#endif
 
#include "itkConvolutionImageFilter.h"

#include "itkImageKernelOperator.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkConstNeighborhoodIterator.h"

#include "vnl/vnl_math.h"

namespace itk {

template<class TInputImage, class TImageKernel, class TOutputImage>
ConvolutionImageFilter<TInputImage, TImageKernel, TOutputImage>
::ConvolutionImageFilter()
{
  this->m_ImageKernel = NULL;
}

template<class TInputImage, class TImageKernel, class TOutputImage>
ConvolutionImageFilter<TInputImage, TImageKernel, TOutputImage>
::~ConvolutionImageFilter()
{
}

template<class TInputImage, class TImageKernel, class TOutputImage>
void 
ConvolutionImageFilter<TInputImage, TImageKernel, TOutputImage>
::GenerateData()
{
  if ( this->m_ImageKernel.IsNull() )
    {
    itkExceptionMacro( "Image kernel is not specified." );  
    }  
  this->GetOutput()->SetRegions( this->GetInput()->GetRequestedRegion() );
  this->GetOutput()->SetOrigin( this->GetInput()->GetOrigin() );
  this->GetOutput()->SetSpacing( this->GetInput()->GetSpacing() );
  this->GetOutput()->Allocate();

  typedef ConstNeighborhoodIterator<InputImageType> NeighborhoodIteratorType;
  typename NeighborhoodIteratorType::RadiusType radius;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    radius[i] = static_cast<unsigned int>( 
      0.5 * this->m_ImageKernel->GetLargestPossibleRegion().GetSize()[i] );
    } 

  typedef typename NeighborhoodAlgorithm
    ::ImageBoundaryFacesCalculator<InputImageType> FaceCalculatorType;
  FaceCalculatorType faceCalculator;

  NeighborhoodInnerProduct<InputImageType> innerProduct;

  ImageKernelOperator<typename ImageKernelType::PixelType, 
    ImageDimension> imageKernelOperator;
  imageKernelOperator.SetImageKernel( this->m_ImageKernel );
  imageKernelOperator.CreateToRadius( radius );

  typename FaceCalculatorType::FaceListType faceList = faceCalculator( 
    this->GetInput(), this->GetInput()->GetRequestedRegion(), radius );
  typename FaceCalculatorType::FaceListType::iterator fit;  

  for ( fit = faceList.begin(); fit != faceList.end(); ++fit )
    {
    NeighborhoodIteratorType inIt( radius, this->GetInput(), *fit );  
    ImageRegionIterator<OutputImageType> outIt( this->GetOutput(), *fit );

    for ( inIt.GoToBegin(), outIt.GoToBegin(); !inIt.IsAtEnd(); ++inIt, ++outIt )
      {
      outIt.Set( innerProduct( inIt, imageKernelOperator ) );
      }
    }
}
}
#endif

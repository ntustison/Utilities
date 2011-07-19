/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAdaptiveFiltering3DImageFilter.hxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:49 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkAdaptiveFiltering3DImageFilter_hxx_
#define __itkAdaptiveFiltering3DImageFilter_hxx_

#include "itkAdaptiveFiltering3DImageFilter.h"

#include "itkMultiplyImageFilter.h"

#include "itkBinaryMagnitudeImageFilter.h"
#include "itkComplexToRealImageFilter.h"
#include "itkComplexToImaginaryImageFilter.h"
#include "itkFFTShiftImageFilter.h"
#include "itkGaussianImageSource.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVnlFFTComplexConjugateToRealImageFilter.h"
#include "itkVnlFFTRealToComplexConjugateImageFilter.h"

namespace itk {

/**
 * Constructor
 */
template <class TInputImage, class TOutputImage>
AdaptiveFiltering3DImageFilter<TInputImage, TOutputImage>
::AdaptiveFiltering3DImageFilter()
{
  RealType a = 2.0;
  RealType b = 1.0 + vcl_sqrt( 5 );
  RealType c = 1.0 / vcl_sqrt( 10.0 + 2.0 * vcl_sqrt( 5.0 ) ); 

  this->m_Directions[0][0] = c * a;
  this->m_Directions[0][1] = 0.0;
  this->m_Directions[0][2] = c * b;

  this->m_Directions[1][0] = -c * a;
  this->m_Directions[1][1] = 0.0;
  this->m_Directions[1][2] = c * b;
    
  this->m_Directions[2][0] = c * b;
  this->m_Directions[2][1] = c * a;
  this->m_Directions[2][2] = 0.0;
    
  this->m_Directions[3][0] = c * b;
  this->m_Directions[3][1] = -c * a;
  this->m_Directions[3][2] = 0.0;
    
  this->m_Directions[4][0] = 0.0;
  this->m_Directions[4][1] = c * b;
  this->m_Directions[4][2] = c * a;
    
  this->m_Directions[5][0] = 0.0;
  this->m_Directions[5][1] = c * b;
  this->m_Directions[5][2] = -c * a;
}

template <class TInputImage, class TOutputImage>
void
AdaptiveFiltering3DImageFilter<TInputImage, TOutputImage>
::GenerateData()
{

  for ( unsigned int i = 0; i < 1; i++ )
    {
    std::cout << i << std::endl;
 
    // Calculate the quadrature filter kernel 
/*   
    typename QuadratureKernelType::Pointer quadratureKernel = 
      QuadratureKernelType::New();
    quadratureKernel->SetDirection( this->m_Directions[i] );
    quadratureKernel->SetSize( this->GetInput()->GetRequestedRegion().GetSize() );
    quadratureKernel->SetOrigin( this->GetInput()->GetOrigin() );
    quadratureKernel->SetSpacing( this->GetInput()->GetSpacing() );
    quadratureKernel->Update();
*/
    typedef GaussianImageSource<OutputImageType> GaussianKernelType;
    typename GaussianKernelType::Pointer gaussianKernel = GaussianKernelType::New();
    gaussianKernel->SetSize( this->GetInput()->GetRequestedRegion().GetSize() );
    gaussianKernel->SetOrigin( this->GetInput()->GetOrigin() );
    gaussianKernel->SetSpacing( this->GetInput()->GetSpacing() ); 
    gaussianKernel->SetNormalized( false );

    typename GaussianKernelType::ArrayType sigma;
    sigma.Fill( 5.0 );
    gaussianKernel->SetSigma( sigma );
    typename GaussianKernelType::ArrayType mean;
    for ( unsigned int j = 0; j < ImageDimension; j++ )
      {
      mean[j] = this->GetInput()->GetOrigin()[j] + 0.5 
        * ( this->GetInput()->GetSpacing()[j] * this->GetInput()->GetLargestPossibleRegion().GetSize()[j ] );
      } 
    gaussianKernel->SetMean( mean );
    gaussianKernel->Update();

/*
    typedef FFTShiftImageFilter<OutputImageType, OutputImageType> ShifterType;
    typename ShifterType::Pointer shifter = ShifterType::New();
    shifter->SetInput( gaussianKernel->GetOutput() );
    shifter->SetInverse( false );
    shifter->Update();
*/
    // Calculate FT of the input image

    typedef VnlFFTRealToComplexConjugateImageFilter
      <typename OutputImageType::PixelType, ImageDimension> FFTFilterType;
    typename FFTFilterType::Pointer fftFilter = FFTFilterType::New();
    fftFilter->SetInput( this->GetInput() );
    fftFilter->Update(); 

    // Multiply the filtered image by the quadrature kernel.

    typedef typename FFTFilterType::OutputImageType ComplexImageType;
    typedef MultiplyImageFilter<ComplexImageType, OutputImageType, ComplexImageType>
      MultiplyImageFilterType;

    typename OutputImageType::Pointer ones = OutputImageType::New();
    ones->SetRegions( this->GetInput()->GetRequestedRegion().GetSize() );
    ones->SetOrigin( this->GetInput()->GetOrigin() );
    ones->SetSpacing( this->GetInput()->GetSpacing() ); 
    ones->Allocate();
    ones->FillBuffer( 1 );

    typename OutputImageType::SizeType size = this->GetInput()->GetLargestPossibleRegion().GetSize();

    ImageRegionIteratorWithIndex<OutputImageType> 
      It( ones, ones->GetLargestPossibleRegion() );
    for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      typename OutputImageType::IndexType idx = It.GetIndex();
      RealType radius = 0.0;
      for ( unsigned int j = 0; j < ImageDimension; j++ )
        {
        radius += ( idx[j] - 0.5 * size[j] ) * ( idx[j] - 0.5 * size[j] ); 
        }
      radius = vcl_sqrt( radius );
      if ( radius > 10 )
        {
        It.Set( 0 );
        } 
      }
    typedef FFTShiftImageFilter<ComplexImageType, ComplexImageType> ShifterType;
    typename ShifterType::Pointer shifter = ShifterType::New();
    shifter->SetInput( fftFilter->GetOutput() );
    shifter->SetInverse( false );
    shifter->Update();

    typename MultiplyImageFilterType::Pointer multiplier = MultiplyImageFilterType::New();
    multiplier->SetInput1( fftFilter->GetOutput() );
    multiplier->SetInput2( ones );
    multiplier->Update();

    typename ShifterType::Pointer inverseShifter = ShifterType::New();
    inverseShifter->SetInput( multiplier->GetOutput() );
    inverseShifter->SetInverse( false );
    inverseShifter->Update();


    typedef ComplexToRealImageFilter<ComplexImageType, OutputImageType> RealImageFilterType;
    typename RealImageFilterType::Pointer realImage = RealImageFilterType::New();
    realImage->SetInput( inverseShifter->GetOutput() );
    realImage->Update();    

    typedef ComplexToImaginaryImageFilter
      <ComplexImageType, OutputImageType> ImaginaryImageFilterType;
    typename ImaginaryImageFilterType::Pointer imaginaryImage = ImaginaryImageFilterType::New();
    imaginaryImage->SetInput( inverseShifter->GetOutput() );
    imaginaryImage->Update();

    typedef BinaryMagnitudeImageFilter<OutputImageType, OutputImageType, OutputImageType>
      MagnitudeImageFilterType;
    typename MagnitudeImageFilterType::Pointer magnitude = MagnitudeImageFilterType::New();
    magnitude->SetInput1( realImage->GetOutput() );
    magnitude->SetInput2( imaginaryImage->GetOutput() );
    magnitude->Update();

    typedef VnlFFTComplexConjugateToRealImageFilter
      <typename OutputImageType::PixelType, ImageDimension> InverseFFTFilterType;
    typename InverseFFTFilterType::Pointer inversefftFilter = InverseFFTFilterType::New();
    inversefftFilter->SetInput( inverseShifter->GetOutput() );
    inversefftFilter->Update(); 
    

/*
    // Calculate the magnitude image 

    typedef ComplexToRealImageFilter<ComplexImageType, OutputImageType> RealImageFilterType;
    typename RealImageFilterType::Pointer realImage = RealImageFilterType::New();
    realImage->SetInput( multiplier->GetOutput() );
    realImage->Update();
  
    typedef ComplexToImaginaryImageFilter
      <ComplexImageType, OutputImageType> ImaginaryImageFilterType;
    typename ImaginaryImageFilterType::Pointer imaginaryImage = ImaginaryImageFilterType::New();
    imaginaryImage->SetInput( multiplier->GetOutput() );
    imaginaryImage->Update();

    typedef BinaryMagnitudeImageFilter<OutputImageType, OutputImageType, OutputImageType>
      MagnitudeImageFilterType;
    typename MagnitudeImageFilterType::Pointer magnitude = MagnitudeImageFilterType::New();
    magnitude->SetInput1( realImage->GetOutput() );
    magnitude->SetInput2( imaginaryImage->GetOutput() );
    magnitude->Update();
*/
   
//    this->GraftNthOutput( i, magnitude->GetOutput() );

    this->GraftNthOutput( i, inversefftFilter->GetOutput() );
    }
} 
 
template <class TInputImage, class TOutputImage>
void
AdaptiveFiltering3DImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

}// end namespace itk

#endif

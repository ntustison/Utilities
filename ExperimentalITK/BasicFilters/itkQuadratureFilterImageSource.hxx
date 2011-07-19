/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkQuadratureFilterImageSource.hxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:53 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkQuadratureFilterImageSource_hxx
#define _itkQuadratureFilterImageSource_hxx

#include "itkQuadratureFilterImageSource.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkProgressReporter.h"
#include "itkObjectFactory.h"
 
namespace itk
{

template <class TOutputImage>
QuadratureFilterImageSource<TOutputImage>
::QuadratureFilterImageSource()
{
  //Initial image is 7 wide in each direction.
  for (unsigned int i = 0; i < ImageDimension; i++)
    {
    this->m_Size[i] = 7;
    this->m_Spacing[i] = 1.0;
    this->m_Origin[i] = 0.0;
    }

  this->m_CenterFrequency = 1.11;
  this->m_RelativeBandwidth = 2.0;

  RealType a = 2.0;
  RealType b = 1.0 + vcl_sqrt( 5 );
  RealType c = 1.0 / vcl_sqrt( 10.0 + 2.0 * vcl_sqrt( 5.0 ) ); 
  this->m_Direction[0] = c * a;
  this->m_Direction[1] = 0.0;
  this->m_Direction[2] = c * b;
}

template <class TOutputImage>
QuadratureFilterImageSource<TOutputImage>
::~QuadratureFilterImageSource()
{
}

//----------------------------------------------------------------------------
template <class TOutputImage>
void 
QuadratureFilterImageSource<TOutputImage>
::GenerateOutputInformation()
{
  OutputImageType *output;
  typename OutputImageType::IndexType index = { { 0 } };
   
  output = this->GetOutput( 0 );

  typename OutputImageType::RegionType largestPossibleRegion;
  largestPossibleRegion.SetSize( this->m_Size );
  largestPossibleRegion.SetIndex( index );
  output->SetLargestPossibleRegion( largestPossibleRegion );

  output->SetSpacing( this->m_Spacing );
  output->SetOrigin( this->m_Origin );
}

template <class TOutputImage>
void 
QuadratureFilterImageSource<TOutputImage>
::GenerateData()
{
  typename OutputImageType::Pointer outputPtr = this->GetOutput();

  // allocate the output buffer
  outputPtr->SetBufferedRegion( outputPtr->GetRequestedRegion() );
  outputPtr->Allocate();
  outputPtr->FillBuffer( 0 );

  // Create an iterator that will walk the output region
  ImageRegionIteratorWithIndex<OutputImageType>
    outIt( outputPtr, outputPtr->GetRequestedRegion() );

  // The position at which the function is evaluated
  Point<double, ImageDimension> evalPoint;

  ProgressReporter progress( this, 0, 
    outputPtr->GetRequestedRegion().GetNumberOfPixels() );

  RealType expFactor = 4.0 * vcl_log( 2.0 ) / vnl_math_sqr( this->m_RelativeBandwidth );
  PointType centerPoint;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    centerPoint[i] = 0.5 * ( this->m_Origin[i] + ( this->m_Size[i] - 1 ) * this->m_Spacing[i] );
    } 

  // Walk the output image, evaluating the spatial function at each pixel
  for ( outIt.GoToBegin(); !outIt.IsAtEnd(); ++outIt )
    {
    typename OutputImageType::IndexType index = outIt.GetIndex();
    outputPtr->TransformIndexToPhysicalPoint( index, evalPoint );

    VectorType frequency = evalPoint - centerPoint;
    RealType frequencyNorm = frequency.GetNorm(); 
    
    RealType dotProduct = ( frequency / frequencyNorm ) * this->m_Direction;
    if ( dotProduct <= 0 )
      {
      continue;
      }

    double value = vnl_math_sqr( dotProduct ) 
      * vcl_exp( expFactor * vnl_math_sqr( vcl_log( frequencyNorm / this->m_CenterFrequency ) ) );

    // Set the pixel value to the function value
    outIt.Set( static_cast<PixelType>( value ) );
    progress.CompletedPixel();
    }
}

template <class TOutputImage>
void 
QuadratureFilterImageSource<TOutputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os,indent );

  os << indent << "Image parameters: " << std::endl;
  os << indent << "  Size: " << this->m_Size << std::endl;
  os << indent << "  Origin: " << this->m_Origin << std::endl;
  os << indent << "  Spacing: " << this->m_Spacing << std::endl;

  os << indent << "QuadratureFilter filter parameters: " << std::endl;
  os << indent << "  Relative Bandwidth: " << this->m_RelativeBandwidth << std::endl;
  os << indent << "  Center Frequency: " << this->m_CenterFrequency << std::endl;
}

} // end namespace itk

#endif

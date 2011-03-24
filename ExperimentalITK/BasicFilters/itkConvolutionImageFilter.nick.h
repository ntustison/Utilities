/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkConvolutionImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/12/01 17:45:39 $
  Version:   $Revision: 1.2 $

  Copyright ( c ) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef _itkConvolutionImageFilter_h_
#define _itkConvolutionImageFilter_h_

#include "itkImageToImageFilter.h"

namespace itk {
template<class TInputImage, 
  class TImageKernel = TInputImage, class TOutputImage = TInputImage> 
class ITK_EXPORT ConvolutionImageFilter 
: public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  typedef ConvolutionImageFilter                               Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage>        Superclass;
  typedef SmartPointer<Self>                                   Pointer;
  typedef SmartPointer<const Self>                             ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );  
  
  /** Run-time type information ( and related methods ) */
  itkTypeMacro( ConvolutionImageFilter, ImageToImageFilter );

  /** Dimensionality of input and output data is assumed to be the same. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TInputImage::ImageDimension );

  
  typedef TInputImage                                          InputImageType;
  typedef TImageKernel                                         ImageKernelType;
  typedef TOutputImage                                         OutputImageType;

  itkSetObjectMacro( ImageKernel, ImageKernelType );

protected:
  /** de/constructor */
  ConvolutionImageFilter(); 
  ~ConvolutionImageFilter();

  void PrintSelf( std::ostream& os, Indent indent ) const {}
  void GenerateData();

private: 
  ConvolutionImageFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented
 
private :

  typename ImageKernelType::Pointer                            m_ImageKernel;

};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkConvolutionImageFilter.txx"
#endif

#endif


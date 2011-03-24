/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBinaryBoundedSpaceDilateImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:49 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBinaryBoundedSpaceDilateImageFilter_h
#define __itkBinaryBoundedSpaceDilateImageFilter_h

#include "itkBinaryDilateImageFilter.h"

namespace itk
{
/**
 * \class BinaryBoundedSpaceDilateImageFilter
 * \brief Fast binary erosion
 *
 * BinaryBoundedSpaceDilateImageFilter is a binary dilation
 * morphologic operation over a bounded space. 
 *
 * \sa ImageToImageFilter BinaryDilateImageFilter BinaryMorphologyImageFilter
 */
template <class TInputImage, class TOutputImage, 
  class TKernel, class TBoundedSpaceImage = TInputImage>
class ITK_EXPORT BinaryBoundedSpaceDilateImageFilter :
    public BinaryDilateImageFilter<TInputImage, TOutputImage, TKernel>
{
public:
  /** Extract dimension from input and output image. */
  itkStaticConstMacro( InputImageDimension, unsigned int,
                       TInputImage::ImageDimension );
  itkStaticConstMacro( OutputImageDimension, unsigned int,
                       TOutputImage::ImageDimension );

  /** Extract the dimension of the kernel */
  itkStaticConstMacro( KernelDimension, unsigned int,
                       TKernel::NeighborhoodDimension );
  
  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage  InputImageType;
  typedef TOutputImage OutputImageType;
  typedef TKernel      KernelType;

  /** Standard class typedefs. */
  typedef BinaryBoundedSpaceDilateImageFilter     Self;
  typedef BinaryDilateImageFilter
    <InputImageType, OutputImageType, KernelType> Superclass;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( BinaryBoundedSpaceDilateImageFilter, 
    BinaryDilateImageFilter )

  /** Image typedef support. */
  typedef TBoundedSpaceImage                         BoundedSpaceImageType;
  typedef typename BoundedSpaceImageType::PixelType  BoundedSpaceValueType;

  itkSetMacro( BoundedSpaceImage, typename BoundedSpaceImageType::Pointer );
  itkGetConstMacro( BoundedSpaceImage, typename BoundedSpaceImageType::Pointer ); 

  itkSetMacro( BoundedSpaceValue, BoundedSpaceValueType );
  itkGetConstMacro( BoundedSpaceValue, BoundedSpaceValueType ); 

  itkSetMacro( Scaling, unsigned int );
  itkGetConstMacro( Scaling, unsigned int );

protected:
  BinaryBoundedSpaceDilateImageFilter();
  virtual ~BinaryBoundedSpaceDilateImageFilter(){}
  void PrintSelf(std::ostream& os, Indent indent) const;

  void GenerateData();

  // type inherited from the superclass
  typedef typename Superclass::NeighborIndexContainer NeighborIndexContainer;

private:
  BinaryBoundedSpaceDilateImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  typename BoundedSpaceImageType::Pointer                 m_BoundedSpaceImage;
  BoundedSpaceValueType                                   m_BoundedSpaceValue;

  unsigned int                                             m_Scaling;
  

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBinaryBoundedSpaceDilateImageFilter.txx"
#endif

#endif

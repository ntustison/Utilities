/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkScalarToFractalImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/11/16 02:05:17 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkScalarToFractalImageFilter_h
#define __itkScalarToFractalImageFilter_h

#include "itkImageToImageFilter.h"

#include "itkConstNeighborhoodIterator.h"

namespace itk {

/** \class ScalarToFractalImageFilter

 */
template<class TInputImage, class TOutputImage>
class ITK_EXPORT ScalarToFractalImageFilter :
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef ScalarToFractalImageFilter                      Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage>   Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;

  /** Standard New method. */
  itkNewMacro( Self );

  /** ImageDimension constants */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TInputImage::ImageDimension );

  /** Some convenient typedefs. */
  typedef float                                   RealType;
  typedef TInputImage                             InputImageType;
  typedef TOutputImage                            OutputImageType;
  typedef Image<char, ImageDimension>             MaskImageType;


  /** Runtime information support. */
  itkTypeMacro( ScalarToFractalImageFilter,
                ImageToImageFilter );

  typedef ConstNeighborhoodIterator<InputImageType>
    ConstNeighborhoodIteratorType;
  typedef typename ConstNeighborhoodIteratorType::RadiusType RadiusType;

  itkSetMacro( NeighborhoodRadius, RadiusType );
  itkGetConstMacro( NeighborhoodRadius, RadiusType );

  itkSetMacro( MaskImage, typename MaskImageType::Pointer );
  itkGetConstMacro( MaskImage, typename MaskImageType::Pointer );

protected:
  ScalarToFractalImageFilter();
  ~ScalarToFractalImageFilter() {};
  void PrintSelf( std::ostream& os, Indent indent ) const;

  void GenerateData();

private:
  ScalarToFractalImageFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  RadiusType                       m_NeighborhoodRadius;
  typename MaskImageType::Pointer  m_MaskImage;

}; // end of class

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkScalarToFractalImageFilter.txx"
#endif

#endif

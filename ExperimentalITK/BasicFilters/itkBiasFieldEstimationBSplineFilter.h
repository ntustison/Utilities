/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBiasFieldEstimationBSplineFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:49 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBiasFieldEstimationBSplineFilter_h
#define __itkBiasFieldEstimationBSplineFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{

template <class TInputImage, class TOutputImage = TInputImage>
class BiasFieldEstimationBSplineFilter :
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef BiasFieldEstimationBSplineFilter   Self;
  typedef ImageToImageFilter<
                       TInputImage, 
                       TOutputImage>        Superclass;

  typedef SmartPointer<Self>                   Pointer;
  typedef SmartPointer<const Self>             ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Extract dimension from input and output image. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TInputImage::ImageDimension );

  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage                           InputImageType;
  typedef TOutputImage                          OutputImageType;

  typedef typename InputImageType::PixelType    InputPixelType;
  typedef typename OutputImageType::PixelType   OutputPixelType;

  typedef double                                RealType;
  typedef Image<RealType, ImageDimension>       RealImageType;

  typedef Vector<RealType, 1>                   VectorType;

  itkSetMacro( ConfidenceImage, typename RealImageType::Pointer );
  itkGetConstMacro( ConfidenceImage, typename RealImageType::Pointer );

  itkSetMacro( MinimumLevel, unsigned int );
  itkGetConstMacro( MinimumLevel, unsigned int );

  itkSetMacro( MaximumLevel, unsigned int );
  itkGetConstMacro( MaximumLevel, unsigned int );

  itkSetMacro( SplineOrder, unsigned int );
  itkGetConstMacro( SplineOrder, unsigned int );

  itkSetMacro( IgnorePixelValue, InputPixelType );
  itkGetConstReferenceMacro( IgnorePixelValue, InputPixelType );
  
protected:

  BiasFieldEstimationBSplineFilter();
  virtual ~BiasFieldEstimationBSplineFilter();
  
  void PrintSelf( std::ostream& os, Indent indent ) const;

  void GenerateData();

private:
   typename RealImageType::Pointer              m_ConfidenceImage;
 
   InputPixelType                               m_IgnorePixelValue;
   unsigned int                                 m_MinimumLevel;
   unsigned int                                 m_MaximumLevel;
   unsigned int                                 m_SplineOrder;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBiasFieldEstimationBSplineFilter.txx"
#endif

#endif

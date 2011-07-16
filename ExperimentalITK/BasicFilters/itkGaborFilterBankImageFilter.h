/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGaborFilterBankImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:52 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGaborFilterBankImageFilter_h
#define __itkGaborFilterBankImageFilter_h

#include "itkImageToImageFilter.h"

#include "itkFixedArray.h"
#include "itkGaussianImageSource.h"


namespace itk
{

/** \class GaborFilterBankImageFilter.h
 * \brief Image filter.
 */

template <class TInputImage, class TOutputImage>
class GaborFilterBankImageFilter 
: public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  typedef GaborFilterBankImageFilter                   Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage>       Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Extract dimension from input image. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TInputImage::ImageDimension );
        
  /** Image typedef support. */
  typedef TInputImage                             InputImageType;
  typedef TOutputImage                            OutputImageType;

  /** Other typedef */
  typedef float                                   RealType;
  typedef Image<RealType, 
    itkGetStaticConstMacro( ImageDimension )>     RealImageType; 
  typedef Image<unsigned int, 
    itkGetStaticConstMacro( ImageDimension )>     LabelImageType;   
  typedef GaussianImageSource<RealImageType>      FourierTransformGaborImageSourceType; 
  typedef FixedArray<RealType, 
    itkGetStaticConstMacro( ImageDimension )>     ArrayType; 
  typedef FixedArray<unsigned int, 
    itkGetStaticConstMacro( ImageDimension )>     UnsignedIntArrayType; 

  /** Helper functions */

  itkSetMacro( NumberOfGaborSpacingSteps, unsigned int );
  itkGetConstMacro( NumberOfGaborSpacingSteps, unsigned int );

  itkSetMacro( GaborSpacingMinimum, RealType );
  itkGetConstMacro( GaborSpacingMinimum, RealType );

  itkSetMacro( GaborSpacingMaximum, RealType );
  itkGetConstMacro( GaborSpacingMaximum, RealType );

  itkSetMacro( RotationAngleMinimum, ArrayType );
  itkGetConstMacro( RotationAngleMinimum, ArrayType );

  itkSetMacro( RotationAngleMaximum, ArrayType );
  itkGetConstMacro( RotationAngleMaximum, ArrayType );

  itkSetMacro( NumberOfRotationAngleSteps, UnsignedIntArrayType );
  itkGetConstMacro( NumberOfRotationAngleSteps, UnsignedIntArrayType );

protected:
  GaborFilterBankImageFilter();
  virtual ~GaborFilterBankImageFilter();
  void PrintSelf( std::ostream& os, Indent indent ) const;

  void GenerateData();

private:
  GaborFilterBankImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  unsigned int                                     m_NumberOfGaborSpacingSteps;
  RealType                                         m_GaborSpacingMinimum;
  RealType                                         m_GaborSpacingMaximum;

  UnsignedIntArrayType                             m_NumberOfRotationAngleSteps;
  ArrayType                                        m_RotationAngleMinimum;
  ArrayType                                        m_RotationAngleMaximum;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGaborFilterBankImageFilter.hxx"
#endif

#endif


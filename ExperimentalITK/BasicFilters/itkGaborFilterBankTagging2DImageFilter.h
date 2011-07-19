/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGaborFilterBankTagging2DImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:52 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGaborFilterBankTagging2DImageFilter_h
#define __itkGaborFilterBankTagging2DImageFilter_h

#include "itkImageToImageFilter.h"

#include "itkFixedArray.h"
#include "itkGaussianImageSource.h"


namespace itk
{

/** \class GaborFilterBankTagging2DImageFilter.h
 * \brief Image filter.
 */

template <class TInputImage, class TOutputImage>
class GaborFilterBankTagging2DImageFilter 
: public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  typedef GaborFilterBankTagging2DImageFilter                   Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage>       Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Extract dimension from input image. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TInputImage::ImageDimension );
        
  /** Image typedef support. */
  typedef TInputImage                                         InputImageType;
  typedef TOutputImage                                        OutputImageType;

  /** Other typedef */
  typedef float                                               RealType;
  typedef Image<RealType, 
    itkGetStaticConstMacro( ImageDimension )>                 RealImageType; 
  typedef Image<unsigned int, 
    itkGetStaticConstMacro( ImageDimension )>                 LabelImageType;   
  typedef GaussianImageSource<RealImageType>                  FourierTransformGaborImageSourceType; 
  typedef FixedArray<RealType, 
    itkGetStaticConstMacro( ImageDimension )>                 ArrayType; 
  typedef FixedArray<unsigned int, 
    itkGetStaticConstMacro( ImageDimension )>                 UnsignedIntArrayType; 

  /** Helper functions */

  itkSetMacro( NumberOfTagSpacingSteps, unsigned int );
  itkGetConstMacro( NumberOfTagSpacingSteps, unsigned int );

  itkSetMacro( TagSpacingMinimum, RealType );
  itkGetConstMacro( TagSpacingMinimum, RealType );

  itkSetMacro( TagSpacingMaximum, RealType );
  itkGetConstMacro( TagSpacingMaximum, RealType );

  itkSetMacro( RotationAngleMinimum, ArrayType );
  itkGetConstMacro( RotationAngleMinimum, ArrayType );

  itkSetMacro( RotationAngleMaximum, ArrayType );
  itkGetConstMacro( RotationAngleMaximum, ArrayType );

  itkSetMacro( NumberOfRotationAngleSteps, UnsignedIntArrayType );
  itkGetConstMacro( NumberOfRotationAngleSteps, UnsignedIntArrayType );

  itkSetClampMacro( ThresholdPercentage, RealType, 0.0001, 0.9999 );
  itkGetConstMacro( ThresholdPercentage, RealType );

  itkSetMacro( UseAutomaticThresholding, bool );
  itkGetConstMacro( UseAutomaticThresholding, bool );
  itkBooleanMacro( UseAutomaticThresholding );

  itkSetMacro( MaskImage, typename LabelImageType::Pointer );
  itkGetConstMacro( MaskImage, typename LabelImageType::Pointer );

  itkGetConstMacro( MaximumResponseImage, typename RealImageType::Pointer );

protected:
  GaborFilterBankTagging2DImageFilter();
  virtual ~GaborFilterBankTagging2DImageFilter();
  void PrintSelf( std::ostream& os, Indent indent ) const;

  void GenerateData();

private:
  GaborFilterBankTagging2DImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  typename LabelImageType::Pointer                           m_MaskImage;

  typename RealImageType::Pointer                            m_MaximumResponseImage;

  unsigned int                                               m_NumberOfTagSpacingSteps;
  RealType                                                   m_TagSpacingMinimum;
  RealType                                                   m_TagSpacingMaximum;

  UnsignedIntArrayType                                       m_NumberOfRotationAngleSteps;
  ArrayType                                                  m_RotationAngleMinimum;
  ArrayType                                                  m_RotationAngleMaximum;

  RealType                                                   m_ThresholdPercentage;
  bool                                                       m_UseAutomaticThresholding;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGaborFilterBankTagging2DImageFilter.hxx"
#endif

#endif


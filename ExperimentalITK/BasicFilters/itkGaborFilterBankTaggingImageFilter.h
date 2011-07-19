/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGaborFilterBankTaggingImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:52 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGaborFilterBankTaggingImageFilter_h
#define __itkGaborFilterBankTaggingImageFilter_h

#include "itkImageToImageFilter.h"

#include "itkFixedArray.h"
#include "itkGaussianImageSource.h"


namespace itk
{

/** \class GaborFilterBankTaggingImageFilter.h
 * \brief Image filter.
 */

template <class TInputImage, class TOutputImage>
class GaborFilterBankTaggingImageFilter 
: public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  typedef GaborFilterBankTaggingImageFilter                   Self;
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
  GaborFilterBankTaggingImageFilter();
  virtual ~GaborFilterBankTaggingImageFilter();
  void PrintSelf( std::ostream& os, Indent indent ) const;

  void GenerateData();

private:
  GaborFilterBankTaggingImageFilter(const Self&); //purposely not implemented
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
#include "itkGaborFilterBankTaggingImageFilter.hxx"
#endif

#endif


/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMultipleLabelToDistanceMapImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:52 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMultipleLabelToDistanceMapImageFilter_h
#define __itkMultipleLabelToDistanceMapImageFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{

/** \class MultipleLabelToDistanceMapImageFilter.h
 * \brief Image filter.
 */

template <class TInputImage, class TOutputImage>
class MultipleLabelToDistanceMapImageFilter 
: public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  typedef MultipleLabelToDistanceMapImageFilter               Self;
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

  /** Image typedef support. */
  typedef typename InputImageType::PixelType                  InputPixelType;
  typedef typename OutputImageType::PixelType                 OutputPixelType;

  typedef typename InputImageType::SizeType                   InputSizeType;
  typedef typename OutputImageType::SizeType                  OutputSizeType;

  typedef typename InputImageType::IndexType                  InputIndexType;
  typedef typename OutputImageType::IndexType                 OutputIndexType;

  typedef typename InputImageType::SpacingType                InputSpacingType;
  typedef typename OutputImageType::SpacingType               OutputSpacingType;

  typedef float                                               RealType;
  typedef Image<RealType, 
    itkGetStaticConstMacro( ImageDimension )>                 RealImageType;


  /** Set/Get if the distance should be squared. */
  itkSetMacro( SquaredDistance, bool );
  itkGetConstMacro( SquaredDistance, bool );
  itkBooleanMacro( SquaredDistance );

  /** Set/Get if image spacing should be used in computing distances. */
  itkSetMacro( UseImageSpacing, bool );
  itkGetConstMacro( UseImageSpacing, bool );
  itkBooleanMacro( UseImageSpacing );

  /** Set/Get  */
  itkSetMacro( NormalizeImage, bool );
  itkGetConstMacro( NormalizeImage, bool );
  itkBooleanMacro( NormalizeImage );

  /** Set/Get  */
  itkSetMacro( Sigma, RealType );
  itkGetConstMacro( Sigma, RealType );

protected:

  MultipleLabelToDistanceMapImageFilter ();

  virtual ~MultipleLabelToDistanceMapImageFilter();
  
  void PrintSelf( std::ostream& os, Indent indent ) const;

  void GenerateData();
  
private:

  MultipleLabelToDistanceMapImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  bool                                                        m_UseImageSpacing;
  bool                                                        m_SquaredDistance;

  bool                                                        m_NormalizeImage;
  RealType                                                    m_Sigma;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultipleLabelToDistanceMapImageFilter.hxx"
#endif

#endif


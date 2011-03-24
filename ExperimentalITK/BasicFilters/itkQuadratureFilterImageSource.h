/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkQuadratureFilterImageSource.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:53 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkQuadratureFilterImageSource_h
#define __itkQuadratureFilterImageSource_h

#include "itkImageSource.h"

#include "itkVector.h"

namespace itk
{

/** \class QuadratureFilterImageSource
 * \brief Generate an 3-dimensional quadrature filter.
 *
 * \ingroup DataSources
 */
template <typename TOutputImage>
class ITK_EXPORT QuadratureFilterImageSource : public ImageSource<TOutputImage>
{
public:

  /** Standard class typedefs. */
  typedef QuadratureFilterImageSource            Self;
  typedef ImageSource<TOutputImage>              Superclass;
  typedef SmartPointer<Self>                     Pointer;
  typedef SmartPointer<const Self>               ConstPointer;

  /** Output image typedefs */
  typedef TOutputImage                           OutputImageType;
  typedef typename OutputImageType::PixelType    PixelType;
  typedef typename OutputImageType::RegionType   RegionType;
  typedef typename OutputImageType::SpacingType  SpacingType;
  typedef typename OutputImageType::PointType    PointType;

  typedef typename RegionType::SizeType          SizeType;

  typedef double                                 RealType;

  /** Run-time type information (and related methods). */
  itkTypeMacro( QuadratureFilterImageSource, ImageSource );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );
  
  /** Dimensionality of the output image */
  itkStaticConstMacro( ImageDimension, unsigned int, 
                       OutputImageType::ImageDimension );

  typedef Vector<RealType, 
    itkGetStaticConstMacro( ImageDimension )>    VectorType;

  /** Gets and sets for filter parameters */
  itkSetMacro( Size, SizeType );
  itkGetConstMacro( Size, SizeType );
  
  itkSetMacro( Spacing, SpacingType );
  itkGetConstMacro( Spacing, SpacingType );

  itkSetMacro( Origin, PointType );
  itkGetConstMacro( Origin, PointType );

  itkSetMacro( RelativeBandwidth, RealType );
  itkGetConstMacro( RelativeBandwidth, RealType );
  
  itkSetMacro( CenterFrequency, RealType );
  itkGetConstMacro( CenterFrequency, RealType );

  itkSetMacro( Direction, VectorType );
  itkGetConstMacro( Direction, VectorType );

protected:
  QuadratureFilterImageSource();
  ~QuadratureFilterImageSource();
  void PrintSelf(std::ostream& os, Indent indent) const;
  void GenerateData();
  virtual void GenerateOutputInformation();

private:
  QuadratureFilterImageSource( const QuadratureFilterImageSource& ); //purposely not implemented
  void operator=( const QuadratureFilterImageSource& ); //purposely not implemented

  SizeType       m_Size;    //size of the output image
  SpacingType    m_Spacing; //spacing
  PointType      m_Origin;  //origin

  /** Parameters for the quadrature filter */
  VectorType                                     m_Direction;
  RealType                                       m_RelativeBandwidth;
  RealType                                       m_CenterFrequency;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkQuadratureFilterImageSource.txx"
#endif

#endif

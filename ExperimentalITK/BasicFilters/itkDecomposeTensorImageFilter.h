/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDecomposeTensorImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:51 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDecomposeTensorImageFilter_h
#define __itkDecomposeTensorImageFilter_h

#include "itkConstNeighborhoodIterator.h"
#include "itkImageToImageFilter.h"
#include "itkVariableSizeMatrix.h"
#include "itkVector.h"

namespace itk
{
/** \class DecomposeTensorImageFilter
 *
 */
template <typename TInputImage,
          typename TRealType = float,
          typename TOutputImage = Image<
            itk::VariableSizeMatrix<TRealType>, 
            ::itk::GetImageDimension<TInputImage>::ImageDimension>
>
class ITK_EXPORT DecomposeTensorImageFilter :
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef DecomposeTensorImageFilter Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods) */
  itkTypeMacro( DecomposeTensorImageFilter, ImageToImageFilter );
  
  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
  typedef TInputImage                                          InputImageType;
  typedef typename InputImageType::PixelType                   InputMatrixType;
  typedef TOutputImage                                         OutputImageType;
  typedef typename OutputImageType::PixelType                  OutputMatrixType;
  typedef typename OutputImageType::RegionType                 OutputImageRegionType;

  /** The dimensionality of the input and output images. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TInputImage::ImageDimension );
  
  /** Dimensions of the row/column matrix type of the input image. */
  itkStaticConstMacro( RowDimensions, unsigned int,
                       InputMatrixType::RowDimensions );
  itkStaticConstMacro( ColumnDimensions, unsigned int,
                       InputMatrixType::ColumnDimensions );

  /** Define the data type and the vector of data type used in calculations. */
  typedef TRealType                                            RealType;
  typedef Image<RealType, 
                itkGetStaticConstMacro( ImageDimension )>      RealImageType;
  typedef Vector<RealType,
                 itkGetStaticConstMacro( ImageDimension )>     VectorType;
  typedef Image<VectorType, 
                itkGetStaticConstMacro( ImageDimension )>      VectorImageType; 

  itkBooleanMacro( SymmetricTensors );
  itkSetMacro( SymmetricTensors, bool );
  itkGetConstReferenceMacro( SymmetricTensors, bool ); 

  itkSetClampMacro( WhichDecomposition, unsigned int, 0, 4 );
  itkGetConstReferenceMacro( WhichDecomposition, unsigned int );

  void SetEigenDecomposition()
    {
    this->SetWhichDecomposition( 0 );
    } 
  void SetRightPolarDecomposition()
    {
    this->SetWhichDecomposition( 1 );
    } 
  void SetLeftPolarDecomposition()
    {
    this->SetWhichDecomposition( 2 );
    } 
  void SetQRDecomposition()
    {
    this->SetWhichDecomposition( 3 );
    } 
  void SetSVDDecomposition()
    {
    this->SetWhichDecomposition( 4 );
    } 

  /**
   * Helper routines to produce eigenvalue and eigenvector images.
   * Eigenvectors and eigenvalues are sorted in ascending order
   * of the eigenvalue.
   */
  RealImageType* GetEigenValueImage( unsigned int );
  VectorImageType* GetEigenVectorImage( unsigned int );

protected:
  DecomposeTensorImageFilter();
  virtual ~DecomposeTensorImageFilter() {}

  void GenerateData();

  void PrintSelf ( std::ostream& os, Indent indent ) const;

private:
  bool m_SymmetricTensors;
  unsigned int m_WhichDecomposition;

  DecomposeTensorImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  void EigenDecomposition( InputMatrixType &, OutputMatrixType &, OutputMatrixType & );
  void RightPolarDecomposition( InputMatrixType &, OutputMatrixType &, OutputMatrixType & );
  void LeftPolarDecomposition( InputMatrixType &, OutputMatrixType &, OutputMatrixType & );
  void QRDecomposition( InputMatrixType &, OutputMatrixType &, OutputMatrixType & );
  void SVDDecomposition( InputMatrixType &, OutputMatrixType &, OutputMatrixType & , OutputMatrixType &);

  typename VectorImageType::Pointer m_EigenVectorImage;
  typename RealImageType::Pointer m_EigenValueImage;
};
  
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDecomposeTensorImageFilter.hxx"
#endif

#endif

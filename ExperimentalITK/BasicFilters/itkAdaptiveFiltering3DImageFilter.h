/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAdaptiveFiltering3DImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:49 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkAdaptiveFiltering3DImageFilter_h
#define __itkAdaptiveFiltering3DImageFilter_h

#include "itkImageToImageFilter.h"

#include "itkFixedArray.h"
#include "itkQuadratureFilterImageSource.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkSymmetricEigenVectorAnalysisImageFilter.h"

namespace itk {
/** \class AdaptiveFiltering3DFunction
 * 
 * \sa MultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter 
 * \sa AdaptiveFiltering3DImageFilter 
 * \ingroup FiniteDifferenceFunctions
 * \ingroup Functions
 */


template <class TInputImage, class TOutputImage>
class ITK_EXPORT AdaptiveFiltering3DImageFilter  
  : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs */
  typedef AdaptiveFiltering3DImageFilter                   Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage>    Superclass;
  typedef SmartPointer<Self>                               Pointer;
  typedef SmartPointer<const Self>                         ConstPointer;
 

  /** Method for creation through the object factory */
  itkNewMacro( Self );

  /** Run-time type information (and related methods) */
  itkTypeMacro( AdaptiveFiltering3DImageFilter, ImageToImageFilter );
  

  /** Convenient typedefs */
  typedef float                                            RealType;
  typedef SymmetricSecondRankTensor<RealType, 3>           TensorType;
  typedef Image<TensorType, 3>                             TensorImageType;
  typedef TInputImage                                      InputImageType;
  typedef TOutputImage                                     OutputImageType;

  /** 
   * Dimensionality of input and output data is 
   * assumed to be the same and equal to 3. 
   */
  itkStaticConstMacro(ImageDimension, 
    unsigned int, TInputImage::ImageDimension);

  typedef Matrix<RealType, 
    itkGetStaticConstMacro( ImageDimension ), 
    itkGetStaticConstMacro( ImageDimension )>              MatrixType;

  // Define image of matrix pixel type 
  typedef Image<MatrixType, 3>                             MatrixImageType;

   // Define the type for storing the eigen-value
  typedef FixedArray<RealType, 3>                          EigenValuesArrayType;
  
  // Declare the types of the output images
  typedef Image<EigenValuesArrayType, 3>                   EigenValuesImageType;
  
  // Declare the type for the filter
  typedef SymmetricEigenVectorAnalysisImageFilter
    <TensorImageType, 
     EigenValuesImageType, 
     MatrixImageType>                                      EigenAnalysisFilterType;

  typedef QuadratureFilterImageSource<OutputImageType>     QuadratureKernelType;
  typedef typename QuadratureKernelType::VectorType        DirectionType;


protected:
  AdaptiveFiltering3DImageFilter();
 ~AdaptiveFiltering3DImageFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  virtual void GenerateData(); 

private:
  //purposely not implemented
  AdaptiveFiltering3DImageFilter(const Self&); 
  void operator=(const Self&); //purposely not implemented

  DirectionType                                          m_Directions[6];
};
  

}// end namespace itk

#if ITK_TEMPLATE_TXX
# include "itkAdaptiveFiltering3DImageFilter.txx"
#endif

#endif

/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkNormalizedCutsSegmentationImageFilter.h,v $
  Language:  C++
  Date:      
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkNormalizedCutsSegmentationImageFilter_h_
#define __itkNormalizedCutsSegmentationImageFilter_h_

#include "itkImageToImageFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkSparseSymmetricMatrixEigenAnalysis.h"

#include "vnl/vnl_sparse_matrix.h"
#include "vnl/vnl_vector.h"

namespace itk 
{
/** \class itkNormalizedCutsSegmentationImageFilter
 * 
 * 
 * \par
 * This class constructs an n-dimensional graph from an image based on
 * the weighting .
 * 
 * \par
 * To construct a graph, we need to create nodes from the appropriate voxels as well
 * as the edges linking the source nodes and the target nodes.
 *
 * \par REFERENCE
 * Y. Boykov and V. Kolmogorov, "An Experimental Comparison of Min-Cut/Max-Flow 
 * Algorithms for Energy Minimization in Vision," IEEE-PAMI, 26(9):1124-1137, 2004.
 *
 * \par INPUT
 *
 *
 *  */

template<class TInputImage, class TLabelImage>
class ITK_EXPORT NormalizedCutsSegmentationImageFilter 
: public ImageToImageFilter<TInputImage, TLabelImage>
{
public:
  /** Standard "Self" typedef. */
  typedef NormalizedCutsSegmentationImageFilter              Self;
  typedef ImageToImageFilter<TInputImage, TLabelImage>       Superclass;
  typedef SmartPointer<Self>                                 Pointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro( Self );  
 
  /** Run-time type information (and related methods). */
  itkTypeMacro( NormalizedCutsSegmentationImageFilter, 
                ImageToImageFilter );

  /** Hold on to the type information specified by the template parameters. */
  typedef TInputImage                                        InputImageType;
  typedef typename InputImageType::PixelType                 InputPixelType;
  typedef typename InputImageType::PointType                 InputPointType;
  typedef typename InputImageType::IndexType                 IndexType;  
  typedef typename InputImageType::SizeType                  SizeType;

  /** Dimensionality of the input image */
  itkStaticConstMacro( ImageDimension, unsigned int, 
    TInputImage::ImageDimension );

  typedef TLabelImage                                        LabelImageType;
  typedef typename LabelImageType::PixelType                 LabelPixelType;
  
  typedef double                                             RealType;
  typedef Image<RealType, 
    itkGetStaticConstMacro(ImageDimension)>                  RealImageType;
  typedef ConstNeighborhoodIterator<InputImageType>          ConstNeighborhoodIteratorType;
  typedef typename ConstNeighborhoodIteratorType::RadiusType RadiusType;    

  typedef vnl_sparse_matrix<RealType>                        SparseMatrixType;
  typedef SparseSymmetricMatrixEigenAnalysis<RealType>       EigenSystemType;

  itkSetMacro( NumberOfClasses, unsigned int );  
  itkGetConstMacro( NumberOfClasses, unsigned int );

  itkSetMacro( NumberOfSplittingPoints, unsigned int );  
  itkGetConstMacro( NumberOfSplittingPoints, unsigned int );

  itkSetMacro( MaskImage, typename LabelImageType::Pointer );
  itkGetConstMacro( MaskImage, typename LabelImageType::Pointer );
  
  itkSetMacro( ForegroundValue, LabelPixelType );
  itkGetConstMacro( ForegroundValue, LabelPixelType );

  typename RealImageType::Pointer GetEigenImage( unsigned int );

//  itkSetMacro( PartitionStrategy, PartitionStrategyTypeEnumeration );
//  itkGetConstMacro( PartitionStrategy, PartitionStrategyTypeEnumeration );
  
protected:  
  NormalizedCutsSegmentationImageFilter(); 
  ~NormalizedCutsSegmentationImageFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  void GenerateData();
  virtual void GenerateOutputInformation(){}; // do nothing

  /** enum to indicate if the gradient image is specified as a single multi-
   * component image or as several separate images */
  typedef enum
    {
    Zero,
    Mean,
    SplittingPoints
    } PartitionStrategyTypeEnumeration;

private:
  NormalizedCutsSegmentationImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  void GenerateEigenSystemFromInputImage();
  void SolveEigensystem();
  void LabelOutputImage();
  RealType EvaluateWeightFunction( IndexType, IndexType );  
  void MulticlassSpectralClustering();
  void EigenVectorClustering();
  
  inline unsigned long 
  IndexToNumber( IndexType index, SizeType size )
    {
    IndexType k;     
    k[0] = 1;    
    
    for ( unsigned int i = 1; i < ImageDimension; i++ )
      {
      k[i] = size[ImageDimension-i-1]*k[i-1];
      }  
    unsigned long number = 0;
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      number += ( k[ImageDimension-i-1] * index[ImageDimension-i-1] );
      }
    return number;
    }

  inline IndexType 
  NumberToIndex( unsigned int number, SizeType size )
    {
    IndexType k;     
    k[0] = 1;    
    
    for ( unsigned int i = 1; i < ImageDimension; i++ )
      {
      k[i] = size[ImageDimension-i-1]*k[i-1];
      }  
    IndexType index;
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      index[ImageDimension-i-1] = static_cast<unsigned int>( number/k[ImageDimension-i-1] );
      number %= k[ImageDimension-i-1];
      }
    return index;
    }

  
  unsigned int                                               m_NumberOfClasses;
  RadiusType                                                 m_Radius;
  typename LabelImageType::Pointer                           m_MaskImage;
  LabelPixelType                                             m_ForegroundValue;
  RealType                                                   m_SpatialSigma;
  RealType                                                   m_DataSigma;

  PartitionStrategyTypeEnumeration                           m_PartitionStrategy;

  unsigned int                                               m_NumberOfSplittingPoints;
   
  SparseMatrixType                                           m_D;
  SparseMatrixType                                           m_E;
  SparseMatrixType                                           m_W;
  typename EigenSystemType::Pointer                          m_EigenSystem;    

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkNormalizedCutsSegmentationImageFilter.hxx"
#endif

#endif

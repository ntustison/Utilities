/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBinaryReinhardtMorphologicalImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:50 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBinaryReinhardtMorphologicalImageFilter_h
#define __itkBinaryReinhardtMorphologicalImageFilter_h

#include "itkBinaryMorphologyImageFilter.h"

namespace itk
{
/**
 * \class BinaryReinhardtMorphologicalImageFilter
 * \brief Fast binary erosion
 *
 * BinaryReinhardtMorphologicalImageFilter . 
 *
 * \sa ImageToImageFilter BinaryDilateImageFilter BinaryMorphologyImageFilter
 */
template <class TInputImage, class TOutputImage, class TKernel>
class ITK_EXPORT BinaryReinhardtMorphologicalImageFilter :
    public BinaryMorphologyImageFilter<TInputImage, TOutputImage, TKernel>
{
public:
  /** Extract dimension from input and output image. */
  itkStaticConstMacro( InputImageDimension, unsigned int,
                       TInputImage::ImageDimension );
  itkStaticConstMacro( OutputImageDimension, unsigned int,
                       TOutputImage::ImageDimension );

  /** Extract the dimension of the kernel */
  itkStaticConstMacro( KernelDimension, unsigned int,
                       TKernel::NeighborhoodDimension );
  
  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage                                 InputImageType;
  typedef TOutputImage                                OutputImageType;
  typedef typename TOutputImage::PixelType            OutputPixelType;
  typedef TKernel                                     KernelType;

  /** Standard class typedefs. */
  typedef BinaryReinhardtMorphologicalImageFilter     Self;
  typedef BinaryMorphologyImageFilter
    <InputImageType, OutputImageType, KernelType>     Superclass;
  typedef SmartPointer<Self>                          Pointer;
  typedef SmartPointer<const Self>                    ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( BinaryReinhardtMorphologicalImageFilter, 
    BinaryMorphologyImageFilter )

  /**
   * salt and pepper repair.
   */
  itkSetMacro( EmploySaltAndPepperRepair, bool );
  itkGetConstMacro( EmploySaltAndPepperRepair, bool );
  itkBooleanMacro( EmploySaltAndPepperRepair );

  itkSetMacro( SaltAndPepperMinimumSizeInPixels, unsigned int );
  itkGetConstMacro( SaltAndPepperMinimumSizeInPixels, unsigned int );

  /**
   * minimum diameter filter
   */
  itkSetMacro( EmployMinimumDiameterFilter, bool );
  itkGetConstMacro( EmployMinimumDiameterFilter, bool );
  itkBooleanMacro( EmployMinimumDiameterFilter );

  itkSetMacro( MinimumDiameterStructuringElementRadius, unsigned int );
  itkGetConstMacro( MinimumDiameterStructuringElementRadius, unsigned int );

  /**
   * unwanted cavity deletion
   */
  itkSetMacro( EmployUnwantedCavityDeletion, bool );
  itkGetConstMacro( EmployUnwantedCavityDeletion, bool );
  itkBooleanMacro( EmployUnwantedCavityDeletion );

  /**
   * minimum size filter
   */
  itkSetMacro( EmployMinimumSizeFilter, bool );
  itkGetConstMacro( EmployMinimumSizeFilter, bool );
  itkBooleanMacro( EmployMinimumSizeFilter );

  itkSetMacro( MinimumSizeStructuringElementRadius, unsigned int );
  itkGetConstMacro( MinimumSizeStructuringElementRadius, unsigned int );

  /**
   * maximum diameter filter
   */
  itkSetMacro( EmployMaximumDiameterFilter, bool );
  itkGetConstMacro( EmployMaximumDiameterFilter, bool );
  itkBooleanMacro( EmployMaximumDiameterFilter );

  itkSetMacro( MaximumDiameterStructuringElementRadius, unsigned int );
  itkGetConstMacro( MaximumDiameterStructuringElementRadius, unsigned int );

  /**
   * connectivity filter
   */
  itkSetMacro( EmployConnectivityFilter, bool );
  itkGetConstMacro( EmployConnectivityFilter, bool );
  itkBooleanMacro( EmployConnectivityFilter );

  itkSetMacro( NumberOfConnectedComponents, unsigned int );
  itkGetConstMacro( NumberOfConnectedComponents, unsigned int );

  /**
   * boundary smoother
   */
  itkSetMacro( EmployBoundarySmoother, bool );
  itkGetConstMacro( EmployBoundarySmoother, bool );
  itkBooleanMacro( EmployBoundarySmoother );

  itkSetMacro( BoundarySmootherStructuringElementRadius, unsigned int );
  itkGetConstMacro( BoundarySmootherStructuringElementRadius, unsigned int );

  /**
   * unclassified pixel processing
   */
  itkSetMacro( EmployUnclassifiedPixelProcessing, bool );
  itkGetConstMacro( EmployUnclassifiedPixelProcessing, bool );
  itkBooleanMacro( EmployUnclassifiedPixelProcessing );

protected:
  BinaryReinhardtMorphologicalImageFilter();
  virtual ~BinaryReinhardtMorphologicalImageFilter(){}
  void PrintSelf(std::ostream& os, Indent indent) const;

  void GenerateData();

  // type inherited from the superclass
  typedef typename Superclass::NeighborIndexContainer     NeighborIndexContainer;

private:
  BinaryReinhardtMorphologicalImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  void SaltAndPepperRepair( typename OutputImageType::Pointer );
  void MinimumDiameterFilter( typename OutputImageType::Pointer );
  void UnwantedCavityDeletion( typename OutputImageType::Pointer );
  void MinimumSizeFilter( typename OutputImageType::Pointer );
  void MaximumDiameterFilter( typename OutputImageType::Pointer );
  void ConnectivityFilter( typename OutputImageType::Pointer );
  void BoundarySmoother( typename OutputImageType::Pointer );
  void UnclassifiedPixelProcessing( typename OutputImageType::Pointer );

  bool                                                  m_EmploySaltAndPepperRepair;
  unsigned int                                          m_SaltAndPepperMinimumSizeInPixels;

  bool                                                  m_EmployMinimumDiameterFilter;
  unsigned int                                          m_MinimumDiameterStructuringElementRadius;
  
  bool                                                  m_EmployUnwantedCavityDeletion;

  bool                                                  m_EmployMinimumSizeFilter;
  unsigned int                                          m_MinimumSizeStructuringElementRadius;
  
  bool                                                  m_EmployMaximumDiameterFilter;
  unsigned int                                          m_MaximumDiameterStructuringElementRadius;

  bool                                                  m_EmployConnectivityFilter;
  unsigned int                                          m_NumberOfConnectedComponents;

  bool                                                  m_EmployBoundarySmoother;
  unsigned int                                          m_BoundarySmootherStructuringElementRadius;
  
  bool                                                  m_EmployUnclassifiedPixelProcessing;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBinaryReinhardtMorphologicalImageFilter.txx"
#endif

#endif

/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDiceMetricImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:52 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkDiceMetricImageFilter_h
#define __itkDiceMetricImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkNumericTraits.h"

namespace itk {

/** \class DiceMetricImageFilter
 * \brief Computes the dice metric between the set same set of labels of
 * non-zero pixels of two images.
 *
 * \sa DiceMetricImageFilter
 *
 * \ingroup MultiThreaded
 */
template<class TInputImage1, class TInputImage2>
class ITK_EXPORT DiceMetricImageFilter :
    public ImageToImageFilter<TInputImage1, TInputImage1>
{
public:
  /** Standard Self typedef */
  typedef DiceMetricImageFilter                          Self;
  typedef ImageToImageFilter<TInputImage1,TInputImage1>  Superclass;
  typedef SmartPointer<Self>                             Pointer;
  typedef SmartPointer<const Self>                       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Runtime information support. */
  itkTypeMacro( DiceMetricImageFilter, ImageToImageFilter );

  /** Image related typedefs. */
  typedef TInputImage1                              InputImage1Type;
  typedef TInputImage2                              InputImage2Type;
  typedef typename TInputImage1::Pointer            InputImage1Pointer;
  typedef typename TInputImage2::Pointer            InputImage2Pointer;
  typedef typename TInputImage1::ConstPointer       InputImage1ConstPointer;
  typedef typename TInputImage2::ConstPointer       InputImage2ConstPointer;

  typedef typename TInputImage1::RegionType         RegionType;
  typedef typename TInputImage1::SizeType           SizeType;
  typedef typename TInputImage1::IndexType          IndexType;

  typedef typename TInputImage1::PixelType          InputImage1PixelType;
  typedef typename TInputImage2::PixelType          InputImage2PixelType;

  /** Image related typedefs. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TInputImage1::ImageDimension );

  /** Type to use form computations. */
  typedef typename NumericTraits<InputImage1PixelType>::RealType RealType;

  /** Set the first input. */
  void SetInput1( const InputImage1Type * image )
  { this->SetInput( image ); }

  /** Set the second input. */
  void SetInput2( const InputImage2Type * image );

  /** Get the first input. */
  const InputImage1Type * GetInput1( void )
  { return this->GetInput(); }

  /** Get the second input. */
  const InputImage2Type * GetInput2( void );

  /** Return the computed dice metric. */
  itkGetMacro( DiceMetric, RealType );

  itkSetMacro( PixelLabel, InputImage1PixelType );
  itkGetMacro( PixelLabel, InputImage1PixelType );

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( Input1HasNumericTraitsCheck,
    ( Concept::HasNumericTraits<InputImage1PixelType> ) );
  /** End concept checking */
#endif

protected:
  DiceMetricImageFilter();
  ~DiceMetricImageFilter(){};
  void PrintSelf( std::ostream& os, Indent indent ) const;

  /** GenerateData. */
  void  GenerateData();

  // Override since the filter needs all the data for the algorithm
  void GenerateInputRequestedRegion();

  // Override since the filter produces all of its output
  void EnlargeOutputRequestedRegion( DataObject *data );

private:
  DiceMetricImageFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  RealType                                          m_DiceMetric;
  InputImage1PixelType                              m_PixelLabel;

} ; // end of class

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiceMetricImageFilter.hxx"
#endif

#endif

/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkScalarImageToRunLengthMatrixFilter.h,v $
  Language:  C++
  Date:      $Date: 2009/05/02 03:00:26 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkScalarImageToRunLengthMatrixFilter_h
#define __itkScalarImageToRunLengthMatrixFilter_h

#include "itkImage.h"
#include "itkHistogram.h"
#include "itkDenseFrequencyContainer2.h"
#include "itkMacro.h"
#include "itkNumericTraits.h"
#include "itkProcessObject.h"
#include "itkVectorContainer.h"

namespace itk {
namespace Statistics {

/** \class ScalarImageToRunLengthMatrixFilter
*
* Author: Nick Tustison
*/

template<class TImageType,
          class THistogramFrequencyContainer = DenseFrequencyContainer2>
class ScalarImageToRunLengthMatrixFilter : public ProcessObject
{
public:
  /** Standard typedefs */
  typedef ScalarImageToRunLengthMatrixFilter        Self;
  typedef ProcessObject                             Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( ScalarImageToRunLengthMatrixFilter, Object );

  /** standard New() method support */
  itkNewMacro( Self );

  typedef TImageType                                      ImageType;
  typedef typename ImageType::Pointer                     ImagePointer;
  typedef typename ImageType::ConstPointer                ImageConstPointer;
  typedef typename ImageType::PixelType                   PixelType;
  typedef typename ImageType::IndexType                   IndexType;
  typedef typename ImageType::RegionType                  RegionType;
  typedef typename ImageType::SizeType                    RadiusType;
  typedef typename ImageType::OffsetType                  OffsetType;
  typedef VectorContainer<unsigned char, OffsetType>      OffsetVector;
  typedef typename OffsetVector::Pointer                  OffsetVectorPointer;
  typedef typename ImageType::PointType                   PointType;

  typedef typename NumericTraits<PixelType>::RealType     MeasurementType;
  typedef typename NumericTraits<PixelType>::RealType     RealType;

  typedef Histogram< MeasurementType, THistogramFrequencyContainer >
                                                          HistogramType;
  typedef typename HistogramType::Pointer                 HistogramPointer;
  typedef typename HistogramType::ConstPointer            HistogramConstPointer;
  typedef typename HistogramType::MeasurementVectorType   MeasurementVectorType;

  /** ImageDimension constants */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TImageType::ImageDimension );

  itkStaticConstMacro( DefaultBinsPerAxis, unsigned int, 256 );

  /** Method to set/get the image */
  void SetInput( const ImageType* image );
  const ImageType* GetInput() const;

  /** Method to set/get the mask image */
  void SetMaskImage( const ImageType* image );
  const ImageType* GetMaskImage() const;

  /** method to get the Histogram */
  const HistogramType * GetOutput() const;

  /** Set the pixel value of the mask that should be considered "inside" the
    object. Defaults to one. */
  itkSetMacro( InsidePixelValue, PixelType );
  itkGetConstMacro( InsidePixelValue, PixelType );

  /** Set number of histogram bins along each axis */
  itkSetMacro( NumberOfBinsPerAxis, unsigned int );
  itkGetConstMacro( NumberOfBinsPerAxis, unsigned int );

  /** Set the min and max (inclusive) pixel
      value that will be placed in the histogram */
  void SetPixelValueMinMax( PixelType min, PixelType max );
  itkGetConstMacro( Min, PixelType );
  itkGetConstMacro( Max, PixelType );

  /** Set the min and max (inclusive) distance
      value that will be placed in the histogram */
  void SetDistanceValueMinMax( RealType min, RealType max );
  itkGetConstMacro( MinDistance, RealType );
  itkGetConstMacro( MaxDistance, RealType );

  /** Set the offset or offsets over which the co-occurrence pairs
   * will be computed. Calling either of these methods clears the
   * previous offsets. */
  itkSetObjectMacro( Offsets, OffsetVector );
  itkGetConstObjectMacro( Offsets, OffsetVector );
  void SetOffset( const OffsetType offset );
  void AddOffset( const OffsetType offset );

protected:
  ScalarImageToRunLengthMatrixFilter();
  virtual ~ScalarImageToRunLengthMatrixFilter() {};
  void PrintSelf( std::ostream& os, Indent indent ) const;
  virtual void FillHistogram( RadiusType radius, RegionType region );
  virtual void FillHistogramWithMask( RadiusType radius, RegionType region,
    const ImageType * maskImage );

  /** Standard itk::ProcessObject subclass method. */
  typedef DataObject::Pointer DataObjectPointer;
  virtual DataObjectPointer MakeOutput(unsigned int idx);

  /** This method causes the filter to generate its output. */
  virtual void GenerateData();

private:
  ScalarImageToRunLengthMatrixFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  OffsetVectorPointer      m_Offsets;
  PixelType                m_Min;
  PixelType                m_Max;
  RealType                 m_MinDistance;
  RealType                 m_MaxDistance;

  unsigned int             m_NumberOfBinsPerAxis;
  MeasurementVectorType    m_LowerBound;
  MeasurementVectorType    m_UpperBound;

  PixelType                m_InsidePixelValue;
};

} // end of namespace Statistics
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkScalarImageToRunLengthMatrixFilter.txx"
#endif

#endif

/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkDynamicHistogramWarpingImageFilter_h
#define __itkDynamicHistogramWarpingImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkHistogram.h"

#include <deque>

namespace itk
{
/**
 * \class DynamicHistogramWarpingImageFilter
 * \brief Normalize the grayscale values between two images by dynamic
 * histogram warping.
 *
 *
 * \author Nicholas J. Tustison
 *
 * \par REFERENCE
 *
 * Patent number: 5,727,080
 * Date of patent: Mar. 10, 1998
 *
 * Cox, I.J.; Roy, S.; Hingorani, S.L.; "Dynamic histogram warping of image
 * pairs for constant image brightness," Proceedings:  International Conference
 * on Image Processing, vol.2, pp.366-369, 1995, doi: 10.1109/ICIP.1995.537491
 *
 * \ingroup IntensityImageFilters Multithreaded
 *
 */
/* THistogramMeasurement -- The precision level for which to do
  HistogramMeasurmenets */
template<class TInputImage, class TOutputImage,
  class THistogramMeasurement = ITK_TYPENAME TInputImage::PixelType>
class ITK_EXPORT DynamicHistogramWarpingImageFilter:
  public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef DynamicHistogramWarpingImageFilter              Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage>   Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( DynamicHistogramWarpingImageFilter,
    ImageToImageFilter );

  /** ImageDimension enumeration. */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TInputImage::ImageDimension );
  itkStaticConstMacro( OutputImageDimension, unsigned int,
    TOutputImage::ImageDimension );

  /** Typedef to describe the output image region type. */
  typedef typename TOutputImage::RegionType OutputImageRegionType;

  /** Inherited typedefs. */
  typedef typename Superclass::InputImageType         InputImageType;
  typedef typename Superclass::InputImagePointer      InputImagePointer;
  typedef typename Superclass::InputImageConstPointer InputImageConstPointer;
  typedef typename Superclass::OutputImageType        OutputImageType;
  typedef typename Superclass::OutputImagePointer     OutputImagePointer;

  /** Pixel related typedefs. */
  typedef typename InputImageType::PixelType  InputPixelType;
  typedef typename OutputImageType::PixelType OutputPixelType;

  /** Histogram related typedefs. */
  typedef Statistics::Histogram<THistogramMeasurement>   HistogramType;
  typedef typename HistogramType::Pointer                HistogramPointer;
  typedef double                                         RealType;

  /** Set/Get the source image. */
  void SetSourceImage( const InputImageType *source )
    {
    this->SetInput( source );
    }
  const InputImageType * GetSourceImage()
    {
    return this->GetInput();
    }

  /** Set/Get the reference image. */
  void SetReferenceImage( const InputImageType *reference );

  const InputImageType * GetReferenceImage();

  /** Set/Get the number of histogram levels used. */
  itkSetMacro( NumberOfHistogramLevels, SizeValueType );
  itkGetConstMacro( NumberOfHistogramLevels, SizeValueType );

  /** Set/Get the size of the maximum compression level for the reference image. */
  itkSetMacro( MaximumReferenceBinCompressionSize, SizeValueType );
  itkGetConstMacro( MaximumReferenceBinCompressionSize, SizeValueType );

  /** Set/Get the size of the maximum compression level for the source image. */
  itkSetMacro( MaximumSourceBinCompressionSize, SizeValueType );
  itkGetConstMacro( MaximumSourceBinCompressionSize, SizeValueType );

  /** Set/Get the threshold at mean intensity flag.
   * If true, only source (reference) pixels which are greater
   * than the mean source (reference) intensity is used in
   * the histogram matching. If false, all pixels are
   * used. */
  itkSetMacro( ThresholdAtMeanIntensity, bool );
  itkGetConstMacro( ThresholdAtMeanIntensity, bool );
  itkBooleanMacro( ThresholdAtMeanIntensity );

  /** This filter requires all of the input to be in the buffer. */
  virtual void GenerateInputRequestedRegion();

  /** Methods to get the histograms of the source, reference, and
   * output. Objects are only valid after Update() has been called
   * on this filter. */
  itkGetObjectMacro( SourceHistogram, HistogramType );
  itkGetObjectMacro( ReferenceHistogram, HistogramType );
  itkGetObjectMacro( OutputHistogram, HistogramType );

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( IntConvertibleToInputCheck,
    ( Concept::Convertible<int, InputPixelType> ) );
  itkConceptMacro( SameDimensionCheck,
    ( Concept::SameDimension<ImageDimension, OutputImageDimension> ) );
  itkConceptMacro( DoubleConvertibleToInputCheck,
    ( Concept::Convertible<double, InputPixelType> ) );
  itkConceptMacro( DoubleConvertibleToOutputCheck,
    ( Concept::Convertible<double, OutputPixelType> ) );
  itkConceptMacro( InputConvertibleToDoubleCheck,
    ( Concept::Convertible< InputPixelType, double> ) );
  itkConceptMacro( OutputConvertibleToDoubleCheck,
    ( Concept::Convertible<OutputPixelType, double> ) );
  itkConceptMacro( SameTypeCheck,
    ( Concept::SameType<InputPixelType, OutputPixelType> ) );
  /** End concept checking */
#endif

protected:
  DynamicHistogramWarpingImageFilter();
  ~DynamicHistogramWarpingImageFilter() {}
  void PrintSelf( std::ostream & os, Indent indent ) const;

  void BeforeThreadedGenerateData();

  void AfterThreadedGenerateData();

  void ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread,
    ThreadIdType threadId );

  /**
   * Override VeriyInputInformation() since this filter does not expect
   * the input images to occupy the same physical space.
   *
			* \sa ProcessObject::VerifyInputInformation
   */
	 virtual void VerifyInputInformation() {}

  /** Compute min, max and mean of an image. */
  void ComputeMinMaxMean( const InputImageType *image,
    THistogramMeasurement & minValue, THistogramMeasurement & maxValue,
    THistogramMeasurement & meanValue );

  /** Construct a histogram from an image. */
  void ConstructHistogram( const InputImageType *image,
    HistogramType *histogram, const THistogramMeasurement minValue,
    const THistogramMeasurement maxValue );

  void ConstructCumulativeHistogram( const HistogramType *histogram,
    HistogramType *cumulativeHistogram );

  void CalculateOptimalHistogramMapping( const HistogramType *referenceHistogram,
    const HistogramType *sourceHistogram );

private:
  DynamicHistogramWarpingImageFilter( const Self & ); //purposely not implemented
  void operator=( const Self & );               //purposely not implemented

  SizeValueType m_NumberOfHistogramLevels;
  SizeValueType m_MaximumSourceBinCompressionSize;
  SizeValueType m_MaximumReferenceBinCompressionSize;

  bool            m_ThresholdAtMeanIntensity;
  InputPixelType  m_SourceIntensityThreshold;
  InputPixelType  m_ReferenceIntensityThreshold;
  OutputPixelType m_OutputIntensityThreshold;

  THistogramMeasurement m_SourceMinValue;
  THistogramMeasurement m_SourceMaxValue;
  THistogramMeasurement m_SourceMeanValue;

  THistogramMeasurement m_ReferenceMinValue;
  THistogramMeasurement m_ReferenceMaxValue;
  THistogramMeasurement m_ReferenceMeanValue;

  THistogramMeasurement m_OutputMinValue;
  THistogramMeasurement m_OutputMaxValue;
  THistogramMeasurement m_OutputMeanValue;

  HistogramPointer m_SourceHistogram;
  HistogramPointer m_ReferenceHistogram;
  HistogramPointer m_OutputHistogram;

  std::deque<RealType> m_WarpedSourceHistogramProfile;
  std::deque<RealType> m_WarpedReferenceHistogramProfile;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDynamicHistogramWarpingImageFilter.hxx"
#endif

#endif

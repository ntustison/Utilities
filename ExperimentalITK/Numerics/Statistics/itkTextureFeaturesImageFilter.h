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
#ifndef __itkTextureFeaturesImageFilter_h
#define __itkTextureFeaturesImageFilter_h

#include "itkConstNeighborhoodIterator.h"
#include "itkDenseFrequencyContainer2.h"
#include "itkHistogram.h"
#include "itkImageToImageFilter.h"
#include "itkVectorContainer.h"
#include "itkVectorImage.h"

namespace itk
{
namespace Statistics
{
/** \class TextureFeaturesImageFilter
 *  \brief This filter computes texture features based on Haralick's
 * cooccurrence matrix.
 *
 * \author
 * \ingroup
 */
template<class TInputImage, class TOutputImage
  = VectorImage<typename TInputImage::PixelType, TInputImage::ImageDimension> >
class ITK_EXPORT TextureFeaturesImageFilter:
  public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef TextureFeaturesImageFilter                      Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage>   Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;

  /** Standard New method. */
  itkNewMacro( Self );

  /** Runtime information support. */
  itkTypeMacro( TextureFeaturesImageFilter, ImageToImageFilter );

  /** ImageDimension constants */
  itkStaticConstMacro( ImageDimension, unsigned int, TInputImage::ImageDimension );

  /** Some convenient typedefs. */
  typedef float                                      RealType;
  typedef TInputImage                                InputImageType;
  typedef typename InputImageType::RegionType        RegionType;
  typedef typename InputImageType::OffsetType        OffsetType;
  typedef std::vector<OffsetType>                    OffsetVectorType;
  typedef typename OffsetType::OffsetValueType       OffsetValueType;
  typedef Image<unsigned int, ImageDimension>        MaskImageType;
  typedef TOutputImage                               OutputImageType;
  typedef typename InputImageType::PixelType         InputPixelType;
  typedef typename OutputImageType::PixelType        OutputPixelType;

  typedef typename NumericTraits<InputPixelType>::RealType                         MeasurementType;
  typedef Statistics::DenseFrequencyContainer2                                     HistogramFrequencyContainerType;
  typedef Statistics::Histogram<MeasurementType, HistogramFrequencyContainerType>  HistogramType;
  typedef typename HistogramType::MeasurementVectorType                            MeasurementVectorType;

  /** Set/Get the input mask image that will constraint the computation to
   * pixels that are on in the mask. This is intended to reduce the computation time.
   */
  void SetMaskImage( const MaskImageType *mask );

  const MaskImageType * GetMaskImage() const;

  /** Type of the neighborhood iterator used to evaluate similarity between the
   * image pixels. */
  typedef ConstNeighborhoodIterator< InputImageType >        ConstNeighborhoodIteratorType;
  typedef typename ConstNeighborhoodIteratorType::RadiusType RadiusType;

  /** Radius defining the local window for evaluating the texture features. */
  itkSetMacro( NeighborhoodRadius, RadiusType );
  itkGetConstMacro( NeighborhoodRadius, RadiusType);

  /** Set the offset or offsets over which the co-occurrence pairs will be computed.
      Calling either of these methods clears the previous offsets. */
  itkGetConstReferenceMacro( Offsets, OffsetVectorType );
  void SetOffsets( const OffsetVectorType &offsets )
  {
    if ( this->m_Offsets != offsets )
      {
      this->m_Offsets.assign( offsets.begin(), offsets.end() );
      this->Modified();
      }
  }
  void SetOffset( const OffsetType &offset );

  /** Set number of histogram bins along each axis */
  itkSetMacro( NumberOfBinsPerAxis, unsigned int );
  itkGetConstMacro( NumberOfBinsPerAxis, unsigned int );

  /**
   * Set the min and max (inclusive) pixel value that will be placed in the
   * histogram.
   */
  void SetPixelValueMinMax( InputPixelType, InputPixelType );

  itkGetConstMacro( Min, InputPixelType );
  itkGetConstMacro( Max, InputPixelType );

  /**
   * Set the calculator to normalize the histogram (divide all bins by the
   * total frequency). Normalization is off by default.
   */
  itkSetMacro( Normalize, bool );
  itkGetConstMacro( Normalize, bool );
  itkBooleanMacro( Normalize );

  itkSetMacro( InsidePixelValue, typename MaskImageType::PixelType );
  itkGetConstMacro( InsidePixelValue, typename MaskImageType::PixelType );

  unsigned int GetNumberOfOutputComponents() { return 8; }

protected:
  TextureFeaturesImageFilter();
  ~TextureFeaturesImageFilter() {};

  virtual void GenerateInputRequestedRegion()
  {
    // currently we require the entire input image to process
    TInputImage *input = const_cast<TInputImage *>( this->GetInput() );
    input->SetRequestedRegionToLargestPossibleRegion();
  }

  virtual void ThreadedGenerateData( const RegionType &, ThreadIdType );

  virtual void BeforeThreadedGenerateData();

  void PrintSelf( std::ostream & os, Indent indent ) const;

  void GenerateOutputInformation();

private:
  TextureFeaturesImageFilter( const Self & );          //purposely not implemented
  void operator=( const Self & );                      //purposely not implemented

  OffsetVectorType                                  m_Offsets;
  std::vector<std::pair<OffsetType, OffsetType> >   m_CooccurenceOffsetVector;

  RadiusType                                        m_NeighborhoodRadius;
  InputPixelType                                    m_Min;
  InputPixelType                                    m_Max;
  unsigned int                                      m_NumberOfBinsPerAxis;
  bool                                              m_Normalize;
  typename MaskImageType::PixelType                 m_InsidePixelValue;

}; // end of class
} // end namespace statistics
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTextureFeaturesImageFilter.hxx"
#endif

#endif

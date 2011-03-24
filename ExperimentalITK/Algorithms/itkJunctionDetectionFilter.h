/*=========================================================================
 * Automatic Junction Detection for Tubular Structures
 * Author:      Guanglei Xiong (guangleixiong at gmail.com)
 * Date:        March 23, 2009
 * Reference:   http://www.insight-journal.org/browse/publication/324
=========================================================================*/

#ifndef __itkJunctionDetectionFilter_h
#define __itkJunctionDetectionFilter_h

#include "itkImageToImageFilter.h"
#include <map>

namespace itk
{

/** \class JunctionDetectionFilter
 * 
 * This class is designed to be used for detection of junctions 
 * in tubular structures. The input is a binary image from the output 
 * of a segmentation method. The output is a label image that have 
 * junctions labeled with numbered spheres of their sizes.
 *
 * The idea of the algorithm is very straightforward. For every pixel
 * inside of a tubular structure, a hollow sphere (defined by InnerRadius 
 * and OuterRadius coefficient times the distance to the surface of 
 * the tubular structure) is placed. We test how many disconnected 
 * components are within the hollow sphere. Clearly, the number of 
 * components need to be more than 3 if the center of the hollow sphere 
 * is located at the junction. Please note that the distance measure 
 * we are using is not Euclidean distance but geodesic distance within 
 * the structure. The method works in 2D, 3D, or even nD.
 *
 * - todo
 *   - optimize for speed. The use of SetLocation of neighborhood 
 *   iterators seems to be slow. Other implementation may be possible. 
 *   Moreover, It seems that testing for every pixel is not necessary 
 *   if neighboring relation between pixels are considered.
 *   - use parallelism/threading
 *
 */
template <class TInputImage>
class ITK_EXPORT JunctionDetectionFilter : 
  public ImageToImageFilter< TInputImage,
    Image<long, ::itk::GetImageDimension<TInputImage>::ImageDimension> >
{
public:
  /** Dimension of the input and output images. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** The type of input image. */
  typedef TInputImage InputImageType;

  /** The type of output image. */
  typedef Image<long, ImageDimension> OutputImageType;

  /** Standard class typedefs. */
  typedef JunctionDetectionFilter                               Self;
  typedef ImageToImageFilter<InputImageType,OutputImageType>    Superclass;
  typedef SmartPointer<Self>                                    Pointer;
  typedef SmartPointer<const Self>                              ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(JunctionDetectionFilter, ImagetoImageFilter);

  /** Image typedef support. */
  typedef typename InputImageType::PixelType          InputPixelType;
  typedef typename InputImageType::IndexType          IndexType;
  typedef typename InputImageType::OffsetType         OffsetType;
  typedef typename InputImageType::SpacingType        SpacingType;

  /** Junction label map typedef support. */
  typedef std::pair< IndexType, float >               JCLabelPairType;
  typedef std::map< long, JCLabelPairType >           JCLabelMapType;

  /** Set/Get the background pixel value of the input binary image. */
  itkSetMacro(BackgroundValue, InputPixelType);
  itkGetConstMacro(BackgroundValue, InputPixelType);

  /** Set/Get the inner radius coefficient. The actual inner radius 
   * is this value multiplied by the distance to the wall. */
  itkSetMacro(InnerRadius, float);
  itkGetConstMacro(InnerRadius, float);

  /** Set/Get the outer radius coefficient. The actual inner radius
   * is this value multiplied by the distance to the wall. */
  itkSetMacro(OuterRadius, float);
  itkGetConstMacro(OuterRadius, float);

  /** Set/Get the minimal number of pixels inside hallow sphere. */
  itkSetMacro(MinNumberOfPixel, unsigned long);
  itkGetConstMacro(MinNumberOfPixel, unsigned long);

  /** Get the junction label map. */
  //itkGetConstReferenceMacro(JCLabelMap, JCLabelMapType);
  const JCLabelMapType & GetJCLabelMap() const { return this->m_JCLabelMap; }

protected:
  JunctionDetectionFilter();
  virtual ~JunctionDetectionFilter();
  void PrintSelf(std::ostream& os, Indent indent) const;
  void GenerateData();

private:
  JunctionDetectionFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** The background pixel value of the input binary image. */
  InputPixelType m_BackgroundValue;

  /** The inner radius coefficient. */
  double m_InnerRadius;

  /** The outer radius coefficient. */
  double m_OuterRadius;

  /** The minimal number of pixels inside hallow sphere. */
  unsigned long m_MinNumberOfPixel;
  
  /** The junction label map. */
  JCLabelMapType m_JCLabelMap;

}; // end of class

} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkJunctionDetectionFilter.txx"
#endif

#endif

#ifndef __itkBranchDecompositionFilter_h
#define __itkBranchDecompositionFilter_h

#include "itkImageToImageFilter.h"
#include <map>
#include <set>

namespace itk
{

/** \class BranchDecompositionFilter
 * 
 * This class is designed to be used for decompsition of tubular 
 * structures into branches. The input is a binary image from the
 * output of a segmentation method and junction information stored 
 * in a map. The output is a label image that have subdomains labeled 
 * with their unique identifiers.
 *
 * The idea of the algorithm is very straightforward.
 * The method works in 2D, 3D, or even nD.
 *
 * - todo
 *   - optimize for speed. The use of SetLocation of neighborhood 
 *   iterators seems to be slow. Other implementation may be possible. 
 *   Moreover, It seems that testing for every pixel is not necessary 
 *   if neighboring relation between pixels are considered.
 *   - use parallelism/threading.
 *   - further division of the junction by assigning its pixels to
 *   neighboring branches. 
 *
 * email: guangleixiong at gmail.com
 *
 */
template <class TInputImage>
class ITK_EXPORT BranchDecompositionFilter : 
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
  typedef BranchDecompositionFilter                             Self;
  typedef ImageToImageFilter<InputImageType,OutputImageType>    Superclass;
  typedef SmartPointer<Self>                                    Pointer;
  typedef SmartPointer<const Self>                              ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(BranchDecompositionFilter, ImagetoImageFilter);

  /** Image typedef support. */
  typedef typename InputImageType::PixelType          InputPixelType;
  typedef typename InputImageType::IndexType          IndexType;
  typedef typename InputImageType::OffsetType         OffsetType;
  typedef typename InputImageType::SpacingType        SpacingType;

  /** Junction label map typedef support. */
  typedef std::pair< IndexType, float >               JCLabelPairType;
  typedef std::map< long, JCLabelPairType >           JCLabelMapType;
  
  /** Branch junction connection map typedef support. */
  typedef std::set< long >                            BRJCSetType;
  typedef std::map< long, BRJCSetType >               BRJCConnectMapType;

  /** Set/Get the inner radius coefficient. The actual inner radius 
   * is this value multiplied by the junction radius. */
  itkSetMacro(InnerRadius, float);
  itkGetConstMacro(InnerRadius, float);

  /** Set/Get the outer radius coefficient. The actual inner radius
   * is this value multiplied by the junction radius. */
  itkSetMacro(OuterRadius, float);
  itkGetConstMacro(OuterRadius, float);

  /** Set/Get the minimal number of pixels to be considered as a subdomain. */
  itkSetMacro(MinNumberOfPixel, unsigned long);
  itkGetConstMacro(MinNumberOfPixel, unsigned long);

  /** Set the junction label map. */
  itkSetMacro(JCLabelMap, JCLabelMapType);

  /** Get the branch junctions connection map. */
  //itkGetConstReferenceMacro(BRJCConnectMap, BRJCConnectMapType);
  const BRJCConnectMapType & GetBRJCConnectMap() const { return this->m_BRJCConnectMap; }

protected:
  BranchDecompositionFilter();
  virtual ~BranchDecompositionFilter();
  void PrintSelf(std::ostream& os, Indent indent) const;
  void GenerateData();

private:
  BranchDecompositionFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** The inner radius coefficient. */
  double m_InnerRadius;

  /** The outer radius coefficient. */
  double m_OuterRadius;

  /** The minimal number of pixels to be considered as a subdomain. */
  unsigned long m_MinNumberOfPixel;
  
  /** The junction label map. */
  JCLabelMapType* m_JCLabelMap;

  /** The branch junction connection map. */
  BRJCConnectMapType m_BRJCConnectMap;

}; // end of class

} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBranchDecompositionFilter.hxx"
#endif

#endif

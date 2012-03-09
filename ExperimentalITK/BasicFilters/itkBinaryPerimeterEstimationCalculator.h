#ifndef __itkBinaryPerimeterEstimationCalculator_h
#define __itkBinaryPerimeterEstimationCalculator_h

#include "itkObject.h"

namespace itk {

/** \class BinaryPerimeterEstimationCalculator
 * \brief TODO
 *
 * \author Gaëtan Lehmann. Biologie du Développement et de la Reproduction, INRA de Jouy-en-Josas, France.
 *
 * \sa 
 */
template<class TInputImage>
class ITK_EXPORT BinaryPerimeterEstimationCalculator : 
    public Object
{
public:
  /** Standard class typedefs. */
  typedef BinaryPerimeterEstimationCalculator Self;
  typedef Object Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage InputImageType;
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename InputImageType::RegionType     InputImageRegionType;
  typedef typename InputImageType::PixelType      InputImagePixelType;
  
  typedef typename InputImageType::RegionType     RegionType;
  typedef typename InputImageType::SizeType       SizeType;
  typedef typename InputImageType::IndexType      IndexType;
  
  /** ImageDimension constants */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Standard New method. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(BinaryPerimeterEstimationCalculator, Object);
  
  /**
   * Set/Get whether the connected components are defined strictly by
   * face connectivity or by face+edge+vertex connectivity.  Default is
   * FullyConnectedOff.  For objects that are 1 pixel wide, use
   * FullyConnectedOn.
   */
  itkSetMacro(FullyConnected, bool);
  itkGetConstReferenceMacro(FullyConnected, bool);
  itkBooleanMacro(FullyConnected);

  /**
   * Set/Get the value used as "foreground" in the output image.
   * Defaults to NumericTraits<PixelType>::max().
   */
  itkSetMacro(ForegroundValue, InputImagePixelType);
  itkGetConstMacro(ForegroundValue, InputImagePixelType);

  itkGetConstMacro(Perimeter, double);

  itkSetObjectMacro( Image, InputImageType );
//  void SetImage( const InputImageType * img )
//    {
//    m_Image = img;
//    }
//  itkGetObjectMacro( Image, InputImageType );
  InputImageType * GetImage()
    {
    return m_Image;
    }
    
  void Compute();
  

protected:
  BinaryPerimeterEstimationCalculator();
  ~BinaryPerimeterEstimationCalculator() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  BinaryPerimeterEstimationCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  bool m_FullyConnected;

  InputImagePixelType m_ForegroundValue;

  double m_Perimeter;
  
  InputImageType * m_Image;

} ; // end of class

} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBinaryPerimeterEstimationCalculator.txx"
#endif

#endif



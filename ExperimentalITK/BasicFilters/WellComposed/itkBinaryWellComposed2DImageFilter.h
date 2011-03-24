#ifndef __itkBinaryWellComposed2DImageFilter_h
#define __itkBinaryWellComposed2DImageFilter_h

#include "itkArray.h"
#include "itkInPlaceImageFilter.h"
#include "itkNeighborhoodIterator.h"

#include <list>

/** \class BinaryWellComposed2DImageFilter
 *  This filter transforms an arbitrary 2D binary image to be well-composed.
 *
 *  \par Inputs and Outputs
 *  This is an in-place-image filter.  The input is a 2D binary image.  
 *  The output is a well-composed version of the input with the border
 *  pixel values assigned to the background value.
 *
 *  \par Parameters
 *  The user simply needs to specify the 
 *  background and foreground values (default = 0 and 1, respectively). 
 *
 *  \cite
 * Longin Jan Latecki,
 * "Discrete representation of spatial objects in computer vision",
 * Kluwer Academic Publishers, 1998, pp. 66-67.
 *
 * \ingroup ImageFeatureExtraction 
 */

namespace itk
{
template <class TImage>
class BinaryWellComposed2DImageFilter 
: public InPlaceImageFilter<TImage>
{
public:

  /** Extract dimension from input image. */
  itkStaticConstMacro( ImageDimension, 
         unsigned int, TImage::ImageDimension);

  /** Convenient typedefs for simplifying declarations. */
  typedef TImage                                ImageType;

  /** Standard class typedefs. */
  typedef BinaryWellComposed2DImageFilter       Self;
  typedef InPlaceImageFilter<ImageType>         Superclass;
  typedef SmartPointer<Self>                    Pointer;
  typedef SmartPointer<const Self>              ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Image typedef support. */
  typedef typename ImageType::PixelType         PixelType;
  typedef typename ImageType::IndexType         IndexType;
  typedef NeighborhoodIterator<ImageType>       NeighborhoodIteratorType;

  /** 
   * Set the background value which defines the object.  Default
   * value is = 0.
   */
  itkSetMacro( BackgroundValue, PixelType );
  itkGetConstMacro( BackgroundValue, PixelType );

  /** 
   * Set the background value which defines the object.  Default
   * value is = 1.
   */
  itkSetMacro( ForegroundValue, PixelType );
  itkGetConstMacro( ForegroundValue, PixelType );

  /** 
   * Set whether full invariance is required.
   */
  itkBooleanMacro( FullInvariance );
  itkSetMacro( FullInvariance, bool );
  itkGetConstMacro( FullInvariance, bool );


protected:
  BinaryWellComposed2DImageFilter();
  virtual ~BinaryWellComposed2DImageFilter();
  void PrintSelf( std::ostream& os, Indent indent ) const;
  void GenerateData();

private:
  BinaryWellComposed2DImageFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented
   
  void InitializeIndices();
  void MakeImageWellComposed();

  void ConvertBoundaryForegroundPixelsToBackgroundPixels();

  bool IsCriticalC1Configuration( Array<PixelType> );
  bool IsCriticalC2Configuration( Array<PixelType> );
  bool IsCriticalC3Configuration( Array<PixelType> );
  bool IsCriticalC4Configuration( Array<PixelType> );
  bool IsSpecialCaseOfC4Configuration( IndexType, IndexType, IndexType );


  bool                                          m_FullInvariance;
  PixelType                                     m_BackgroundValue;  
  PixelType                                     m_ForegroundValue;
  unsigned long                                 m_NumberOfC1Configurations;
  unsigned long                                 m_NumberOfC2Configurations; 
  unsigned long                                 m_NumberOfC3Configurations;
  unsigned long                                 m_NumberOfC4Configurations; 
  Array<unsigned int>                           m_RotationIndices[4];
  Array<unsigned int>                           m_ReflectionIndices[2];
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBinaryWellComposed2DImageFilter.txx"
#endif

#endif

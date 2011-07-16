#ifndef __itkBinaryWellComposed3DImageFilter_h
#define __itkBinaryWellComposed3DImageFilter_h

#include "itkArray.h"
#include "itkInPlaceImageFilter.h"
#include "itkNeighborhoodIterator.h"

#include <list>

/** \class BinaryWellComposed3DImageFilter
 *  This filter transforms an arbitrary 3D binary image to be well-composed.
 *
 *  \par Inputs and Outputs
 *  This is an in-place-image filter.  The input is a 3D binary image.  
 *  The output is a well-composed version of the input with the border
 *  pixel values assigned to the background value.
 *
 *  \par Parameters
 *  The user simply needs to specify the 
 *  background and foreground values (default = 0 and 1, respectively). 
 *
 *  \cite
 * Marcelo Siqueira and Longin Jan Latecki and Jean Gallier,
 * "Making 3D Binary Digital Images Well-Composed",
 * Proceedings of the IS&T/SPIE Conference on Vision Geometry XIII, 2005.
 *
 * \ingroup ImageFeatureExtraction 
 */

namespace itk
{
template <class TImage>
class BinaryWellComposed3DImageFilter 
: public InPlaceImageFilter<TImage>
{
public:

  /** Extract dimension from input image. */
  itkStaticConstMacro(ImageDimension, unsigned int, TImage::ImageDimension);

  /** Convenient typedefs for simplifying declarations. */
  typedef TImage                                ImageType;

  /** Standard class typedefs. */
  typedef BinaryWellComposed3DImageFilter       Self;
  typedef InPlaceImageFilter<ImageType>         Superclass;
  typedef SmartPointer<Self>                    Pointer;
  typedef SmartPointer<const Self>              ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Image typedef support. */
  typedef typename ImageType::PixelType         PixelType;
  typedef typename ImageType::IndexType         IndexType;
  typedef std::list<IndexType>                  IndexContainerType;
  typedef NeighborhoodIterator<ImageType>       NeighborhoodIteratorType;
  typedef typename NeighborhoodIteratorType::
                     OffsetType                 OffsetType;

  /** 
   * Set the background value which defines the object.  Default
   * value is = 0.
   */
  itkSetMacro( BackgroundValue, PixelType );
  itkGetConstReferenceMacro( BackgroundValue, PixelType );

  /** 
   * Set the background value which defines the object.  Default
   * value is = 1.
   */
  itkSetMacro( ForegroundValue, PixelType );
  itkGetConstReferenceMacro( ForegroundValue, PixelType );

protected:
  BinaryWellComposed3DImageFilter();
  virtual ~BinaryWellComposed3DImageFilter();
  void PrintSelf( std::ostream& os, Indent indent ) const;
  void GenerateData();

private:
  BinaryWellComposed3DImageFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented
   
  void InitializeOffsetsAndIndices();
  void ConvertBoundaryForegroundPixelsToBackgroundPixels();
  void CountCriticalConfigurations();
  void LocateCriticalConfigurations();
  void MakeImageWellComposed();

  void InsertCriticalConfiguration( IndexType );  
  bool IsCriticalC1Configuration( Array<unsigned char> );
  bool IsCriticalC2Configuration( Array<unsigned char> );
  bool IsChangeSafe( IndexType );
  void MakeRandomChange( unsigned char &, IndexType, 
                         unsigned char &, IndexType );
  void MakeRandomChange( unsigned char &, IndexType, 
                         unsigned char &, IndexType, 
                         unsigned char &, IndexType, 
                         unsigned char &, IndexType, 
                         unsigned char &, IndexType, 
                         unsigned char &, IndexType );

  PixelType              m_BackgroundValue;  
  PixelType              m_ForegroundValue;
  unsigned long          m_NumberOfC1Configurations;
  unsigned long          m_NumberOfC2Configurations; 
  OffsetType             m_Offsets8[8];
  OffsetType             m_Offsets27[27];
  Array<unsigned int>    m_C1IndicesI[6];
  Array<unsigned int>    m_C2IndicesI[6];
  Array<unsigned int>    m_C1IndicesII[12];
  Array<unsigned int>    m_C2IndicesII[8];
  IndexContainerType     m_CriticalConfigurationIndices;
  IndexContainerType     m_NewCriticalConfigurationIndices;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBinaryWellComposed3DImageFilter.hxx"
#endif

#endif

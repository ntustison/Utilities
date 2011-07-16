#ifndef __itkBinaryWellComposedImageFilter_h
#define __itkBinaryWellComposedImageFilter_h

#include "itkArray.h"
#include "itkInPlaceImageFilter.h"
#include "itkNeighborhoodIterator.h"

#include <list>

/** \class BinaryWellComposedImageFilter
 *  This filter transforms an arbitrary  binary image to be well-composed.
 *
 *  \par Inputs and Outputs
 *  This is an in-place-image filter.  The input is a  binary image.  
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
 * Longin Jan Latecki,
 * "Discrete representation of spatial objects in computer vision",
 * Kluwer Academic Publishers, 1998, pp. 66-67.
 *
 * \ingroup ImageFeatureExtraction 
 */

namespace itk
{
template <class TImage>
class BinaryWellComposedImageFilter 
: public InPlaceImageFilter<TImage>
{
public:

  /** Extract dimension from input image. */
  itkStaticConstMacro( ImageDimension, 
         unsigned int, TImage::ImageDimension);

  /** Convenient typedefs for simplifying declarations. */
  typedef TImage                                ImageType;

  /** Standard class typedefs. */
  typedef BinaryWellComposedImageFilter         Self;
  typedef InPlaceImageFilter<ImageType>         Superclass;
  typedef SmartPointer<Self>                    Pointer;
  typedef SmartPointer<const Self>              ConstPointer;
  
  /** Run-time type information (and related methods) */
  itkTypeMacro( BinaryWellComposedImageFilter, InPlaceImageFilter );
  
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
  itkGetConstMacro( BackgroundValue, PixelType );

  /** 
   * Set the background value which defines the object.  Default
   * value is = 1.
   */
  itkSetMacro( ForegroundValue, PixelType );
  itkGetConstMacro( ForegroundValue, PixelType );

  /** 
   * Set whether full invariance is required (2-D only).
   */
  itkBooleanMacro( FullInvariance );
  itkSetMacro( FullInvariance, bool );
  itkGetConstMacro( FullInvariance, bool );


protected:
  BinaryWellComposedImageFilter();
  virtual ~BinaryWellComposedImageFilter();
  void PrintSelf( std::ostream& os, Indent indent ) const;
  void GenerateData();

private:
  BinaryWellComposedImageFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented
   
  void ConvertBoundaryForegroundPixelsToBackgroundPixels();

  /**
   * Functions/data for the 2-D case
   */
  void Make2DImageWellComposed();
  void InitializeIndices();
  bool IsCriticalC1Configuration( Array<PixelType> );
  bool IsCriticalC2Configuration( Array<PixelType> );
  bool IsCriticalC3Configuration( Array<PixelType> );
  bool IsCriticalC4Configuration( Array<PixelType> );
  bool IsSpecialCaseOfC4Configuration( IndexType, IndexType, IndexType );

  bool                   m_FullInvariance;
  PixelType              m_BackgroundValue;  
  PixelType              m_ForegroundValue;
  unsigned long          m_NumberOfC1Configurations;
  unsigned long          m_NumberOfC2Configurations; 
  unsigned long          m_NumberOfC3Configurations;
  unsigned long          m_NumberOfC4Configurations; 
  Array<unsigned int>    m_RotationIndices[4];
  Array<unsigned int>    m_ReflectionIndices[2];

  /**
   * Functions/data for the 3-D case
   */
  void Make3DImageWellComposed();
  void InitializeOffsetsAndIndices();
  void CountCriticalConfigurations();
  void LocateCriticalConfigurations();
  void InsertCriticalConfiguration( IndexType );  
  bool IsCriticalC1Configuration3D( Array<unsigned char> );
  bool IsCriticalC2Configuration3D( Array<unsigned char> );
  bool IsChangeSafe( IndexType );
  void MakeRandomChange( unsigned char &, IndexType, 
                         unsigned char &, IndexType );
  void MakeRandomChange( unsigned char &, IndexType, 
                         unsigned char &, IndexType, 
                         unsigned char &, IndexType, 
                         unsigned char &, IndexType, 
                         unsigned char &, IndexType, 
                         unsigned char &, IndexType );
                         
  unsigned long          m_NumberOfC1Configurations3D;
  unsigned long          m_NumberOfC2Configurations3D; 
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
#include "itkBinaryWellComposedImageFilter.hxx"
#endif

#endif

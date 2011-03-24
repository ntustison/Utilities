#ifndef __itkWellComposedImageFilter_h
#define __itkWellComposedImageFilter_h

#include "itkArray.h"
#include "itkInPlaceImageFilter.h"
#include "itkNeighborhoodIterator.h"

#include <deque>

/** \class WellComposedImageFilter
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
class WellComposedImageFilter 
: public InPlaceImageFilter<TImage>
{
public:

  /** Extract dimension from input image. */
  itkStaticConstMacro( ImageDimension, 
         unsigned int, TImage::ImageDimension);

  /** Convenient typedefs for simplifying declarations. */
  typedef TImage                                ImageType;

  /** Standard class typedefs. */
  typedef WellComposedImageFilter         Self;
  typedef InPlaceImageFilter<ImageType>         Superclass;
  typedef SmartPointer<Self>                    Pointer;
  typedef SmartPointer<const Self>              ConstPointer;
  
  /** Run-time type information (and related methods) */
  itkTypeMacro( WellComposedImageFilter, InPlaceImageFilter );
  
  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Image typedef support. */
  typedef typename ImageType::PixelType         PixelType;
  typedef typename ImageType::IndexType         IndexType;
  typedef std::deque<IndexType>                 IndexContainerType;
  typedef NeighborhoodIterator<ImageType>       NeighborhoodIteratorType;

  /** 
   * Set the total number of labels.  
   */
  itkSetMacro( TotalNumberOfLabels, unsigned long );
  itkGetConstMacro( TotalNumberOfLabels, unsigned long );

  /** 
   * Set whether full invariance is required (2-D only).
   */
  itkBooleanMacro( FullInvariance );
  itkSetMacro( FullInvariance, bool );
  itkGetConstMacro( FullInvariance, bool );

protected:
  WellComposedImageFilter();
  virtual ~WellComposedImageFilter();
  void PrintSelf( std::ostream& os, Indent indent ) const;
  void GenerateData();

private:
  WellComposedImageFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented
   
  void ConvertBoundaryPixels( PixelType );
  unsigned int                                  m_TotalNumberOfLabels;
  IndexContainerType                            m_CriticalConfigurationIndices;
  IndexContainerType                            m_NewCriticalConfigurationIndices;

  /**
   * Functions/data for the 2-D case
   */
  void MakeImageWellComposed2D();
  void InitializeIndices2D();
  bool IsChangeSafe2D( PixelType, IndexType );
  bool IsCriticalC1Configuration2D( Array<char> );
  bool IsCriticalC2Configuration2D( Array<char> );
  bool IsCriticalC3Configuration2D( Array<char> );
  bool IsCriticalC4Configuration2D( Array<char> );
  bool IsSpecialCaseOfC4Configuration2D( PixelType, IndexType, 
                                         IndexType, IndexType );

  bool                            m_FullInvariance;
  Array<unsigned int>             m_RotationIndices[4];
  Array<unsigned int>             m_ReflectionIndices[2];

  /**
   * Functions/data for the 3-D case
   */
  void MakeImageWellComposed3D();
  void InitializeIndices3D();
  void LocateCriticalConfigurations3D( PixelType );
  void InsertCriticalConfiguration3D( PixelType, IndexType );  
  bool IsCriticalC1Configuration3D( Array<char> );
  unsigned int IsCriticalC2Configuration3D( Array<char> );
  bool IsChangeSafe3D( PixelType, IndexType );
  bool RemoveCriticalC1Configuration3D( int, PixelType, IndexType );
  bool RemoveCriticalC2Configuration3D( int, PixelType, IndexType );

  Array<unsigned int>             m_C1Indices[12];
  Array<unsigned int>             m_C2Indices[8];
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkWellComposedImageFilter.txx"
#endif

#endif

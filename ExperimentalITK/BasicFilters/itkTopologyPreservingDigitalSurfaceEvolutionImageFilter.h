#ifndef __itkTopologyPreservingDigitalSurfaceEvolutionImageFilter_h
#define __itkTopologyPreservingDigitalSurfaceEvolutionImageFilter_h

#include "itkArray.h"
#include "itkImageToImageFilter.h"
#include "itkNeighborhoodIterator.h"

#include <vector>

/** \class TopologyPreservingDigitalSurfaceEvolutionImageFilter
 *
 */

namespace itk
{
template <class TImage>
class TopologyPreservingDigitalSurfaceEvolutionImageFilter
: public ImageToImageFilter<TImage, TImage>
{
public:

  /** Extract dimension from input image. */
  itkStaticConstMacro( ImageDimension,
         unsigned int, TImage::ImageDimension );


  /** Standard class typedefs. */
  typedef TopologyPreservingDigitalSurfaceEvolutionImageFilter    Self;
  typedef ImageToImageFilter<TImage, TImage>                      Superclass;
  typedef SmartPointer<Self>                                      Pointer;
  typedef SmartPointer<const Self>                                ConstPointer;

  /** Run-time type information (and related methods) */
  itkTypeMacro( TopologyPreservingDigitalSurfaceEvolutionImageFilter,
    ImageToImageFilter );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Image typedef support. */
  typedef TImage                                ImageType;
  typedef typename ImageType::PixelType         PixelType;
  typedef typename ImageType::IndexType         IndexType;
  typedef std::vector<IndexType>                IndexContainerType;
  typedef std::vector<PixelType>                LabelArrayType;
  typedef NeighborhoodIterator<ImageType>       NeighborhoodIteratorType;

  typedef float                                 RealType;
  typedef Image<RealType,
    itkGetStaticConstMacro( ImageDimension )>   RealImageType;

  void SetSourceImage( const ImageType *source )
    {
    itkDebugMacro( "setting input SourceImage to " << source );
    if( source != static_cast<ImageType *>(
      this->ProcessObject::GetInput( 0 ) ) )
      {
      this->ProcessObject::SetNthInput( 0, const_cast<ImageType *>( source ) );
      this->Modified();
      }
    }
  const ImageType* GetSourceImage()
    {
    itkDebugMacro( "returning input SourceImage of " <<
      static_cast<const ImageType *>( this->ProcessObject::GetInput() ) );
    return static_cast<const ImageType *>( this->ProcessObject::GetInput() );
    }

  itkSetObjectMacro( TargetImage, RealImageType );
  itkGetObjectMacro( TargetImage, RealImageType );

  itkSetMacro( NumberOfIterations, unsigned int );
  itkGetConstMacro( NumberOfIterations, unsigned int );

  itkSetClampMacro( ThresholdValue, RealType, 0.0, 1.0 );
  itkGetConstMacro( ThresholdValue, RealType );

  itkSetClampMacro( GluingStrategy, unsigned int, 0, 2 );
  itkGetConstMacro( GluingStrategy, unsigned int );

  itkSetMacro( GrowOnly, bool );
  itkGetConstMacro( GrowOnly, bool );
  itkBooleanMacro( GrowOnly );

  itkSetMacro( GlueObjects, bool );
  itkGetConstMacro( GlueObjects, bool );
  itkBooleanMacro( GlueObjects );

  itkSetMacro( FillCavities, bool );
  itkGetConstMacro( FillCavities, bool );
  itkBooleanMacro( FillCavities );

  LabelArrayType GetOrderedObjectLabels()
    {
    return this->m_OrderedObjectLabels;
    }
  LabelArrayType GetCavityLabels()
    {
    return this->m_CavityLabels;
    }
  LabelArrayType GetGluingLabels()
    {
    return this->m_GluingLabels;
    }

protected:
  TopologyPreservingDigitalSurfaceEvolutionImageFilter();
  virtual ~TopologyPreservingDigitalSurfaceEvolutionImageFilter();
  void PrintSelf( std::ostream& os, Indent indent ) const;
  void GenerateData();

private:
  //purposely not implemented
  TopologyPreservingDigitalSurfaceEvolutionImageFilter( const Self& );
  void operator=( const Self& );

  typename RealImageType::Pointer            m_TargetImage;
  typename ImageType::Pointer                m_LabelSurfaceImage;

  /**
   * user-selected parameters
   */
  unsigned int                               m_NumberOfIterations;
  RealType                                   m_ThresholdValue;
  bool                                       m_GrowOnly;
  bool                                       m_GlueObjects;
  bool                                       m_FillCavities;
  unsigned int                               m_GluingStrategy;

  /**
   * label arrays created and used internally.  Perhaps of use to the user
   * following execution.
   */
  LabelArrayType                             m_OrderedObjectLabels;
  LabelArrayType                             m_GluingLabels;
  LabelArrayType                             m_CavityLabels;

  /**
   * variables for internal use
   */
  bool                                       m_IsInversionStep;
  PixelType                                  m_SurfaceLabel;
  PixelType                                  m_ForegroundValue;
  PixelType                                  m_BackgroundValue;

  /**
   * Functions/data for the 2-D case
   */
  void InitializeIndices2D();
  bool IsChangeWellComposed2D( IndexType );
  bool IsCriticalC1Configuration2D( Array<short> );
  bool IsCriticalC2Configuration2D( Array<short> );
  bool IsCriticalC3Configuration2D( Array<short> );
  bool IsCriticalC4Configuration2D( Array<short> );
  bool IsSpecialCaseOfC4Configuration2D( PixelType, IndexType,
                                         IndexType, IndexType );

  Array<unsigned int>                        m_RotationIndices[4];
  Array<unsigned int>                        m_ReflectionIndices[2];

  /**
   * Functions/data for the 3-D case
   */
  void InitializeIndices3D();
  bool IsCriticalC1Configuration3D( Array<short> );
  unsigned int IsCriticalC2Configuration3D( Array<short> );
  bool IsChangeWellComposed3D( IndexType );

  bool IsChangeSafe( IndexType );
  bool IsCriticalTopologicalConfiguration( IndexType );

  Array<unsigned int>                        m_C1Indices[12];
  Array<unsigned int>                        m_C2Indices[8];


  /**
   * Other internal functions
   */
  void CreateLabelSurfaceImage();
  void FindLabelOrderingForObjectGluing();

  void GlueObjects();
  void FillCavities();
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTopologyPreservingDigitalSurfaceEvolutionImageFilter.hxx"
#endif

#endif

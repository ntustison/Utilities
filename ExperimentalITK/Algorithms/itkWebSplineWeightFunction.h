#ifndef _itkWebSplineWeightFunction_h_
#define _itkWebSplineWeightFunction_h_

#include "itkImageFunction.h"
#include "itkSmoothRepresentationFunction.h"

namespace itk {

template <class TInputImage, class TOutput = double, class TCoordRep = double>
class ITK_EXPORT WebSplineWeightFunction 
: public ImageFunction<TInputImage, TOutput, TCoordRep>
{
public:
  typedef WebSplineWeightFunction                            Self;
  typedef ImageFunction<TInputImage, TOutput, TCoordRep>     Superclass;
  typedef SmartPointer<Self>                                 Pointer;
  typedef SmartPointer<const Self>                           ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  
  
  /** Run-time type information (and related methods) */
  itkTypeMacro(WebSplineWeightFunction, ImageFunction);
  
  typedef typename Superclass::IndexType                     IndexType;
  typedef typename Superclass::ContinuousIndexType           ContinuousIndexType;     
  typedef typename Superclass::PointType                     PointType;

  typedef TInputImage                                        ImageType;
  typedef typename ImageType::PixelType                      PixelType;
  typedef TOutput                                            OutputType;
  typedef TCoordRep                                          CoordRepType;
  
  /** Dimensionality of input and output data is assumed to be the same. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       ImageType::ImageDimension );
                       
  typedef double                                             RealType;
  typedef Image<RealType, 
     itkGetStaticConstMacro( ImageDimension )>               RealImageType;
  typedef SmoothRepresentationFunction<RealImageType, 
                     OutputType, CoordRepType>               InterpolateFunctionType;                     

  itkSetMacro( BackgroundValue, PixelType );
  itkGetConstMacro( BackgroundValue, PixelType );

  itkSetMacro( Delta, RealType );
  itkGetConstMacro( Delta, RealType );

  itkSetMacro( Gamma, RealType );
  itkGetConstMacro( Gamma, RealType );

  /** Set the input image.
   * \warning this method caches BufferedRegion information.
   * If the BufferedRegion has changed, user must call
   * SetInputImage again to update cached values. */
  virtual void SetInputImage( const ImageType * );

  /** Evaluate the function at specified Point position. */
  virtual TOutput Evaluate( const PointType & ) const;

  /** Evaluate the function at specified Index position. */
  virtual TOutput EvaluateAtIndex( const IndexType &idx ) const
    {
    PointType pt;
    this->m_Image->TransformIndexToPhysicalPoint( idx, pt );
    return this->Evaluate( pt );
    }

  /** Evaluate the function at specified ContinousIndex position. */
  virtual TOutput EvaluateAtContinuousIndex( const ContinuousIndexType &cidx ) const
    {
    PointType pt;
    this->m_Image->TransformContinuousIndexToPhysicalPoint( cidx, pt );
    return this->Evaluate( pt );
    }

  /** Evaluate the gradient at specified Point position. */
  PointType EvaluateGradient( const PointType & ) const;

  /** Evaluate the gradient at specified Index position. */
  PointType EvaluateGradientAtIndex( const IndexType &idx ) const
    {
    PointType pt;
    this->m_Image->TransformIndexToPhysicalPoint( idx, pt );
    return this->EvaluateGradient( pt );
    }

  /** Evaluate the gradient at specified ContinousIndex position. */
  PointType EvaluateGradientAtContinuousIndex( const ContinuousIndexType &cidx ) const
    {
    PointType pt;
    this->m_Image->TransformContinuousIndexToPhysicalPoint( cidx, pt );
    return this->EvaluateGradient( pt );
    }
  

protected :
  /** de/constructor */
  WebSplineWeightFunction(); 
  virtual ~WebSplineWeightFunction(); 

  void PrintSelf( std::ostream& os, Indent indent ) const;
  void GenerateData();

private : 
  
  WebSplineWeightFunction( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented
  
  typename InterpolateFunctionType::Pointer   m_Interpolator; 
  PixelType                                   m_BackgroundValue;
  RealType                                    m_Delta;
  RealType                                    m_Gamma;

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkWebSplineWeightFunction.hxx"
#endif

#endif



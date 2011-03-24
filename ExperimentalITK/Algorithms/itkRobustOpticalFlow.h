/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRobustOpticalFlow.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:13:43 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkRobustOpticalFlow_h_
#define _itkRobustOpticalFlow_h_

#include "itkAvantsPDEDeformableRegistrationFunction.h"
#include "itkPoint.h"
#include "itkCovariantVector.h"
#include "itkInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkCentralDifferenceImageFunction.h"
#include "itkGradientRecursiveGaussianImageFilter.h"

namespace itk {

/**
 * \class RobustOpticalFlow
 *
 * This class encapsulate the PDE which drives the demons registration 
 * algorithm. It is used by CrossCorrelationRegistrationFilter to compute the 
 * output deformation field which will map a moving image onto a
 * a fixed image.
 *
 * Non-integer moving image values are obtained by using
 * interpolation. The default interpolator is of type
 * LinearInterpolateImageFunction. The user may set other
 * interpolators via method SetMovingImageInterpolator. Note that the input
 * interpolator must derive from baseclass InterpolateImageFunction.
 *
 * This class is templated over the fixed image type, moving image type,
 * and the deformation field type.
 *
 * \warning This filter assumes that the fixed image type, moving image type
 * and deformation field type all have the same number of dimensions.
 *
 * \sa CrossCorrelationRegistrationFilter
 * \ingroup FiniteDifferenceFunctions
 */
template<class TFixedImage, class TMovingImage, class TDeformationField>
class ITK_EXPORT RobustOpticalFlow : 
  public AvantsPDEDeformableRegistrationFunction< TFixedImage,
    TMovingImage, TDeformationField>
{
public:
  /** Standard class typedefs. */
  typedef RobustOpticalFlow    Self;
  typedef AvantsPDEDeformableRegistrationFunction< TFixedImage,
    TMovingImage, TDeformationField >    Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( RobustOpticalFlow, 
    AvantsPDEDeformableRegistrationFunction );

  /** MovingImage image type. */
  typedef typename Superclass::MovingImageType     MovingImageType;
  typedef typename Superclass::MovingImagePointer  MovingImagePointer;

  /** FixedImage image type. */
  typedef typename Superclass::MetricImageType     MetricImageType;
  typedef typename Superclass::MetricImageType::Pointer     MetricImagePointer;
  typedef typename Superclass::FixedImageType     FixedImageType;
  typedef typename Superclass::FixedImagePointer  FixedImagePointer;
  typedef typename FixedImageType::IndexType      IndexType;
  typedef typename FixedImageType::SizeType       SizeType;
  
  /** Deformation field type. */
  typedef typename Superclass::DeformationFieldType    DeformationFieldType;
  typedef typename Superclass::DeformationFieldTypePointer   
    DeformationFieldTypePointer;
  typedef typename TDeformationField::PixelType VectorType;
  
  typedef CovariantVector<float,
          itkGetStaticConstMacro(ImageDimension)> GradientPixelType;
  typedef Image<GradientPixelType,
               itkGetStaticConstMacro(ImageDimension)> GradientImageType;
  typedef SmartPointer<GradientImageType>     GradientImagePointer;
  typedef GradientRecursiveGaussianImageFilter< MetricImageType,GradientImageType > 
          GradientImageFilterType;   
  typedef typename GradientImageFilterType::Pointer GradientImageFilterPointer;
  typedef Image<float,itkGetStaticConstMacro(ImageDimension)> BinaryImageType;
  typedef typename BinaryImageType::Pointer BinaryImagePointer;

  /** Inherit some enums from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

  /** Inherit some enums from the superclass. */
  typedef typename Superclass::PixelType     PixelType;
  typedef typename Superclass::RadiusType    RadiusType;
  typedef typename Superclass::NeighborhoodType    NeighborhoodType;
//  typedef typename Superclass::NeighborhoodType    BoundaryNeighborhoodType;
  typedef typename Superclass::FloatOffsetType  FloatOffsetType;
  typedef typename Superclass::TimeStepType TimeStepType;

  /** Interpolator type. */
  typedef double CoordRepType;
  typedef InterpolateImageFunction<MovingImageType,CoordRepType> InterpolatorType;
  typedef typename InterpolatorType::Pointer         InterpolatorPointer;
  typedef typename InterpolatorType::PointType       PointType;
  typedef LinearInterpolateImageFunction<MovingImageType,CoordRepType>
    DefaultInterpolatorType;

  /** Covariant vector type. */
  typedef CovariantVector<double,itkGetStaticConstMacro(ImageDimension)> CovariantVectorType;

  /** Gradient calculator type. */
  typedef CentralDifferenceImageFunction<FixedImageType> GradientCalculatorType;
  typedef typename GradientCalculatorType::Pointer   GradientCalculatorPointer;

  /** Set the moving image interpolator. */
  void SetMovingImageInterpolator( InterpolatorType * ptr )
    { m_MovingImageInterpolator = ptr; }

  /** Get the moving image interpolator. */
  InterpolatorType * GetMovingImageInterpolator(void)
    { return m_MovingImageInterpolator; }

  double ComputeMetricAtPair(IndexType fixedindex, typename TDeformationField::PixelType vec);

  /** This class uses a constant timestep of 1. */
  virtual TimeStepType ComputeGlobalTimeStep(void *GlobalData) const
    { return m_TimeStep; }

  /** Return a pointer to a global data structure that is passed to
   * this object from the solver at each calculation.  */
  virtual void *GetGlobalDataPointer() const
    {
    GlobalDataStruct *global = new GlobalDataStruct();
    return global;
    }

  /** Release memory for global data structure. */
  virtual void ReleaseGlobalDataPointer( void *GlobalData ) const
    { delete (GlobalDataStruct *) GlobalData;  }

  /** Set the object's state before each iteration. */
  virtual void InitializeIteration();


  virtual VectorType  OpticalFlowUpdate(const NeighborhoodType &neighborhood)
  {
  // Get fixed image related information
  IndexType index=neighborhood.GetIndex();    
		typename TDeformationField::PixelType vec;
		if ( Superclass::m_DeformationField )
				{
				vec = Superclass::m_DeformationField->GetPixel(index);
				}
		else
				{
				vec.Fill( 0 );
				} 
  VectorType update;
  update.Fill(0.0);
  double fixedValue;
  CovariantVectorType fixedGradient;
  double fixedGradientSquaredMagnitude = 0;
  fixedValue = (double) Superclass::m_FixedImage->GetPixel( index );
  fixedGradient = m_FixedImageGradientCalculator->EvaluateAtIndex( index );
  for( unsigned int j = 0; j < ImageDimension; j++ )
    {
    fixedGradientSquaredMagnitude += vnl_math_sqr( fixedGradient[j] );
    } 
  double movingValue;
  int j;
  PointType mappedPoint;
  for( j = 0; j < ImageDimension; j++ )
    {
    mappedPoint[j] = double( index[j] ) * m_FixedImageSpacing[j] + 
      m_FixedImageOrigin[j];
    mappedPoint[j] += vec[j];
    }
  if( m_MovingImageInterpolator->IsInsideBuffer( mappedPoint ) )
    {
    movingValue = m_MovingImageInterpolator->Evaluate( mappedPoint );
    }
  else
    {
    for( j = 0; j < ImageDimension; j++ )
      {
      update[j] = 0.0;
      }
    return update;
    }
  double speedValue = fixedValue - movingValue;
  double denominator = vnl_math_sqr( speedValue ) / m_Normalizer + 
    fixedGradientSquaredMagnitude;
  double m_DenominatorThreshold = 1e-9;
  double m_IntensityDifferenceThreshold = 0.001;  
  if ( vnl_math_abs(speedValue) < m_IntensityDifferenceThreshold || 
       denominator < m_DenominatorThreshold )
    {
    for( j = 0; j < ImageDimension; j++ )
      {
      update[j] = 0.0;
      }
    return update;
    }
  for( j = 0; j < ImageDimension; j++ )
    {
    update[j] = speedValue * fixedGradient[j] / denominator;
    }
    return update;
  }

  virtual VectorType  ComputeUpdate2(const NeighborhoodType &neighborhood,
                     void *globalData,
                     const FloatOffsetType &offset = FloatOffsetType(0.0));
					 

  virtual VectorType ComputeUpdate(const NeighborhoodType &neighborhood,
                                   void *globalData,
                                   const FloatOffsetType &offset = FloatOffsetType(0.0))
  {
    bool m_Use1SidedDiff=false;
    VectorType update;
    update.Fill(0.0);   
    IndexType oindex = neighborhood.GetIndex();

    FixedImageType* img =const_cast<FixedImageType *>(Superclass::m_FixedImage.GetPointer()); 
    if (!img) return update;
    typename FixedImageType::SpacingType spacing=img->GetSpacing();
    typename FixedImageType::SizeType imagesize=img->GetLargestPossibleRegion().GetSize();
    bool inimage=true;
    
    for (unsigned int dd=0; dd<ImageDimension; dd++)
    {
      if ( oindex[dd] <= 3 || 
           oindex[dd] >= static_cast<typename IndexType::IndexValueType>(imagesize[dd]-3) ) 
           return update;
    }    

    VectorType fixedGradient;
    fixedGradient=this->OpticalFlowUpdate(neighborhood);		
    typename TDeformationField::PixelType vec = Superclass::m_DeformationField->GetPixel(oindex);
    float loce=0.0;
    //    for (int imd=0; imd<ImageDimension; imd++)
    {
      double nccp1=0,nccm1=0;
      typename TDeformationField::PixelType fdvec1 = Superclass::m_DeformationField->GetPixel(oindex);
      typename TDeformationField::PixelType fdvec2 = Superclass::m_DeformationField->GetPixel(oindex);
      //      float step=spacing[imd]*0.5;
      for (int imd=0; imd<ImageDimension; imd++)
	{
	  fdvec1[imd]=fixedGradient[imd];
	  fdvec2[imd]=fixedGradient[imd]*(-1.);
	}
      if (binaryimage->GetPixel(oindex) > m_Thresh)
	{
	  nccp1=this->ComputeMetricAtPair(oindex,fdvec1);
	  nccm1=this->ComputeMetricAtPair(oindex,fdvec2);
	}
      float sign=nccp1-nccm1; 
      if (sign < 0 ) sign=(0.0);  
      if (sign > 0.0 ) sign=1.0;
      for (int imd=0; imd<ImageDimension; imd++) update[imd]=sign*fixedGradient[imd];
      loce+=(nccp1+nccm1);
    }
	
	
    loce/=(2.0*(float)ImageDimension);//this->ComputeMetricAtPair(oindex,vec);
     Superclass:: Superclass::m_Energy+=loce;
    float mag=0;
    for (int imd=0; imd<ImageDimension; imd++) mag+=update[imd]*update[imd];
    mag=sqrt(mag);
    if (mag >  1.) 
    {
      update/=mag;
      mag=1;
    }
    m_AvgMag+=mag;
    if (mag>m_MaxMag) m_MaxMag=mag;
    if (mag<m_MinMag) m_MinMag=mag;
    binaryimage->SetPixel(oindex,mag);
    return update* Superclass:: Superclass::m_GradientStep;
  }

protected:
  RobustOpticalFlow();
  ~RobustOpticalFlow() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** FixedImage image neighborhood iterator type. */
  typedef ConstNeighborhoodIterator<FixedImageType> FixedImageNeighborhoodIteratorType;

  /** A global data type for this class of equation. Used to store
   * iterators for the fixed image. */
  struct GlobalDataStruct
   {
   FixedImageNeighborhoodIteratorType   m_FixedImageIterator;
   };

  MetricImagePointer MakeImage() 
  {
    typedef ImageRegionIteratorWithIndex<MetricImageType> ittype;
    typedef ImageRegionIteratorWithIndex<BinaryImageType> ittype2;
    FixedImageType* img =const_cast<FixedImageType *>(Superclass::m_FixedImage.GetPointer()); 
    typename FixedImageType::SizeType imagesize=img->GetLargestPossibleRegion().GetSize();
  
    {
  //    ittype it(m_MetricImage,m_MetricImage->GetLargestPossibleRegion().GetSize());
  //    for( it.GoToBegin(); !it.IsAtEnd(); ++it ) it.Set(0);
    }
    bool makebinimg=false;
    if (!binaryimage) makebinimg=true;
    else if ( binaryimage->GetLargestPossibleRegion().GetSize()[0] != 
                      img->GetLargestPossibleRegion().GetSize()[0] )makebinimg=true;

    if (makebinimg)
    { 
      m_Iteration=0;
      binaryimage = BinaryImageType::New();
      binaryimage->SetLargestPossibleRegion(img->GetLargestPossibleRegion()  );
      binaryimage->SetBufferedRegion(img->GetLargestPossibleRegion());
      binaryimage->SetSpacing(img->GetSpacing());
      binaryimage->SetOrigin(img->GetOrigin());
      binaryimage->Allocate();   
      ittype2 it(binaryimage,binaryimage->GetLargestPossibleRegion().GetSize());
      for( it.GoToBegin(); !it.IsAtEnd(); ++it ) it.Set(1);
    }
    return binaryimage;
  }


private:
  RobustOpticalFlow(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  /** Cache fixed image information. */
  typename TFixedImage::SpacingType                m_FixedImageSpacing;
  typename TFixedImage::PointType                  m_FixedImageOrigin;

  /** Function to compute derivatives of the fixed image. */
  GradientCalculatorPointer       m_FixedImageGradientCalculator;

  GradientCalculatorPointer       m_MetricImageGradientCalculator;

  /** Function to interpolate the moving image. */
  InterpolatorPointer             m_MovingImageInterpolator;

  /** The global timestep. */
  TimeStepType                    m_TimeStep;

  /** Threshold below which the denominator term is considered zero. */
  double                          m_DenominatorThreshold;

  /** Threshold below which two intensity value are assumed to match. */
  double                           m_IntensityDifferenceThreshold;

  mutable double                   m_MetricTotal;
  mutable float                    m_MinMag;
  mutable float                    m_MaxMag;
  mutable float                    m_AvgMag;
  mutable float                    m_Thresh;
  
  GradientImagePointer             m_MetricGradientImage;


  MetricImagePointer               finitediffimages[ImageDimension*2];
  BinaryImagePointer               binaryimage;

  unsigned int   m_Iteration;
  float m_Normalizer;
};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRobustOpticalFlow.cxx"
#endif

#endif

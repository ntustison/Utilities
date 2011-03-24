/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRobustOpticalFlow.cxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:13:43 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkRobustOpticalFlow_txx_
#define _itkRobustOpticalFlow_txx_

#include "itkRobustOpticalFlow.h"
#include "itkExceptionObject.h"
#include "vnl/vnl_math.h"
#include "itkImageFileWriter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkMeanImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkImageFileWriter.h"
namespace itk {

/*
 * Default constructor
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
RobustOpticalFlow<TFixedImage,TMovingImage,TDeformationField>
::RobustOpticalFlow()
{

  m_Iteration=0;
  RadiusType r;
  unsigned int j;
  for( j = 0; j < ImageDimension; j++ )
    {
    r[j] = 2;
    }
  this->SetRadius(r);
   Superclass:: Superclass::m_Energy = 0.0;
  m_TimeStep = 1.0;
  m_DenominatorThreshold = 1e-9;
  m_IntensityDifferenceThreshold = 0.001;
  Superclass::m_MovingImage = NULL;
  Superclass::m_FixedImage = NULL;
  m_FixedImageSpacing.Fill( 1.0 );
  m_FixedImageOrigin.Fill( 0.0 );
  m_FixedImageGradientCalculator = GradientCalculatorType::New();
  binaryimage=NULL;

  m_MetricImageGradientCalculator = GradientCalculatorType::New();

  typename DefaultInterpolatorType::Pointer interp =
    DefaultInterpolatorType::New();

  m_MovingImageInterpolator = static_cast<InterpolatorType*>(
    interp.GetPointer() );



}


/*
 * Standard "PrintSelf" method.
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
RobustOpticalFlow<TFixedImage,TMovingImage,TDeformationField>
::PrintSelf(std::ostream& os, Indent indent) const
{
  
  Superclass::PrintSelf(os, indent);
/*
  os << indent << "MovingImageIterpolator: ";
  os << m_MovingImageInterpolator.GetPointer() << std::endl;
  os << indent << "FixedImageGradientCalculator: ";
  os << Superclass::m_FixedImageGradientCalculator.GetPointer() << std::endl;
  os << indent << "DenominatorThreshold: ";
  os << m_DenominatorThreshold << std::endl;
  os << indent << "IntensityDifferenceThreshold: ";
  os << m_IntensityDifferenceThreshold << std::endl;
*/
}


/*
 * Set the function state values before each iteration
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
RobustOpticalFlow<TFixedImage,TMovingImage,TDeformationField>
::InitializeIteration()
{
  typedef ImageRegionIteratorWithIndex<MetricImageType> ittype;
  if( !this->Superclass::m_MovingImage)
    {
      itkExceptionMacro( << "MovingImage not set" ); 
      throw ExceptionObject(__FILE__,__LINE__);
    }
  if( !this->Superclass::m_FixedImage)
    {
      itkExceptionMacro( << "FixedImage not set" ); 
      throw ExceptionObject(__FILE__,__LINE__);
    }
  if(  !this->m_MovingImageInterpolator)
    {
      itkExceptionMacro( << "MovingImageInterp not set" ); 
      throw ExceptionObject(__FILE__,__LINE__);
    }

  // cache fixed image information
  m_FixedImageSpacing    = Superclass::m_FixedImage->GetSpacing();
  m_FixedImageOrigin     = Superclass::m_FixedImage->GetOrigin();

  // setup gradient calculator
  m_FixedImageGradientCalculator->SetInputImage( Superclass::m_FixedImage);
  
  // setup moving image interpolator
  m_MovingImageInterpolator->SetInputImage( Superclass::m_MovingImage);


  unsigned long numpix=1;
  for (int i=0; i<ImageDimension;i++) numpix*=Superclass::m_DeformationField->GetLargestPossibleRegion().GetSize()[i];
//    Superclass::m_DeformationField->GetLargestPossibleRegion().GetSize()<< " image size " <<
//    Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize() << std::endl;
  m_MetricTotal=0.0;
    Superclass::m_Energy=0.0;
 
 
  // compute the normalizer
  m_Normalizer      = 0.0;
  for( unsigned int k = 0; k < ImageDimension; k++ )
    {
    m_Normalizer += m_FixedImageSpacing[k] * m_FixedImageSpacing[k];
    }
  m_Normalizer /= static_cast<double>( ImageDimension );

  m_AvgMag/=(float)(numpix+1);
  m_Thresh=m_AvgMag*1.e-6;
  if (m_Iteration==0) {m_AvgMag=0.0;m_Thresh=0.0;}
  typename FixedImageType::SpacingType spacing=this->GetFixedImage()->GetSpacing();
//  for (int i=0; i<ImageDimension*2; i++) finitediffimages[i]=this->MakeImage();
  finitediffimages[0]=this->MakeImage();
  finitediffimages[1]=this->MakeImage();
  
  typedef itk::DiscreteGaussianImageFilter<BinaryImageType, BinaryImageType> dgf1;
//  typedef itk::DiscreteGaussianImageFilter<BinaryImageType, BinaryImageType> dgf;
  typedef itk::MeanImageFilter<BinaryImageType, BinaryImageType> dgf;
  typedef itk::MedianImageFilter<BinaryImageType, BinaryImageType> dgf2;
  float sig=15.;
  
  RadiusType r;
  for( int j = 0; j < ImageDimension; j++ )    r[j] = 3;

//  if( m_Iteration <= 1 || m_Iteration % 4 == 0)
  {
      typename dgf::Pointer filter = dgf::New();
//      filter->SetVariance(sig);
//      filter->SetUseImageSpacingOn();
      filter->SetRadius(this->GetRadius());
  //    filter->SetRadius(r);
      filter->SetInput(this->GetFixedImage());
      filter->Update();
      finitediffimages[0]=filter->GetOutput(); 
      typedef ImageFileWriter<BinaryImageType> writertype;
      typename writertype::Pointer w= writertype::New();
      w->SetInput( filter->GetOutput() );
      w->SetFileName("sigf.hdr");
      //  w->Write();
  }  
  //if( m_Iteration <= 1 || m_Iteration % 4 == 0)
  {
      typename dgf::Pointer filter = dgf::New();
//      filter->SetVariance(sig);
//      filter->SetUseImageSpacingOn();
      filter->SetRadius(this->GetRadius());
    //  filter->SetRadius(r);
      filter->SetInput(this->GetMovingImage());
      filter->Update();
      finitediffimages[1]=filter->GetOutput(); 
      typedef ImageFileWriter<BinaryImageType> writertype;
      typename writertype::Pointer w= writertype::New();
      w->SetInput( filter->GetOutput() );
      w->SetFileName("sigm.hdr");
      //w->Write();
  }

  {
      typename dgf1::Pointer filter = dgf1::New();
      filter->SetVariance(1.0);
      filter->SetUseImageSpacingOn();
      filter->SetMaximumError(.01f);
      filter->SetInput(binaryimage);
      filter->Update();
      binaryimage=filter->GetOutput(); 
  }

  m_MaxMag=0.0;
  m_MinMag=9.e9;
  m_AvgMag=0.0;  
  m_Iteration++;
}


/*
 * Compute the ncc metric everywhere
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
double
RobustOpticalFlow<TFixedImage,TMovingImage,TDeformationField>
::ComputeMetricAtPair(IndexType fixedindex, typename TDeformationField::PixelType vec) 
{
  double localCrossCorrelation=1;

  bool m_UseNormalizedCrossCorrelation=false;
  typedef ImageRegionIteratorWithIndex<MetricImageType> ittype;
  unsigned int j;
  typename FixedImageType::SizeType hradius=this->GetRadius(); 
  typename FixedImageType::SpacingType spacing=this->GetMovingImage()->GetSpacing();
  FixedImageType* img =const_cast<FixedImageType *>(Superclass::m_FixedImage.GetPointer()); 
  typename FixedImageType::SizeType imagesize=img->GetLargestPossibleRegion().GetSize();
  bool inimage=true;
  
    double sff=0.0;
    double smm=0.0;
    double sfm=0.0;
    double fixedValue;
    double movingValue;

  unsigned long ct = 0;
  NeighborhoodIterator<FixedImageType>    hoodIt( hradius , img, img->GetLargestPossibleRegion());
  IndexType oindex = fixedindex;
  hoodIt.SetLocation(oindex);
  sff=0.0;
  smm=0.0;
  sfm=0.0;

  double fixedMean=0;
  double movingMean=0;

  PointType mappedPoint;
  unsigned int indct;
  unsigned int hoodlen=hoodIt.Size();
/*
    unsigned int inct=0;
    for(indct=0; indct<hoodlen; indct++)
    {

      IndexType index=hoodIt.GetIndex(indct);
      inimage=true;
      for (unsigned int dd=0; dd<ImageDimension; dd++)
      {
        if ( index[dd] < 0 || index[dd] > static_cast<typename IndexType::IndexValueType>(imagesize[dd]-1) ) inimage=false;
      }
    
      if (inimage)
      {
        fixedValue = (double) Superclass::m_FixedImage->GetPixel( index );
        for( j = 0; j < ImageDimension; j++ )
        {
          mappedPoint[j] = double( index[j] ) * m_FixedImageSpacing[j] + 
          m_FixedImageOrigin[j];
          mappedPoint[j] += vec[j];
        }
        if( m_MovingImageInterpolator->IsInsideBuffer( mappedPoint ) )
        {
          movingValue = m_MovingImageInterpolator->Evaluate( mappedPoint );
          movingMean+=movingValue;
          fixedMean+=fixedValue;
          inct++;
        }
        else
        {
          movingValue = 0.0;
        }
      }
    }
    
    if (inct > 0)fixedMean/=(float)inct;
    if (inct > 0)movingMean/=(float)inct;
    
*/  

    indct=0;
    IndexType index=oindex;
    inimage=true;
    for (unsigned int dd=0; dd<ImageDimension; dd++)
    {
      if ( index[dd] < 0 || index[dd] > static_cast<typename IndexType::IndexValueType>(imagesize[dd]-1) ) inimage=false;
    }
    if (inimage )
    {
        fixedMean = (double) finitediffimages[0]->GetPixel(oindex);
        for( j = 0; j < ImageDimension; j++ )
        {
          mappedPoint[j] = double( index[j] ) * m_FixedImageSpacing[j] + m_FixedImageOrigin[j];
          mappedPoint[j] += vec[j];
          index[j] = (unsigned long)(mappedPoint[j]/spacing[j]+0.5);
        }
        if( m_MovingImageInterpolator->IsInsideBuffer( mappedPoint ) )
        {
          movingMean=(double) finitediffimages[1]->GetPixel(index);
        }
    }

 
  
        
    for(indct=0; indct<hoodlen; indct++)
    {

      IndexType index=hoodIt.GetIndex(indct);
      bool inimage=true;
      float dist=0.0;
      for (unsigned int dd=0; dd<ImageDimension; dd++)
      {
        if ( index[dd] < 0 || index[dd] > static_cast<typename IndexType::IndexValueType>(imagesize[dd]-1) ) inimage=false;
        float temp=(float)oindex[dd] -  (float)index[dd];
        dist+=temp*temp;
      }
      
      if (inimage )
      {
        fixedValue = (double) Superclass::m_FixedImage->GetPixel( index )-fixedMean;
        for( j = 0; j < ImageDimension; j++ )
        {
          mappedPoint[j] = double( index[j] ) * m_FixedImageSpacing[j] + 
          m_FixedImageOrigin[j];
          mappedPoint[j] += vec[j];
        }
        if( m_MovingImageInterpolator->IsInsideBuffer( mappedPoint ) )
        {
          movingValue = m_MovingImageInterpolator->Evaluate( mappedPoint )-movingMean;
        }
        else
        {
          movingValue = movingMean;
        }

        sff+=fixedValue*fixedValue;
        smm+=movingValue*movingValue;
        sfm+=fixedValue*movingValue;
      }
    
    
    if(m_UseNormalizedCrossCorrelation)
    {
      if (sff*smm !=0.0) localCrossCorrelation = sfm / sqrt( sff * smm );
      else localCrossCorrelation = 1.0;
    }
    else
    {
      if (sff*smm !=0.0) localCrossCorrelation = sfm*sfm / ( sff * smm );
      else if (sff == 0.0 && smm == 0) localCrossCorrelation = 1.0;
      else localCrossCorrelation = 1.0;
    }
  }
        
  return localCrossCorrelation;

}



template <class TFixedImage, class TMovingImage, class TDeformationField>
typename RobustOpticalFlow<TFixedImage,TMovingImage,TDeformationField>
::VectorType
RobustOpticalFlow<TFixedImage,TMovingImage,TDeformationField>
::ComputeUpdate2(const NeighborhoodType &it, void * itkNotUsed(globalData),
                const FloatOffsetType& itkNotUsed(offset)) 
{

  PixelType update;
  update.Fill(0.0);
  unsigned int j;

  double fixedMean=0;
  double movingMean=0;

  IndexType oindex = it.GetIndex();

  if (false)
  {
    typename TDeformationField::PixelType vec = Superclass::m_DeformationField->GetPixel(oindex);
    float locncc=0.0;
    double nccp1,nccm1;
//    nccp1=this->ComputeMetricAtPair(oindex,vec);
    for (int imd=0; imd<ImageDimension; imd++)
    {
      typename TDeformationField::PixelType fdvec1 = Superclass::m_DeformationField->GetPixel(oindex);
      typename TDeformationField::PixelType fdvec2 = Superclass::m_DeformationField->GetPixel(oindex);
      float step=1.0;
      fdvec1[imd]=vec[imd]+step;
      nccp1=this->ComputeMetricAtPair(oindex,fdvec1);
      fdvec2[imd]=vec[imd]-step;
      nccm1=this->ComputeMetricAtPair(oindex,fdvec2);
      update[imd]=nccp1-nccm1;
      locncc+=(nccp1+nccm1);
    }
      CovariantVectorType fixedGradient;
      fixedGradient = m_FixedImageGradientCalculator->EvaluateAtIndex( oindex ); 
      PixelType update2=update;

    float mag1=0,mag2=0;
    for (int imd=0; imd<ImageDimension; imd++)
    { 
      mag1+=update[imd]*fixedGradient[imd];
      mag2+=fixedGradient[imd]*fixedGradient[imd];            
    }     
    for (int imd=0; imd<ImageDimension; imd++)
    {
      update2[imd]=mag1 / mag2 *fixedGradient[imd];
      if ( update2[imd] < 0.0 && update[imd] > 0.0) update2[imd]*=(-1.0);
      if ( update2[imd] > 0.0 && update[imd] < 0.0) update2[imd]*=(-1.0);
    }



    locncc/=(2.0*(float)ImageDimension);//this->ComputeMetricAtPair(oindex,vec);
      Superclass::m_Energy+=locncc;
    m_MetricTotal+=locncc;
    return update*  Superclass::m_GradientStep;
  }

  RadiusType r;
 
  for( j = 0; j < ImageDimension; j++ )
    {
    r[j] = 1;
    }
  typename FixedImageType::SizeType hradius=r; //it.GetRadius();

 
  FixedImageType* img =const_cast<FixedImageType *>(Superclass::m_FixedImage.GetPointer()); 

  typename FixedImageType::SizeType imagesize=img->GetLargestPossibleRegion().GetSize();
  
//  if (img->GetPixel(oindex) <= 0) return update;

  NeighborhoodIterator<FixedImageType> 
    hoodIt( hradius , img, img->GetRequestedRegion());
  hoodIt.SetLocation(oindex);



  double sff=0.0;
  double smm=0.0;
  double sfm=0.0;
  double fixedValue;
  double movingValue;

  double derivativeF[ImageDimension];
  double derivativeM[ImageDimension];
  for (j=0; j<ImageDimension;j++){
    derivativeF[j]=0;
    derivativeM[j]=0;
  }

  PointType mappedPoint;

  unsigned int indct;
  unsigned int hoodlen=hoodIt.Size();
  
  bool m_UseNormalizedCrossCorrelation=false;
  if (!m_UseNormalizedCrossCorrelation)
  {
      IndexType index=oindex;
        typename TDeformationField::PixelType vec = Superclass::m_DeformationField->GetPixel(index);   
        fixedMean = (double) finitediffimages[0]->GetPixel( index );
        movingMean = (double) finitediffimages[1]->GetPixel( index );
  }
  
 
  for(indct=0; indct<hoodlen; indct++)
  {

    IndexType index=hoodIt.GetIndex(indct);
 

    bool inimage=true;
    for (unsigned int dd=0; dd<ImageDimension; dd++)
    {
      if ( index[dd] < 0 || index[dd] > static_cast<typename IndexType::IndexValueType>(imagesize[dd]-1) ) inimage=false;
    }
    
    if (inimage)
    {
        // Get fixed image related information   
      double fixedGradientSquaredMagnitude = 0;
      
      // Note: no need to check the index is within
      // fixed image buffer. This is done by the external filter.
      fixedValue = (double) Superclass::m_FixedImage->GetPixel( index )-fixedMean;
      CovariantVectorType fixedGradient;
      fixedGradient = m_FixedImageGradientCalculator->EvaluateAtIndex( index ); 
      for( j = 0; j < ImageDimension; j++ )
      {
        fixedGradientSquaredMagnitude += vnl_math_sqr( fixedGradient[j] ) * m_FixedImageSpacing[j];
      } 

      // Get moving image related information
      
      typename TDeformationField::PixelType vec = Superclass::m_DeformationField->GetPixel(index);

      for( j = 0; j < ImageDimension; j++ )
      {
        mappedPoint[j] = double( index[j] ) * m_FixedImageSpacing[j] + 
        m_FixedImageOrigin[j];
        mappedPoint[j] += vec[j];
      }
      if( m_MovingImageInterpolator->IsInsideBuffer( mappedPoint ) )
      {
        movingValue = m_MovingImageInterpolator->Evaluate( mappedPoint )-movingMean;
      }
      else
      {
        movingValue = 0.0;
      }

      sff+=fixedValue*fixedValue;
      smm+=movingValue*movingValue;
      sfm+=fixedValue*movingValue;

      for(unsigned int dim=0; dim<ImageDimension; dim++)
      {
        double differential = fixedGradient[dim];
        derivativeF[dim]+= fixedValue  * differential;
        derivativeM[dim]+= movingValue * differential;
      }
    } 
    else return update;
   
  }

  double localCrossCorrelation;
  double factor;
  
    if(m_UseNormalizedCrossCorrelation)
    {
      if (sff*smm !=0.0) 
      {
        localCrossCorrelation = sfm / sqrt( sff * smm );
        factor = 1. / sqrt( sff * smm );
      }
      else if (sff == 0.0 && smm == 0) 
      {
        factor = 0.0;
        localCrossCorrelation = 0.0;
      }
      else 
      {
        factor = 1.0;
        localCrossCorrelation = 1.0;
      }
    }
    else
    {
      if (sff*smm !=0.0) 
      {
        localCrossCorrelation = sfm*sfm / ( sff * smm );
        factor = 1. / ( sff * smm );
      }
      else if (sff == 0.0 && smm == 0) 
      {
        factor = 0.0;
        localCrossCorrelation = 0.0;
      }
      else 
      {
        factor = 1.0;
        localCrossCorrelation = 1.0;
      }
    }
    m_MetricTotal+=localCrossCorrelation;
      Superclass::m_Energy+=localCrossCorrelation;

//      derivative[i] = factor * ( derivativeF[i] - (sfm/smm)*derivativeM[i]);

    double updatenorm = 0.0;
      CovariantVectorType fixedGradient;
      fixedGradient = m_FixedImageGradientCalculator->EvaluateAtIndex( oindex ); 
    for(unsigned int i=0; i<ImageDimension; i++)
    {
//      update[i] = ( derivativeF[i]*factor - localCrossCorrelation*derivativeM[i]);
      update[i] = ( factor - localCrossCorrelation)*fixedGradient[i];
      updatenorm+=update[i]*update[i];
    }
    updatenorm=sqrt(updatenorm);
    if (sqrt(updatenorm) > 1.0)
    {
      update.Fill(0);
    }
//    m_MetricImage->SetPixel(oindex,localCrossCorrelation);
  return update*  Superclass::m_GradientStep;
}




} // end namespace itk

#endif

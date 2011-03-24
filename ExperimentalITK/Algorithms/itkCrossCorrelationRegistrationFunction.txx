/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkCrossCorrelationRegistrationFunction.txx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:13:39 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkCrossCorrelationRegistrationFunction_txx_
#define _itkCrossCorrelationRegistrationFunction_txx_

#include "itkCrossCorrelationRegistrationFunction.h"
#include "itkExceptionObject.h"
#include "vnl/vnl_math.h"
#include "itkImageFileWriter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkMeanImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkWarpImageFilter.h"
namespace itk {

/*
 * Default constructor
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
CrossCorrelationRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::CrossCorrelationRegistrationFunction()
{
  m_AvgMag=0;
  m_Iteration=0;
  RadiusType r;
  unsigned int j;
  for( j = 0; j < ImageDimension; j++ )
    {
    r[j] = 2;
    }
  this->SetRadius(r);
  this->m_Energy = 0.0;
  m_TimeStep = 1.0;
  m_DenominatorThreshold = 1e-9;
  m_IntensityDifferenceThreshold = 0.001;
  Superclass::m_MovingImage = NULL;
  m_MetricGradientImage = NULL;
  Superclass::m_FixedImage = NULL;
  m_FixedImageSpacing.Fill( 1.0 );
  m_FixedImageOrigin.Fill( 0.0 );
  m_FixedImageGradientCalculator = GradientCalculatorType::New();
  binaryimage=NULL;

  m_MovingImageGradientCalculator = GradientCalculatorType::New();

  typename DefaultInterpolatorType::Pointer interp =
    DefaultInterpolatorType::New();

  m_MovingImageInterpolator = static_cast<InterpolatorType*>(
    interp.GetPointer() );

  for (int i=0; i<5; i++) finitediffimages[i]=NULL;

}


/*
 * Standard "PrintSelf" method.
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
void
CrossCorrelationRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::PrintSelf(std::ostream& os, Indent indent) const
{
  
  Superclass::PrintSelf(os, indent);
/*
  os << indent << "MovingImageIterpolator: ";
  os << m_MovingImageInterpolator.GetPointer() << std::endl;
  os << indent << "FixedImageGradientCalculator: ";
  os << m_FixedImageGradientCalculator.GetPointer() << std::endl;
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
CrossCorrelationRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::InitializeIteration()
{
  typedef ImageRegionIteratorWithIndex<MetricImageType> ittype;
  if( !Superclass::m_MovingImage || !Superclass::m_FixedImage || !m_MovingImageInterpolator )
    {
    itkExceptionMacro( << "MovingImage, FixedImage and/or Interpolator not set" );
    throw ExceptionObject(__FILE__,__LINE__);
    }

  // cache fixed image information
  m_FixedImageSpacing    = Superclass::m_FixedImage->GetSpacing();
  m_FixedImageOrigin     = Superclass::m_FixedImage->GetOrigin();

  // Warp the moving image according to the deformation field
		typedef WarpImageFilter<MovingImageType, 
																										MovingImageType, 
																										DeformationFieldType> WarperType;
		typename WarperType::Pointer warper = WarperType::New();

		warper->SetInput( this->m_MovingImage );
		warper->SetDeformationField( this->m_DeformationField );
		warper->SetOutputSpacing( this->m_MovingImage->GetSpacing() );
		warper->SetOutputOrigin( this->m_MovingImage->GetOrigin() );
		warper->Update();
  this->SetMovingImage( warper->GetOutput() );


  // setup gradient calculator
  m_FixedImageGradientCalculator->SetInputImage( Superclass::m_FixedImage );
  m_MovingImageGradientCalculator->SetInputImage( Superclass::m_MovingImage  );
  
  // setup moving image interpolator
  m_MovingImageInterpolator->SetInputImage( Superclass::m_MovingImage );


  unsigned long numpix=1;
  for (int i=0; i<ImageDimension;i++) numpix*=Superclass::m_FixedImage->GetLargestPossibleRegion().GetSize()[i];
  m_MetricTotal=0.0;
  this->m_Energy=0.0;

  typedef itk::MeanImageFilter<MovingImageType, BinaryImageType> MovingFilterType;
  typedef itk::MeanImageFilter<FixedImageType, BinaryImageType> FixedFilterType;
 
  // compute the normalizer
  m_Normalizer      = 0.0;
  for( unsigned int k = 0; k < ImageDimension; k++ )
    {
    m_Normalizer += m_FixedImageSpacing[k] * m_FixedImageSpacing[k];
    }
  m_Normalizer /= static_cast<double>( ImageDimension );
  
  typename FixedImageType::SpacingType spacing=this->GetFixedImage()->GetSpacing();

  bool makeimg=false;
  if ( m_Iteration==0 ) makeimg=true;
  else if (!finitediffimages[0] ) makeimg=true;
  else if ( finitediffimages[0]->GetLargestPossibleRegion().GetSize()[0] != 
	    this->GetFixedImage()->GetLargestPossibleRegion().GetSize()[0] ) makeimg=true;
  
  if (makeimg)
    {  
      finitediffimages[2]=this->MakeImage();
      finitediffimages[3]=this->MakeImage();
      finitediffimages[4]=this->MakeImage();
    }
  
  float sig=15.;
  
  RadiusType r;
  for( int j = 0; j < ImageDimension; j++ )    r[j] = 2;
  
  {
    typename FixedFilterType::Pointer fixedfilter = FixedFilterType::New();
    fixedfilter->SetRadius(this->GetRadius());
    fixedfilter->SetInput(this->GetFixedImage());
    fixedfilter->Update();
    finitediffimages[0]=fixedfilter->GetOutput();       
  }  
  {
    typename MovingFilterType::Pointer movingfilter = MovingFilterType::New();
    movingfilter->SetRadius(this->GetRadius());
    movingfilter->SetInput(this->GetMovingImage());
    movingfilter->Update();
    finitediffimages[1]=movingfilter->GetOutput();       
  }
 

  typedef itk::ImageRegionIteratorWithIndex<MetricImageType> Iterator;
  Iterator outIter(this->finitediffimages[0],this->finitediffimages[0]->GetLargestPossibleRegion() );


  for( outIter.GoToBegin(); !outIter.IsAtEnd(); ++outIter )
    { 
      IndexType index=outIter.GetIndex();
      double fixedValue = this->GetFixedImage()->GetPixel(index)-(double)this->finitediffimages[0]->GetPixel( index );
      double movingValue=this->GetMovingImage()->GetPixel(index)-(double)this->finitediffimages[1]->GetPixel( index );
      this->finitediffimages[0]->SetPixel( index,fixedValue );
      this->finitediffimages[1]->SetPixel( index,movingValue);
    }

  for( outIter.GoToBegin(); !outIter.IsAtEnd(); ++outIter )
    {

      NeighborhoodIterator<MetricImageType> 
	hoodIt( this->GetRadius() ,this->finitediffimages[0] , this->finitediffimages[0]->GetLargestPossibleRegion());
      IndexType oindex = outIter.GetIndex();
      hoodIt.SetLocation(oindex);
      double sff=0.0;
      double smm=0.0;
      double sfm=0.0;
     
      double fixedMean=0;
      double movingMean=0;
     
      PointType mappedPoint;
      unsigned int indct;
      unsigned int hoodlen=hoodIt.Size();
      
      unsigned int inct=0;
      //     double sumj=0,sumi=0;
      typename FixedImageType::SizeType imagesize=this->finitediffimages[0]->GetLargestPossibleRegion().GetSize();
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
	      double fixedValue =(double)this->finitediffimages[0]->GetPixel( index );
	      double movingValue=(double)this->finitediffimages[1]->GetPixel( index );
	      
	      sff+=fixedValue*fixedValue;
	      smm+=movingValue*movingValue;
	      sfm+=fixedValue*movingValue;
	      // sumj+=movingValue;
	      //	     sumi+=fixedValue;
	    }
	}
      
      this->finitediffimages[2]->SetPixel( oindex, sfm );//A
      this->finitediffimages[3]->SetPixel( oindex, sff );//B
      this->finitediffimages[4]->SetPixel( oindex, smm );//C
      //         this->finitediffimages[5]->SetPixel( oindex , sumi);//B*C
      //    this->finitediffimages[6]->SetPixel( oindex , sumj);//B*C
      
    }
  
  //m_FixedImageGradientCalculator->SetInputImage(finitediffimages[0]); 
  
  m_MaxMag=0.0;
  m_MinMag=9.e9;
  m_AvgMag=0.0;  
  m_Iteration++;

}


/*
 * Compute the ncc metric everywhere
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
typename TDeformationField::PixelType
CrossCorrelationRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::ComputeMetricAtPairB(IndexType oindex, typename TDeformationField::PixelType vec)
{
  
  typename TDeformationField::PixelType deriv;
  deriv.Fill(0.0);
  double localCrossCorrelation=1;
  double sff=0.0;
  double smm=0.0;
  double sfm=0.0;
  double fixedValue;
  double movingValue;
  sff=0.0;
  smm=0.0;
  sfm=0.0;
  PointType mappedPoint;
  CovariantVectorType gradI,gradJ;
  
  sfm=finitediffimages[2]->GetPixel(oindex);
  sff=finitediffimages[3]->GetPixel(oindex);
  smm=finitediffimages[4]->GetPixel(oindex);
  
  //double sumi=finitediffimages[5]->GetPixel(oindex);
  //double sumj=finitediffimages[6]->GetPixel(oindex);

  //  for(indct=0; indct<hoodlen; indct++)
    {
      
      IndexType index=oindex;//hoodIt.GetIndex(indct);
      bool inimage=true;
      //      float dist=0.0;
      //for (unsigned int dd=0; dd<ImageDimension; dd++)
      //	{
      //  if ( index[dd] < 0 || index[dd] > static_cast<typename IndexType::IndexValueType>(imagesize[dd]-1) ) inimage=false;
      //  float temp=(float)oindex[dd] -  (float)index[dd];
      //  dist+=temp*temp;
      // }
      if (sff == 0.0) sff=1.0;
      if (smm == 0.0) smm=1.0;
      //if (inimage && sff > 0 & smm > 0 )
      {
	
	gradI = m_FixedImageGradientCalculator->EvaluateAtIndex( index ); 
	//	gradJ = m_MovingImageGradientCalculator->EvaluateAtIndex( index ); 

	float  Ji=finitediffimages[1]->GetPixel(index);
	float  Ii=finitediffimages[0]->GetPixel(index);
	m_TEMP=2.0*sfm/(sff*smm)*( Ji - sfm/sff*Ii );
	for (int qq=0; qq<ImageDimension; qq++) 
	  {
	    deriv[qq]   -=2.0*sfm/(sff*smm)*( Ji - sfm/sff*Ii )*gradI[qq];
	    //	    derivinv[qq]-=2.0*sfm/(sff*smm)*( Ii - sfm/smm*Ji )*gradJ[qq];
	  }
	//for (int qq=0; qq<ImageDimension; qq++) deriv[qq]+=2.0*sfm/(sff*smm)*( sumj - sfm/sff*sumi )*gradI[qq];
	
      }
    
    }
  if (sff*smm !=0.0) localCrossCorrelation = sfm*sfm / ( sff * smm );
  else if (sff == 0.0 && smm == 0) localCrossCorrelation = 1.0;
  else localCrossCorrelation = 1.0;
  
  this->m_Energy-=localCrossCorrelation;
  return deriv;//localCrossCorrelation;

}


/*
 * Compute the ncc metric everywhere
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
typename TDeformationField::PixelType
CrossCorrelationRegistrationFunction<TFixedImage,TMovingImage,TDeformationField>
::ComputeMetricAtPairC(IndexType oindex, typename TDeformationField::PixelType vec)
{
  
  typename TDeformationField::PixelType deriv;
  deriv.Fill(0.0);
  double localCrossCorrelation=1;
  double sff=0.0;
  double smm=0.0;
  double sfm=0.0;
  double fixedValue;
  double movingValue;
  sff=0.0;
  smm=0.0;
  sfm=0.0;
  PointType mappedPoint;
  CovariantVectorType gradI,gradJ;
  
  sfm=finitediffimages[2]->GetPixel(oindex);
  sff=finitediffimages[3]->GetPixel(oindex);
  smm=finitediffimages[4]->GetPixel(oindex);
  
      
  IndexType index=oindex;//hoodIt.GetIndex(indct);
  if (sff == 0.0) sff=1.0;
  if (smm == 0.0) smm=1.0;
  
  ///gradI = m_FixedImageGradientCalculator->EvaluateAtIndex( index ); 
  gradJ = m_MovingImageGradientCalculator->EvaluateAtIndex( index ); 

  float  Ji=finitediffimages[1]->GetPixel(index);
  float  Ii=finitediffimages[0]->GetPixel(index);
  for (int qq=0; qq<ImageDimension; qq++) 
    {
      //deriv[qq]   -=2.0*sfm/(sff*smm)*( Ji - sfm/sff*Ii )*gradI[qq];
      deriv[qq]-=2.0*sfm/(sff*smm)*( Ii - sfm/smm*Ji )*gradJ[qq];
    }


  if (sff*smm !=0.0) localCrossCorrelation = sfm*sfm / ( sff * smm );
  else if (sff == 0.0 && smm == 0) localCrossCorrelation = 1.0;
  else localCrossCorrelation = 1.0;
    
  this->m_Energy-=localCrossCorrelation;
  return deriv;//localCrossCorrelation;

}





} // end namespace itk

#endif

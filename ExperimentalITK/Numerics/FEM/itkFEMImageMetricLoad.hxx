/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: itkFEMImageMetricLoad.hxx,v $
Language:  C++

Date:      $Date: 2008/10/18 00:22:49 $
Version:   $Revision: 1.1.1.1 $

Copyright (c) Insight Software Consortium. All rights reser
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for detail.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkFEMImageMetricLoad_hxx_
#define _itkFEMImageMetricLoad_hxx_

#include "itkFEMImageMetricLoad.h"

namespace itk {
namespace fem {

template<class TMovingImage, class TFixedImage>
ImageMetricLoad<TMovingImage, TFixedImage>
::ImageMetricLoad()
{
  m_MetricRadius.Fill(1);
  m_NumberOfIntegrationPoints = 0;

  // Assign the default transfrom
  m_Transform = DeformationFieldTransformType::New();

  // Set up the default image interpolator
  typedef LinearInterpolateImageFunction
          <MovingImageType, RealType> DefaultImageInterpolatorType;
  typename DefaultImageInterpolatorType::Pointer ImageInterpolator
          = DefaultImageInterpolatorType::New();
  m_ImageInterpolator = ImageInterpolator;           
}

template<class TMovingImage, class TFixedImage>
void
ImageMetricLoad<TMovingImage, TFixedImage>
::InitializeMetric()
{ 
  m_Energy = 0.0;
  m_Gamma = 1.0;

  m_Metric->SetMovingImage(m_MovingImage);
  m_Metric->SetFixedImage(m_FixedImage);

  typename FixedImageType::RegionType region;
  typename FixedImageType::SizeType  size;
  typename FixedImageType::IndexType index;

  size.Fill(1);
  index.Fill(0);
  region.SetSize(size);
  region.SetIndex(index);

  m_Metric->SetFixedImageRegion(region);

//  m_Transform->SetDeformationField(m_DeformationField);
  m_Metric->SetTransform(m_Transform);
  m_ImageInterpolator->SetInputImage(m_MovingImage);
  m_Metric->SetInterpolator(m_ImageInterpolator);

  try 
  {
    m_Metric->Initialize();
  } 
  catch( ExceptionObject & e )
  {
    std::cout << "Metric initialization failed" << std::endl;
    std::cout << "Reason " << e.GetDescription() << std::endl;
  }
}

template<class TMovingImage, class TFixedImage>
typename ImageMetricLoad<TMovingImage , TFixedImage>::FEMVectorType 
ImageMetricLoad<TMovingImage, TFixedImage>
::Fe(FEMVectorType Gpos, FEMVectorType Gsol) 
{
  FEMVectorType OutVec(ImageDimension, 0.0);
 
  for (unsigned int k = 0; k < ImageDimension; k++) 
  {
    if (vnl_math_isnan(Gpos[k]) || vnl_math_isinf(Gpos[k]) ||
        vnl_math_isnan(Gsol[k]) || vnl_math_isinf(Gsol[k]) ||
        fabs(Gpos[k]) > 1.e33   || fabs(Gsol[k]) > 1.e33 ) 
    {
      return OutVec;
    }
  }

  std::cout << "HERE 0" << std::endl;

  typename FixedImageType::RegionType requestedRegion;
  RadiusType regionRadius;
  typename FixedImageType::IndexType tindex;
  typename MovingImageType::IndexType rindex; 
  typename DeformationFieldTransformType::ParametersType parameters(m_Transform->GetNumberOfParameters());

  int lobordercheck, hibordercheck;
  for (unsigned int k = 0; k < ImageDimension; k++ )
  { 
    rindex[k] = (int)(Gpos[k] + Gsol[k] + 0.5);  // where the piece of reference image currently lines up under the above translation
    tindex[k] = (int)(Gpos[k] + 0.5) - (int)(0.5*m_MetricRadius[k]);  // position in reference image
    parameters[k]= Gsol[k]; // this gives the translation by the vector field 

    hibordercheck=(int)tindex[k]+(int)m_MetricRadius[k]-(int)m_MovingImage->GetLargestPossibleRegion().GetSize()[k];
    lobordercheck=(int)tindex[k]-(int)m_MetricRadius[k];
    if (hibordercheck >= 0) regionRadius[k]=m_MetricRadius[k]-(long)hibordercheck-1;
    else if (lobordercheck < 0) regionRadius[k]=m_MetricRadius[k]+(long)lobordercheck;
    else regionRadius[k]=m_MetricRadius[k];
    tindex[k]= (long)(Gpos[k]+0.5)-(long)regionRadius[k]/2;  // position in reference image
  }

// Set the associated region

  requestedRegion.SetSize(regionRadius);
  requestedRegion.SetIndex(tindex);
  m_Metric->SetFixedImageRegion(requestedRegion);

//--------------------------------------------------------
// Get metric values

  typename MetricBaseType::MeasureType     measure;
  typename MetricBaseType::DerivativeType  derivative;

  try
  { 
    m_Metric->GetValueAndDerivative(parameters, measure, derivative);
  }
  catch( ... )
  {
  }


  m_Energy += static_cast<RealType>(measure);
  RealType gmag = 0.0;
  for (unsigned int k = 0; k < ImageDimension; k++ )
  {
    if (lobordercheck < 0 || hibordercheck >=0 
        || vnl_math_isnan(derivative[k]) 
        || vnl_math_isinf(derivative[k])) 
    {
      OutVec[k] = 0.0;
    } 
    else OutVec[k]= (m_MaximizeMetric) ?  m_Gamma*derivative[k]
                                       : -m_Gamma*derivative[k];
    gmag += OutVec[k]*OutVec[k];
  }
  if (gmag == 0.0) gmag = 1.0;
  return OutVec/sqrt(gmag);
}

template<class TMovingImage, class TFixedImage>
typename ImageMetricLoad<TMovingImage , TFixedImage>::RealType 
ImageMetricLoad<TMovingImage, TFixedImage>
::EvaluateMetricGivenSolution(Element::ArrayType* el, RealType step)
{
  RealType energy = 0.0; 
  vnl_vector_fixed<RealType, 2*ImageDimension> InVec;
   
  Element::VectorType ip,shapef;
  Element::MatrixType solmat;
  Element::Float w;
  Element::ArrayType::iterator it;

  for (it = el->begin(); it != el->end(); it++) 
  {
    for (unsigned int i = 0; i < m_NumberOfIntegrationPoints; i++)
    {
      dynamic_cast<Element*>(&*(*it))->GetIntegrationPointAndWeight(i, ip, w, m_NumberOfIntegrationPoints); 
      shapef = (*it)->ShapeFunctions(ip);

      RealType solval, posval;
      Float detJ = (*it)->JacobianDeterminant(ip);
        
      InVec.fill(0.0);
      solmat.set_size((*it)->GetNumberOfNodes()*ImageDimension, 1);
      for (unsigned int j = 0; j < ImageDimension; j++)
      {
        for (unsigned int n = 0; n < (*it)->GetNumberOfNodes(); n++)
        {
          InVec[j] += shapef[n]*(((*it)->GetNodeCoordinates(n))[j]);
          solmat[(n*ImageDimension)+j][0] = 
             ((m_Solution)->GetSolutionValue( (*it)->GetNode(n)->GetDegreeOfFreedom(j), 1)
            +(m_Solution)->GetSolutionValue((*it)->GetNode(n)->GetDegreeOfFreedom(j), 0)*step);
          InVec[j+ImageDimension] += shapef[n] * solmat[(n*ImageDimension)+j][0];   
        }
      }
      try
      {
        RealType val = fabs(this->GetMetricValue(InVec));
        for (unsigned int n = 0; n < (*it)->GetNumberOfNodes(); n++)
        {
          energy += shapef[n]*val*w*detJ;
        }
      }
      catch(itk::ExceptionObject &)
      { 
        energy = 0.0;
      }
    }  
  }
  return fabs(energy*m_Gamma);
}

template<class TMovingImage, class TFixedImage>
typename ImageMetricLoad<TMovingImage, TFixedImage>::RealType
ImageMetricLoad<TMovingImage, TFixedImage>
::GetMetricValue(FEMVectorType InVec) 
{
  typename FixedImageType::RegionType requestedRegion;
  typename FixedImageType::IndexType tindex;
  typename MovingImageType::IndexType rindex;
  RadiusType regionRadius;

  FEMVectorType OutVec(ImageDimension, 0.0); 

  for( unsigned int k = 0; k < ImageDimension; k++ )
  { 
    rindex[k] =(long)(InVec[k]+InVec[k+ImageDimension]+0.5);  // where the piece of reference image currently lines up under the above translation
    tindex[k]= (long)(InVec[k]+0.5)-(long)m_MetricRadius[k]/2;  // position in reference image
    int hibordercheck=(int)tindex[k]+(int)m_MetricRadius[k]-(int)m_MovingImage->GetLargestPossibleRegion().GetSize()[k];
    int lobordercheck=(int)tindex[k]-(int)m_MetricRadius[k];
    if (hibordercheck > 0) regionRadius[k]=m_MetricRadius[k]-(long)hibordercheck-1;
    else if (lobordercheck < 0) regionRadius[k]=m_MetricRadius[k]+(long)lobordercheck;
    else regionRadius[k]=m_MetricRadius[k];  
    tindex[k]= (long)(InVec[k]+0.5)-(long)regionRadius[k]/2;  // position in reference image
  }

// Set the associated region

  requestedRegion.SetSize(regionRadius);
  requestedRegion.SetIndex(tindex);
  m_Metric->SetFixedImageRegion(requestedRegion);

  typename DeformationFieldTransformType::ParametersType params;

  try
  { 
    return static_cast<RealType>(m_Metric->GetValue(params));
  }
  catch( ... )
  {
  }
  return 0.0;
}

template<class TMovingImage, class TFixedImage> 
int ImageMetricLoad<TMovingImage,TFixedImage>::CLID()
{
  static const int CLID_ = FEMOF::Register( ImageMetricLoad::NewB,(std::string("ImageMetricLoad(")
                +typeid(TMovingImage).name()+","+typeid(TFixedImage).name()+")").c_str());
  return CLID_;
}

template<class TMovingImage, class TFixedImage> 
const int ImageMetricLoad<TMovingImage,TFixedImage>
::DummyCLID = ImageMetricLoad<TMovingImage,TFixedImage>::CLID();

} // end namespace fem
} // end namespace itk

#endif

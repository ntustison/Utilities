/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: itkFEMFiniteDifferenceFunctionLoad.hxx,v $
Language:  C++

Date:      $Date: 2008/10/18 00:22:49 $
Version:   $Revision: 1.1.1.1 $

Copyright (c) Insight Software Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for detail.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkFEMFiniteDifferenceFunctionLoad_hxx_
#define _itkFEMFiniteDifferenceFunctionLoad_hxx_

#include "itkFEMFiniteDifferenceFunctionLoad.h"

#include "itkNeighborhoodIterator.h"

namespace itk {
namespace fem {

template<class TMovingImage,class TFixedImage>
FiniteDifferenceFunctionLoad<TMovingImage , TFixedImage>
::FiniteDifferenceFunctionLoad()
{
  m_Metric = NULL;
  m_MetricRadius.Fill(1);
}

template<class TMovingImage,class TFixedImage>
typename FiniteDifferenceFunctionLoad<TMovingImage , TFixedImage>::RealType 
FiniteDifferenceFunctionLoad<TMovingImage , TFixedImage>
::EvaluateMetricGivenSolution(Element::ArrayType* el, RealType step)
{
  return 10.0;  //FIXME
  Float energy=0.0,defe=0.0;

  vnl_vector_fixed<Float,2*ImageDimension> InVec(0.0);
   
  typename Element::VectorType ip,shapef;
  typename Element::MatrixType solmat;
  typename Element::Float w;
 
  typedef typename Element::ArrayType ArrayType;

  ArrayType::iterator elt = el->begin();

  const unsigned int Nnodes=(*elt)->GetNumberOfNodes();

  FEMVectorType Gpos,Gsol;
  Gpos.set_size(ImageDimension); Gpos.fill(0.0);
  Gsol.set_size(ImageDimension); Gsol.fill(0.0);

  solmat.set_size(Nnodes*ImageDimension,1);

  for(  ; elt != el->end(); elt++) 
  {
    for(unsigned int i=0; i<m_NumberOfIntegrationPoints; i++)
    {
      dynamic_cast<Element*>(&*(*elt))->GetIntegrationPointAndWeight(i,ip,w,m_NumberOfIntegrationPoints); // FIXME REMOVE WHEN ELEMENT NEW IS BASE CLASS
      shapef = (*elt)->ShapeFunctions(ip);

      float solval,posval;
      Float detJ=(*elt)->JacobianDeterminant(ip);
        
      for(unsigned int f=0; f<ImageDimension; f++)
      {
        solval=0.0;
        posval=0.0;
        for(unsigned int n=0; n<Nnodes; n++)
        {
          posval+=shapef[n]*(((*elt)->GetNodeCoordinates(n))[f]);
          float nodeval=( (m_Solution)->GetSolutionValue( (*elt)->GetNode(n)->GetDegreeOfFreedom(f), 1)
            +(m_Solution)->GetSolutionValue( (*elt)->GetNode(n)->GetDegreeOfFreedom(f), 0)*step);
      
          solval+=shapef[n] * nodeval;   
          solmat[(n*ImageDimension)+f][0]=nodeval;
        }
        InVec[f]=posval;
        Gpos[f]=posval;
        InVec[f+ImageDimension]=solval;
        Gsol[f]=solval;
      }

      float tempe=0.0;
      try
      {
        this->Fe(Gpos, Gsol); // FIXME
        tempe = fabs(0.0);
      }
      catch( ... )
      { 
      }
      for(unsigned int n=0; n<Nnodes; n++)
      {
        itk::fem::Element::Float temp=shapef[n]*tempe*w*detJ;
        energy += temp;
      }
    }  
  }
  return fabs((double)energy*(double)m_Gamma-(double)defe);
}

template<class TMovingImage,class TFixedImage>
typename FiniteDifferenceFunctionLoad<TMovingImage , TFixedImage>::FEMVectorType 
FiniteDifferenceFunctionLoad<TMovingImage , TFixedImage>
::Fe(FEMVectorType Gpos, FEMVectorType Gsol) 
{
  // We assume the vector input is of size 2*ImageDimension.
  // The 0 to ImageDimension-1 elements contain the position, p,
  // in the reference image.  The next ImageDimension to 2*ImageDimension-1
  // elements contain the value of the vector field at that point, v(p).
  //
  // Thus, we evaluate the derivative at the point p+v(p) with respect to
  // some region of the target (fixed) image by calling the metric with 
  // the translation parameters as provided by the vector field at p.
  //------------------------------------------------------------

  VectorType OutVec;
  FEMVectorType femVec;
  femVec.set_size(ImageDimension);
  femVec.fill(0.0);

  if (!m_Metric || !m_DeformationField || !m_FixedImage || !m_MovingImage)
  { 
    m_Metric->InitializeIteration();
    if (!m_DeformationField || !m_FixedImage || !m_MovingImage)
    {
      return femVec;
    }
  }

  typename DeformationFieldType::IndexType idx;
  for (unsigned int k = 0; k < ImageDimension; k++) 
  {    
    if (vnl_math_isnan(Gpos[k]) || vnl_math_isinf(Gpos[k]) ||
        vnl_math_isnan(Gsol[k]) || vnl_math_isinf(Gsol[k]) ||
        fabs(Gpos[k]) > 1.e33   || fabs(Gsol[k]) > 1.e33   ||
        (Gpos[k]+0.5) < 0.0     ||
        (Gpos[k]+0.5) > m_FixedImage->GetLargestPossibleRegion().GetSize()[k]-1)
    {
      return femVec;
    }
    idx[k] = static_cast<unsigned int>(Gpos[k]+0.5);
  }

  NeighborhoodIterator<DeformationFieldType> nD(m_MetricRadius, 
          m_DeformationField, m_DeformationField->GetLargestPossibleRegion());
  nD.SetLocation(idx);
 
  void* globalData = NULL;
  OutVec = m_Metric->ComputeUpdate(nD, globalData);
  for (unsigned int k = 0; k < ImageDimension; k++) 
  {
    femVec[k] = (vnl_math_isnan(OutVec[k]) || vnl_math_isinf(OutVec[k])) 
              ? 0.0 : OutVec[k];
  }
  return femVec;
}

template<class TMovingImage,class TFixedImage> 
int FiniteDifferenceFunctionLoad<TMovingImage,TFixedImage>
::CLID()
{
  std::string clsnm = std::string("FiniteDifferenceFunctionLoad(")+typeid(TMovingImage).name()+","+typeid(TFixedImage).name()+")";
  static const int CLID_ = FEMOF::Register(FiniteDifferenceFunctionLoad::NewB,clsnm.c_str());
  return CLID_;
}


template<class TMovingImage,class TFixedImage> 
const int 
FiniteDifferenceFunctionLoad<TMovingImage,TFixedImage>
::DummyCLID = FiniteDifferenceFunctionLoad<TMovingImage,TFixedImage>::CLID();


} // end namespace fem
} // end namespace itk

#endif

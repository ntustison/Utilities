/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMFiniteDifferenceFunctionLoad.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:22:49 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkFEMFiniteDifferenceFunctionLoad_h_
#define _itkFEMFiniteDifferenceFunctionLoad_h_

#include "itkFEMLoadElementBase.h"

#include "itkImage.h"
#include "itkPDEDeformableRegistrationFunction.h"

#include "vnl/vnl_math.h"

namespace itk 
{
namespace fem
{

/**
 * \class FiniteDifferenceFunctionLoad
 * \brief General image pair load that uses the itkFiniteDifferenceFunctions.
 *
 * This load computes FEM gravity loads by using derivatives provided 
 * by itkFiniteDifferenceFunctions (e.g. mean squares intensity difference.)
 * The function responsible for this is called Fg, as required by the FEMLoad
 * standards.  It takes a vnl_vector as input.
 * We assume the vector input is of size 2*ImageDimension.
 * The 0 to ImageDimension-1 elements contain the position, p,
 * in the reference (moving) image.  The next ImageDimension to 2*ImageDimension-1
 * elements contain the value of the vector field at that point, v(p).
 * The metrics return both a scalar similarity value and vector-valued derivative.  
 * The derivative is what gives us the force to drive the FEM registration.
 * These values are computed with respect to some region in the Fixed image.
 * This region size may be set by the user by calling SetMetricRadius.
 * As the metric derivative computation evolves, performance should improve
 * and more functionality will be available (such as scale selection).
 */ 
template<class TMovingImage, class TFixedImage> 
class FiniteDifferenceFunctionLoad 
: public LoadElement
{
FEM_CLASS(FiniteDifferenceFunctionLoad, LoadElement)

public:
  typedef TMovingImage                                    MovingImageType;
  typedef TFixedImage                                     FixedImageType;
  typedef typename MovingImageType::Pointer               MovingImageTypePointer;
  typedef typename FixedImageType::Pointer                FixedImageTypePointer;

  itkStaticConstMacro(ImageDimension, unsigned int,
                      MovingImageType::ImageDimension);

  typedef NeighborhoodIterator<MovingImageType>           NeighborhoodIteratorType; 
  typedef typename NeighborhoodIteratorType::RadiusType   RadiusType;

  typedef typename LoadElement::Float                     Float;
  typedef vnl_vector<Float>                               FEMVectorType;

  typedef double                                          RealType;
  typedef Vector<RealType,
                 itkGetStaticConstMacro(ImageDimension)>  VectorType;
  typedef Image<VectorType, 
                itkGetStaticConstMacro(ImageDimension)>   DeformationFieldType;
  typedef typename DeformationFieldType::Pointer          DeformationFieldTypePointer;

  typedef PDEDeformableRegistrationFunction
          <FixedImageType, MovingImageType, 
           DeformationFieldType>                          MetricType; 
  typedef typename MetricType::Pointer                    MetricTypePointer;


  /** Main functions */

  FiniteDifferenceFunctionLoad(); 
  RealType EvaluateMetricGivenSolution(Element::ArrayType*, RealType);
  FEMVectorType Fe(FEMVectorType, FEMVectorType);
 
  static Baseclass* NewFiniteDifferenceFunctionLoad()
    { 
      return new FiniteDifferenceFunctionLoad; 
    }

  /** Set/Get functions */

  void SetMetric(MetricTypePointer metric)
    {
      metric->SetFixedImage(m_FixedImage);
      metric->SetMovingImage(m_MovingImage);
      metric->SetRadius(m_MetricRadius);
      metric->SetDeformationField(m_DeformationField);
      metric->InitializeIteration();
      m_Metric = metric; 
    }
  MetricTypePointer GetMetric()
    {
      return m_Metric;
    }
  void SetMovingImage(MovingImageType* I)
    { 
      m_MovingImage = I; 
      if (m_Metric) m_Metric->SetMovingImage(m_MovingImage);
    }
  MovingImageTypePointer GetMovingImage() 
    { 
      return m_MovingImage; 
    }
  void SetFixedImage(FixedImageType* T)
    { 
      m_FixedImage = T; 
      if (m_Metric) m_Metric->SetFixedImage(m_MovingImage); 
    }
  FixedImageTypePointer GetFixedImage() 
    { 
      return m_FixedImage; 
    }
  void SetMetricRadius(RadiusType R) 
    { 
      m_MetricRadius = R; 
    }    
  RadiusType GetMetricRadius() 
    { 
      return m_MetricRadius; 
    }        
  void SetNumberOfIntegrationPoints(unsigned int N) 
    { 
      m_NumberOfIntegrationPoints = N; 
    }
  unsigned int GetNumberOfIntegrationPoints()
    { 
      return m_NumberOfIntegrationPoints; 
    }
  void SetGamma(RealType s) 
    { 
      m_Gamma = s; 
    }
  RealType GetGamma(RealType s) 
    { 
      return m_Gamma; 
    }
  void SetSolution(Solution::ConstPointer ptr) 
    {  
      m_Solution = ptr; 
    }
  Solution::ConstPointer GetSolution() 
    { 
      return m_Solution; 
    }
  RealType GetSolution(unsigned int i, unsigned int which = 0)
    {  
      return m_Solution->GetSolutionValue(i, which); 
    }
  void SetDeformationField(DeformationFieldTypePointer F) 
    { 
      m_DeformationField = F; 
    }
  DeformationFieldTypePointer GetDeformationField() 
    { 
      return m_DeformationField; 
    }
  RealType GetCurrentEnergy() 
    { 
      return ((m_Metric) ? m_Metric->GetEnergy() : 0.0); 
    }
  void SetCurrentEnergy(RealType e) 
    { 
      if (m_Metric) m_Metric->SetEnergy(e); 
    }

protected:

private:
  MovingImageTypePointer                                  m_MovingImage;
  FixedImageTypePointer                                   m_FixedImage;
  DeformationFieldTypePointer                             m_DeformationField;
  MetricTypePointer                                       m_Metric;
  RadiusType                                              m_MetricRadius; 
  unsigned int                                            m_NumberOfIntegrationPoints;
  typename Solution::ConstPointer                         m_Solution;
  RealType                                                m_Gamma;

  /** Dummy static int that enables automatic registration
      with FEMObjectFactory. */
  static const int DummyCLID;
};
}} // end namespace fem/itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFEMFiniteDifferenceFunctionLoad.txx"
#endif

#endif

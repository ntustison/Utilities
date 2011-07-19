/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBSplineKernelTransform.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:20:03 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBSplineKernelTransform_h
#define __itkBSplineKernelTransform_h

#include "itkKernelTransform.h"

#include "itkBSplineScatteredDataPointSetToImageFilter.h"

namespace itk
{

/** \class BSplineKernelTransform
 *
 * \ingroup Transforms
 *
 */
template <class TScalarType, unsigned int NDimensions> 
class ITK_EXPORT BSplineKernelTransform 
: public KernelTransform<TScalarType, NDimensions>
{
public:
  /** Standard class typedefs. */
  typedef BSplineKernelTransform                       Self;
  typedef KernelTransform<TScalarType, NDimensions>    Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro( BSplineKernelTransform, KernelTransform );

  /** New macro for creation of through a Smart Pointer */
  itkNewMacro( Self );

  /** Dimension of the domain space. */
  itkStaticConstMacro( SpaceDimension, 
    unsigned int, NDimensions );

  /** Scalar type. */
  typedef typename Superclass::ScalarType              ScalarType;

  /** Parameters type. */
  typedef typename Superclass::ParametersType          ParametersType;

  /** Jacobian type. */
  typedef typename Superclass::JacobianType            JacobianType;

  /** Standard coordinate point type for this class. */
  typedef typename Superclass::InputPointType          InputPointType;
  typedef typename Superclass::OutputPointType         OutputPointType;
  
  typedef typename Superclass::PointSetType            PointSetType;
  typedef typename Superclass::PointSetPointer         PointSetPointer;
  typedef typename Superclass::PointsContainer         PointsContainer;
  typedef typename Superclass::PointsIterator          PointsIterator;
  typedef typename Superclass::PointsConstIterator     PointsConstIterator;
  
  
  /** Standard vector type for this class. */
  typedef typename Superclass::InputVectorType         InputVectorType;
  typedef typename Superclass::OutputVectorType        OutputVectorType;
  typedef typename Superclass::VectorSetType           VectorSetType;
  
  /**
   * typedefs for the B-Spline transform generator
   */
  typedef Image<InputVectorType, 
    itkGetStaticConstMacro( SpaceDimension )>          TransformationDomainType;
  typedef PointSet<InputVectorType, 
    itkGetStaticConstMacro( SpaceDimension )>          BSplinePointSetType;
  typedef BSplineScatteredDataPointSetToImageFilter
    <BSplinePointSetType, TransformationDomainType>    BSplineTransformType;
  typedef typename BSplineTransformType::ArrayType     ArrayType;

  typedef typename TransformationDomainType::PointType OriginType;
  typedef typename TransformationDomainType::PointType PointType;
  typedef typename TransformationDomainType::SizeType  SizeType;
  typedef typename 
    TransformationDomainType::SpacingType              SpacingType;

  
  /** Compute the position of point in the new space */
  virtual OutputPointType TransformPoint( 
    const InputPointType &thisPoint ) const;

  /** Compute the Jacobian Matrix of the transformation at one point */
  virtual const JacobianType & GetJacobian( 
    const InputPointType &point ) const;
    
  /** Get/Set functions */
  itkSetMacro( NumberOfControlPoints, ArrayType );
  itkGetConstReferenceMacro( NumberOfControlPoints, ArrayType );  
    
  itkSetMacro( SplineOrder, unsigned int );  
  itkGetConstMacro( SplineOrder, unsigned int );

  itkSetMacro( NumberOfLevels, unsigned int );  
  itkGetConstMacro( NumberOfLevels, unsigned int );

  itkSetMacro( DomainExpansionFactor, float );  
  itkGetConstMacro( DomainExpansionFactor, float );

  itkSetMacro( Size, SizeType );
  itkGetConstMacro( Size, SizeType );

  itkSetMacro( Spacing, SpacingType );
  itkGetConstMacro( Spacing, SpacingType );

  itkSetMacro( Origin, OriginType );
  itkGetConstMacro( Origin, OriginType );

  VectorSetType * GetDisplacements();

protected:
  BSplineKernelTransform();
  virtual ~BSplineKernelTransform();
  void PrintSelf(std::ostream& os, Indent indent) const;

//protected:
// 
//  /** The list of displacements.
//   * d[i] = q[i] - p[i]; */
//  VectorSetPointer m_Displacements;
//
//  /** The list of source landmarks, denoted 'p'. */
//  PointSetPointer m_SourceLandmarks;
//  
//  /** The list of target landmarks, denoted 'q'. */
//  PointSetPointer m_TargetLandmarks;
//
 private:
  BSplineKernelTransform(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  void GenerateBSplineTransformation();
  
  
  typename BSplineTransformType::Pointer               m_BSplineTransform;

  bool                                                 m_BSplineTransformIsCalculated;
  ArrayType                                            m_NumberOfControlPoints;
  unsigned int                                         m_SplineOrder;  
  unsigned int                                         m_NumberOfLevels;

  /**
   * Unlike TPS, B-splines have a finite domain so that the transformation
   * is also defined only over a finite domain.  We provide two options for
   * specifying the domain:
   *   1. spacing, size, and origin are used to define the transformation 
   *      domain.
   *   2. from the target and source landmarks we can calculated a bounding box
   *      encompassing both point sets.  The expansion factor determines how 
   *      much padding is wanted around the boundaries.  
   */  

  float                                                m_DomainExpansionFactor;
  OriginType                                           m_Origin;
  SizeType                                             m_Size;
  SpacingType                                          m_Spacing;

};

} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_BSplineKernelTransform(_, EXPORT, x, y) namespace itk { \
  _(2(class EXPORT BSplineKernelTransform< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef BSplineKernelTransform< ITK_TEMPLATE_2 x > \
                                                  BSplineKernelTransform##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkBSplineKernelTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkBSplineKernelTransform.hxx"
#endif

#endif // __itkBSplineKernelTransform_h

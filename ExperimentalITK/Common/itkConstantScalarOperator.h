/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkConstantScalarOperator.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:20:03 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkConstantScalarOperator_h
#define __itkConstantScalarOperator_h

#include "itkNeighborhoodOperator.h"
namespace itk {

/**
 * \class ConstantScalarOperator
 * \brief A NeighborhoodOperator whose coefficients are a one
 * dimensional, discrete ConstantScalar kernel.
 *
 * ConstantScalarOperator can be used to perform ConstantScalar blurring
 * by taking its inner product with to a Neighborhood
 * (NeighborhooIterator) that is swept across an image region.
 * It is a directional operator. 
 
 * (1) The floating-point scalar of the desired ConstantScalar function.
 
 * \sa NeighborhoodOperator
 * \sa NeighborhoodIterator
 * \sa Neighborhood
 *
 * \ingroup Operators
 */
template<class TPixel,unsigned int VDimension=2,
  class TAllocator = NeighborhoodAllocator<TPixel> >
class ITK_EXPORT ConstantScalarOperator
  : public NeighborhoodOperator<TPixel, VDimension, TAllocator>
{
public:
  /** Standard class typedefs. */
  typedef ConstantScalarOperator Self;
  typedef NeighborhoodOperator<TPixel, VDimension, TAllocator>  Superclass;
  
  /** Constructor. */
  ConstantScalarOperator() : m_Scalar(1) { }

  /** Copy constructor */
  ConstantScalarOperator(const Self &other)
    : NeighborhoodOperator<TPixel, VDimension, TAllocator>(other)
  {
    m_Scalar = other.m_Scalar;
  }

  /** Assignment operator */
  Self &operator=(const Self &other)
  {
    Superclass::operator=(other);
    m_Scalar = other.m_Scalar;
    return *this;
  }
  
  /** Sets the desired scalar of the ConstantScalar kernel. */
  void SetScalar(const double &scalar)
  {  m_Scalar = scalar;  }

  /** Returns the scalar of the ConstantScalar (scale) for the operator. */
  double GetScalar()
    {  return m_Scalar;  }

  /** Prints some debugging information. */
  virtual void PrintSelf(std::ostream &os, Indent i) const
  {
    os << i << "ConstantScalarOperator { this=" << this
       << ", m_Scalar = " << m_Scalar
       << "} "  << std::endl;
    Superclass::PrintSelf(os, i.GetNextIndent());
  }
  
protected:
  typedef typename Superclass::CoefficientVector CoefficientVector;

  /** Calculates operator coefficients. */
  CoefficientVector GenerateCoefficients();

  /**
   * Arranges coefficients spatially in the memory buffer.
   */
  void Fill(const CoefficientVector &c) {};


private:
  /** Desired variance of the discrete ConstantScalar function. */
  double m_Scalar;

  /** For compatibility with itkWarningMacro */
  const char *GetNameOfClass()
    { return "itkConstantScalarOperator"; }
  
};

} // namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkConstantScalarOperator.hxx"
#endif

#endif

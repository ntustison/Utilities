/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkPointSetToPointSetFilter.h,v $
  Language:  C++
  Date:      $$
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkPointSetToPointSetFilter_h
#define __itkPointSetToPointSetFilter_h

#include "itkPointSetSource.h"

namespace itk
{

/** \class PointSetToPointSetFilter
 * \brief Base class for filters that take a graph as an input and output 
 * another graph.
 *
 * PointSetToPointSetFilter is the base class for all process objects that output
 * graph data, and require graph data as input. Specifically, this class
 * defines the SetInput() method for defining the input to a filter.
 * 
 * \ingroup PointSetFilters
 *
 */
template <class TInputPointSet, class TOutputPointSet>
class ITK_EXPORT PointSetToPointSetFilter : public PointSetSource<TOutputPointSet>
{
public:
  /** Standard class typedefs. */
  typedef PointSetToPointSetFilter  Self;
  typedef PointSetSource<TOutputPointSet> Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(PointSetToPointSetFilter,PointSetSource);

  /** Some convenient typedefs. */
  typedef TInputPointSet InputPointSetType;
  typedef typename InputPointSetType::Pointer InputPointSetPointer;
  
  /** Set the point set input of this process object.  */
  virtual void SetInput( const InputPointSetType * pointSet );
  virtual void SetInput( unsigned int, const TInputPointSet * pointSet );
  const InputPointSetType * GetInput(void);
  const InputPointSetType * GetInput(unsigned int idx);

protected:
  PointSetToPointSetFilter();
  ~PointSetToPointSetFilter() {};
  
private:
  PointSetToPointSetFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPointSetToPointSetFilter.txx"
#endif

#endif

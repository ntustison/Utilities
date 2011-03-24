/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGraphDuplicator.h,v $
  Language:  C++
  Date:      $Date: 2008/11/11 03:08:24 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkGraphDuplicator_h
#define __itkGraphDuplicator_h

#include "itkObject.h"
#include "itkGraph.h"

namespace itk
{

/**
 * This helper class creates a herlper graph which is perfect copy of 
 * the input graph.
 */
template <class TInputGraph>            
class ITK_EXPORT GraphDuplicator : public Object 
{
public:
  /** Standard class typedefs. */
  typedef GraphDuplicator Self;
  typedef Object  Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(GraphDuplicator, Object);

  /** Type definitions for the input Graph. */
  typedef TInputGraph  GraphType;
  typedef typename TInputGraph::Pointer  GraphPointer;
  typedef typename TInputGraph::ConstPointer GraphConstPointer;

  /** Set the input Graph. */
  itkSetConstObjectMacro(Input, GraphType);
  
  /** Get the output Graph. */
  itkGetObjectMacro(Output, GraphType);

  /** Compute of the input Graph. */
  void Update(void);

protected:
  GraphDuplicator();
  virtual ~GraphDuplicator() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  GraphDuplicator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  GraphConstPointer   m_Input;
  GraphPointer        m_Output;
  unsigned long       m_InternalGraphTime;
  
};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGraphDuplicator.txx"
#endif

#endif /* __itkGraphDuplicator_h */

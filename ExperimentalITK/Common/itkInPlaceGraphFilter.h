/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkInPlaceGraphFilter.h,v $
  Language:  C++
  Date:      $Date: 2008/11/11 03:08:24 $
  Version:   $Revision: 1.2 $

  Copyright ( c ) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkInPlaceGraphFilter_h
#define __itkInPlaceGraphFilter_h

#include "itkGraphToGraphFilter.h"

namespace itk
{
  
/** \class InPlaceGraphFilter
 * \brief Base class for filters that take an graph as input and overwrite 
 *        that graph as the output.
 *
 * InPlaceGraphFilter is the base class for all process objects whose
 * output graph data is constructed by overwriting the input graph
 * data. In other words, the output bulk data is the same block of
 * memory as the input bulk data.  This filter provides the mechanisms
 * for in place graph processing while maintaining general pipeline
 * mechanics. InPlaceGraphFilters use less memory than standard
 * GraphToGraphFilters because the input buffer is reused as the
 * output buffer.  However, this benefit does not come without a cost.
 * Since the filter overwrites its input, the ownership of the bulk
 * data is transitioned from the input data object to the output data
 * object.  When a data object has multiple consumers with one
 * of the consumers being an in place filter, the in place filter
 * effectively destroys the bulk data for the data object. Upstream
 * filters will then have to re-execute to regenerate the data object's
 * bulk data for the remaining consumers.
 *
 * Since an InPlaceGraphFilter reuses the input bulk data memory for the
 * output bulk data memory, the input graph type must match the output
 * graph type.  If the input and output graph types are not identical,
 * the filter reverts to a traditional GraphToGraphFilter behaviour
 * where an output graph is allocated.  In place operation can also be
 * controlled ( when the input and output graph type match ) via the
 * methods InPlaceOn() and InPlaceOff().
 *
 * Subclasses of InPlaceGraphFilter must take extra care in how they
 * manage memory using ( and perhaps overriding ) the implementations of
 * ReleaseInputs() and AllocateOutputs() provided here.
 *
 * \ingroup GraphFilters
 */
template <class TInputGraph, class TOutputGraph=TInputGraph>
class ITK_EXPORT InPlaceGraphFilter 
: public GraphToGraphFilter<TInputGraph, TOutputGraph>
{
public:
  /** Standard class typedefs. */
  typedef InPlaceGraphFilter  Self;
  typedef GraphToGraphFilter<TInputGraph, TOutputGraph>  Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  
  /** Run-time type information ( and related methods ). */
  itkTypeMacro( InPlaceGraphFilter, GraphToGraphFilter );

  /** Superclass typedefs. */
  typedef typename Superclass::OutputGraphType    OutputGraphType;
  typedef typename Superclass::OutputGraphPointer OutputGraphPointer;

  /** Some convenient typedefs. */
  typedef TInputGraph InputGraphType;
  typedef typename InputGraphType::Pointer        InputGraphPointer;
  typedef typename InputGraphType::ConstPointer   InputGraphConstPointer;

  /** In place operation can be turned on and off. This only has an
   * effect when the input and output graph type match. */
  itkSetMacro( InPlace, bool );
  itkGetMacro( InPlace, bool );
  itkBooleanMacro( InPlace );

  /** Can the filter run in place? To do so, the filter's first input
   * and output must have the same dimension and pixel type. This
   * method can be used in conjunction with the InPlace ivar to
   * determine whether a particular use of the filter is really
   * running in place. Some filters may be able to optimize their
   * operation if the InPlace is true and CanRunInPlace is true. */
   bool CanRunInPlace() const
     {
       return ( typeid( TInputGraph ) == typeid( TOutputGraph ) );
     };

 protected:
  InPlaceGraphFilter();
  ~InPlaceGraphFilter();

  virtual void PrintSelf( std::ostream& os, Indent indent ) const;

  /** The GenerateData method normally allocates the buffers for all
   * of the outputs of a filter. Since InPlaceGraphFilter's can use an
   * overwritten version of the input for its output, the output
   * buffer should not be allocated. When possible, we graft the input
   * to the filter to the output.  If an InPlaceFilter has multiple
   * outputs, then it would need to override this method to graft one
   * of its outputs and allocate the remaining. If a filter is
   * threaded ( i.e. it provides an implementation of
   * ThreadedGenerateData() ), this method is called automatically. If
   * an InPlaceFilter is not threaded ( i.e. it provides an
   * implementation of GenerateData() ), then this method ( or
   * equivalent ) must be called in GenerateData(). */
  virtual void AllocateOutputs();

  /** InPlaceGraphFilter may transfer ownership of the input bulk data
   * to the output object.  Once the output object owns the bulk data
   * ( done in AllocateOutputs() ), the input object must release its
   * hold on the bulk data.  ProcessObject::ReleaseInputs() only
   * releases the input bulk data when the user has set the
   * ReleaseDataFlag.  InPlaceGraphFilter::ReleaseInputs() also
   * releases the input that it has overwritten.
   *
   * \sa ProcessObject::ReleaseInputs() */
  virtual void ReleaseInputs(); 

private:
  InPlaceGraphFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  bool m_InPlace;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkInPlaceGraphFilter.txx"
#endif

#endif

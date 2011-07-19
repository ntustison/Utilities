/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBinaryBoxStructuringElement.h,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:20:03 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBinaryBoxStructuringElement_h
#define __itkBinaryBoxStructuringElement_h

#include "itkNeighborhood.h"

namespace itk {
  
/** \class BinaryBoxStructuringElement
 * \brief A Neighborhood that represents a box structuring element 
 *        with binary elements.
 *
 *
 * \sa Neighborhood
 * \sa MorphologyImageFilter
 * \sa BinaryDilateImageFilter
 * \sa BinaryErodeImageFilter
 * 
 * \ingroup Operators
 * \ingroup ImageIterators
 */

template<class TPixel, unsigned int VDimension = 2,
         class TAllocator = NeighborhoodAllocator<TPixel> >
class ITK_EXPORT BinaryBoxStructuringElement
  : public Neighborhood<TPixel, VDimension, TAllocator>
{
public:
  /** Standard class typedefs. */
  typedef BinaryBoxStructuringElement                Self;
  typedef Neighborhood<TPixel, VDimension, TAllocator> Superclass;

  /** External support for allocator type. */
  typedef TAllocator AllocatorType;

  /** External support for dimensionality. */
  itkStaticConstMacro(NeighborhoodDimension, unsigned int, VDimension);
  
  /** External support for pixel type. */
  typedef TPixel PixelType;
  
  /** Iterator typedef support. Note the naming is intentional, i.e.,
  * ::iterator and ::const_iterator, because the allocator may be a
  * vnl object or other type, which uses this form. */
  typedef typename AllocatorType::iterator       Iterator;
  typedef typename AllocatorType::const_iterator ConstIterator;
  
  /** Size and value typedef support. */
  typedef typename Superclass::SizeType      SizeType;
  typedef typename Superclass::SizeValueType SizeValueType;

  /** Offset and value typedef support. */
  typedef typename Superclass::OffsetType      OffsetType;
  typedef typename OffsetType::OffsetValueType OffsetValueType;
  
  /** Radius typedef support. */
  typedef typename Superclass::RadiusType RadiusType;

  /** External slice iterator type typedef support. */
  typedef SliceIterator<TPixel, Self> SliceIteratorType;
  
  /** Default constructor. */
  BinaryBoxStructuringElement() {}

  /** Default destructor. */
  virtual ~BinaryBoxStructuringElement() {}
    
  /** Copy constructor. */
  BinaryBoxStructuringElement(const Self& other)
    : Neighborhood<TPixel, VDimension, TAllocator>(other)
    {
    }

  /** Assignment operator. */
  Self &operator=(const Self& other)
    {
    Superclass::operator=(other);
    return *this;
    }

  /** Build the structuring element */
  void CreateStructuringElement();   
  
protected:
  
private:

};

} // namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_BinaryBoxStructuringElement(_, EXPORT, x, y) namespace itk { \
  _(2(class EXPORT BinaryBoxStructuringElement< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef BinaryBoxStructuringElement< ITK_TEMPLATE_2 x > \
                                                  BinaryBoxStructuringElement##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkBinaryBoxStructuringElement+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkBinaryBoxStructuringElement.hxx"
#endif

#endif

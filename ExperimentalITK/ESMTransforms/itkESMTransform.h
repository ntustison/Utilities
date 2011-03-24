#ifndef MZ_ESMTransform_H_
#define MZ_ESMTransform_H_

#include <iostream>

#include "itkMatrix.h"
#include "itkTransform.h"
#include "itkExceptionObject.h"
#include "itkMacro.h"

namespace itk
{


/**
 *
 * \ingroup Transforms
 *
 */

template <
   class TScalarType=double,   // Data type for scalars
   unsigned int NDimensions=2> // Number of dimensions in the space
class ESMTransform
   : public Transform< TScalarType, NDimensions, NDimensions >
{
public:
   /** Standard typedefs   */
   typedef ESMTransform                          Self;
   typedef Transform< TScalarType,
                      NDimensions,
                      NDimensions >              Superclass;
   typedef SmartPointer<Self>                    Pointer;
   typedef SmartPointer<const Self>              ConstPointer;

   /** Run-time type information (and related methods).   */
   itkTypeMacro( ESMTransform, Transform );

   /** Dimension of the domain space. */
   itkStaticConstMacro(SpaceDimension, unsigned int, NDimensions);

   /** Parameters Type   */
   typedef typename Superclass::ParametersType                  ParametersType;

   /** Jacobian Type   */
   typedef typename Superclass::JacobianType                    JacobianType;

   /** Standard scalar type for this class */
   typedef typename Superclass::ScalarType                      ScalarType;

   /** Standard vector type for this class   */
   typedef Vector<TScalarType,
                  itkGetStaticConstMacro(SpaceDimension)>       VectorType;

   /** Standard covariant vector type for this class   */
   typedef CovariantVector<TScalarType,
                           itkGetStaticConstMacro(SpaceDimension)>
                                                                CovariantVectorType;

   /** Standard vnl_vector type for this class   */
   typedef vnl_vector_fixed<TScalarType,
                            itkGetStaticConstMacro(SpaceDimension)>
                                                                VnlVectorType;

   /** Standard coordinate point type for this class   */
   typedef Point<TScalarType,
                 itkGetStaticConstMacro(SpaceDimension)>        PointType;

   /** Base inverse ESM transform type. This type should not be changed to the
    * concrete inverse transform type or inheritance would be lost.*/
   typedef Self                                                 InverseESMTransformBaseType;
   typedef typename InverseESMTransformBaseType::Pointer        InverseESMTransformBasePointer;

   /** Return an inverse of this transform.
    */
   virtual InverseESMTransformBasePointer GetInverseESMTransform() const = 0;

   /** Use the given update vector to update the current transform
    *  Typically this will be done with    this <- this o exp ( update )
    */
   virtual void IncrementalUpdate( const ParametersType & updatevector ) = 0;

   /** Get the derivative of exp ( update ) ( point ) w.r.t the update at zero.
    *  This is a SpaceDimension x ParametersDimension matrix.
    */
   virtual const JacobianType & GetIncrementalUpdateJacobian(const PointType  & pt) const = 0;

   /** Get the derivative of s ( point ) w.r.t the point coordinates.
    *  This is a SpaceDimension x SpaceDimension matrix.
    */
   virtual const JacobianType & GetSpatialJacobian(const PointType  & pt) const = 0;

protected:
   ESMTransform(unsigned int dimension,unsigned int numberOfParameters)
      : Superclass(dimension,numberOfParameters)
      , m_IncrementalUpdateJacobian(dimension,numberOfParameters)
      , m_SpatialJacobian(dimension,dimension)
   {
   }

   ESMTransform()
      : Superclass()
      , m_IncrementalUpdateJacobian(NDimensions,1)
      , m_SpatialJacobian(NDimensions,NDimensions)
   {
   }

   /** Print contents of an ESMTransform */
   void PrintSelf(std::ostream &os, Indent indent) const
   {
      Superclass::PrintSelf(os,indent);

      os << indent << "m_IncrementalUpdateJacobian: " << m_IncrementalUpdateJacobian << std::endl;
      os << indent << "m_SpatialJacobian: " << m_SpatialJacobian << std::endl;
   }

   mutable JacobianType       m_IncrementalUpdateJacobian;
   mutable JacobianType       m_SpatialJacobian;

private:

   ESMTransform(const Self & other);
   const Self & operator=( const Self & );

}; //class ESMTransform

}  // namespace itk

#endif /* MZ_ESMTransform_H_ */


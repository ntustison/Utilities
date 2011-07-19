#ifndef MZ_ESMRigid2DTransform_H_
#define MZ_ESMRigid2DTransform_H_

#include "itkESMMatrixOffsetTransformBase.h"

#include <iostream>
#include <itkExceptionObject.h>

namespace itk
{

/** \class ESMRigid2DTransform
 * \brief ESMRigid2DTransform of a vector space (e.g. space coordinates)
 *
 * This transform applies a rigid transformation in 2D space.
 * The transform is specified as a rotation around a arbitrary center
 * and is followed by a translation.
 *
 * The parameters for this transform can be set either using
 * individual Set methods or in serialized form using
 * SetParameters() and SetFixedParameters().
 *
 * The serialization of the optimizable parameters is an array of 3 elements
 * ordered as follows:
 * p[0] = angle
 * p[1] = x component of the translation
 * p[2] = y component of the translation
 *
 * The serialization of the fixed parameters is an array of 2 elements
 * ordered as follows:
 * p[0] = x coordinate of the center
 * p[1] = y coordinate of the center
 *
 * Access methods for the center, translation and underlying matrix
 * offset vectors are documented in the superclass ESMMatrixOffsetTransformBase.
 *
 * \sa Transfrom
 * \sa ESMMatrixOffsetTransformBase
 *
 * \ingroup Transforms
 */
template < class TScalarType=double >    // Data type for scalars (float or double)
class ITK_EXPORT ESMRigid2DTransform :
      public ESMMatrixOffsetTransformBase< TScalarType, 2> // Dimensions of input and output spaces
{
public:
   /** Standard class typedefs. */
   typedef ESMRigid2DTransform                            Self;
   typedef ESMMatrixOffsetTransformBase< TScalarType, 2 > Superclass;
   typedef SmartPointer<Self>                             Pointer;
   typedef SmartPointer<const Self>                       ConstPointer;

   /** Run-time type information (and related methods). */
   itkTypeMacro( ESMRigid2DTransform, ESMMatrixOffsetTransformBase );

   /** New macro for creation of through a Smart Pointer */
   itkNewMacro( Self );

   /** Dimension of the space. */
   itkStaticConstMacro(SpaceDimension, unsigned int, 2);
   itkStaticConstMacro(ParametersDimension, unsigned int, 3);

   /** Scalar type. */
   typedef typename Superclass::ScalarType      ScalarType;

   /** Parameters type. */
   typedef typename Superclass::ParametersType  ParametersType;

   /** Jacobian type. */
   typedef typename Superclass::JacobianType    JacobianType;

   /// Standard matrix type for this class
   typedef typename Superclass::MatrixType      MatrixType;

   /// Standard vector type for this class
   typedef typename Superclass::OffsetType      OffsetType;

   /// Standard vector type for this class
   typedef typename Superclass::VectorType      VectorType;

   /// Standard covariant vector type for this class
   typedef typename Superclass::CovariantVectorType CovariantVectorType;

   /// Standard vnl_vector type for this class
   typedef typename Superclass::VnlVectorType   VnlVectorType;

   /// Standard coordinate point type for this class
   typedef typename Superclass::PointType       PointType;

   /** Base inverse transform type. This type should not be changed to the
    * concrete inverse transform type or inheritance would be lost.*/
   typedef typename Superclass::InverseESMTransformBaseType InverseESMTransformBaseType;
   typedef typename InverseESMTransformBaseType::Pointer    InverseESMTransformBasePointer;

   /**
    * Set the rotation Matrix of a Rigid2D Transform
    *
    * This method sets the 2x2 matrix representing the rotation
    * in the transform.  The Matrix is expected to be orthogonal
    * with a certain tolerance.
    *
    * \warning This method will throw an exception is the matrix
    * provided as argument is not orthogonal.
    *
    * \sa ESMMatrixOffsetTransformBase::SetMatrix()
    */
   virtual void SetMatrix( const MatrixType & matrix );

   /**
    * Set/Get the rotation matrix. These methods are old and are
    * retained for backward compatibility. Instead, use SetMatrix()
    * GetMatrix().
    */
   virtual void SetRotationMatrix(const MatrixType &matrix)
   {
      this->SetMatrix( matrix );
   }
   const MatrixType & GetRotationMatrix() const
   {
      return this->GetMatrix();
   }


   /**
    * Compose the transformation with a translation
    *
    * This method modifies self to include a translation of the
    * origin.  The translation is precomposed with self if pre is
    * true, and postcomposed otherwise.
    */
   void Translate(const OffsetType &offset, bool pre=false);

   /** Set/Get the angle of rotation in radians */
   void SetAngle(TScalarType angle);
   itkGetConstReferenceMacro( Angle, TScalarType );

   /** Set the angle of rotation in degrees. */
   void SetAngleInDegrees(TScalarType angle);

   /** Set/Get the angle of rotation in radians. These methods
    * are old and are retained for backward compatibility.
    * Instead, use SetAngle() and GetAngle(). */
   void SetRotation(TScalarType angle)
   {
      this->SetAngle(angle);
   }

   virtual const TScalarType & GetRotation() const
   {
      return m_Angle;
   }

   /** Set the transformation from a container of parameters
    * This is typically used by optimizers.
    * There are 3 parameters. The first one represents the
    * angle of rotation in radians and the last two represents the translation.
    * The center of rotation is fixed.
    *
    * \sa Transform::SetParameters()
    * \sa Transform::SetFixedParameters() */
   void SetParameters( const ParametersType & parameters );

   /** Get the parameters that uniquely define the transform
    * This is typically used by optimizers.
    * There are 3 parameters. The first one represents the
    * angle or rotation in radians and the last two represents the translation.
    * The center of rotation is fixed.
    *
    * \sa Transform::GetParameters()
    * \sa Transform::GetFixedParameters() */
   const ParametersType & GetParameters( void ) const;

   /** This method computes the Jacobian matrix of the transformation
    * at a given input point.
    *
    * \sa Transform::GetJacobian() */
   const JacobianType & GetJacobian(const PointType  &point ) const;

   /** Return an inverse of this transform. */
   virtual InverseESMTransformBasePointer GetInverseESMTransform() const;

   /**
    * This method creates and returns a new ESMRigid2DTransform object
    * which has the same parameters.
    */
   void CloneTo( Pointer & clone ) const;

   /** Reset the parameters to create and identity transform. */
   virtual void SetIdentity(void);

   /** Use the given update vector to update the current transform
    *  Typically this will be done with    this <- this o exp ( update )
    */
   virtual void IncrementalUpdate( const ParametersType & updatevector );

   virtual const JacobianType & GetIncrementalUpdateJacobian(const PointType  & point ) const;

protected:
   ESMRigid2DTransform();
   ESMRigid2DTransform( unsigned int outputSpaceDimension,
                        unsigned int parametersDimension);

   ~ESMRigid2DTransform();

   /**
    * Print contents of an ESMRigid2DTransform
    */
   void PrintSelf(std::ostream &os, Indent indent) const;

   /** Compute the matrix from angle. This is used in Set methods
    * to update the underlying matrix whenever a transform parameter
    * is changed. */
   virtual void ComputeMatrix(void);

   /** Compute the angle from the matrix. This is used to compute
    * transform parameters from a given matrix. This is used in
    * ESMMatrixOffsetTransformBase::Compose() and
    * ESMMatrixOffsetTransformBase::GetInverse(). */
   virtual void ComputeMatrixParameters(void);

   /** Update angle without recomputation of other internal variables. */
   void SetVarAngle( TScalarType angle )
   {
      m_Angle = angle;
   }

private:
   ESMRigid2DTransform(const Self&); //purposely not implemented
   void operator=(const Self&); //purposely not implemented

   TScalarType         m_Angle;

}; //class ESMRigid2DTransform

}  // namespace itk

#include "itkESMRigid2DTransform.hxx"

#endif /* MZ_ESMRigid2DTransform_H_ */

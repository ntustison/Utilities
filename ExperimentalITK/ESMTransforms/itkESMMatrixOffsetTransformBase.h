#ifndef MZ_ESMMatrixOffsetTransformBase_H_
#define MZ_ESMMatrixOffsetTransformBase_H_

#include <iostream>

#include "itkMatrix.h"
#include "itkESMTransform.h"
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
class ESMMatrixOffsetTransformBase
   : public ESMTransform< TScalarType, NDimensions >
{
public:
   /** Standard typedefs   */
   typedef ESMMatrixOffsetTransformBase          Self;
   typedef ESMTransform< TScalarType,
                         NDimensions >           Superclass;
   typedef SmartPointer<Self>                    Pointer;
   typedef SmartPointer<const Self>              ConstPointer;

   /** Run-time type information (and related methods).   */
   itkTypeMacro( ESMMatrixOffsetTransformBase, ESMTransform );

   /** New macro for creation of through a Smart Pointer   */
   itkNewMacro( Self );

   /** Dimension of the domain space. */
   itkStaticConstMacro(SpaceDimension, unsigned int, NDimensions);
   itkStaticConstMacro(ParametersDimension, unsigned int,
                       NDimensions*(NDimensions+1));


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

   /** Standard matrix type for this class   */
   typedef Matrix<TScalarType, itkGetStaticConstMacro(SpaceDimension),
                  itkGetStaticConstMacro(SpaceDimension)>       MatrixType;

   /** Standard inverse matrix type for this class   */
   typedef Matrix<TScalarType, itkGetStaticConstMacro(SpaceDimension),
                  itkGetStaticConstMacro(SpaceDimension)>       InverseMatrixType;

   typedef vnl_matrix_fixed<TScalarType,
                            itkGetStaticConstMacro(SpaceDimension)+1,
                            itkGetStaticConstMacro(SpaceDimension)+1>
                                                                VnlIncrementMatrixType;

   typedef PointType                                            CenterType;

   typedef VectorType                                           OffsetType;

   typedef VectorType                                           TranslationType;

   /** Base inverse ESM transform type. This type should not be changed to the
    * concrete inverse transform type or inheritance would be lost.*/
   typedef typename Superclass::InverseESMTransformBaseType     InverseESMTransformBaseType;
   typedef typename Superclass::InverseESMTransformBasePointer  InverseESMTransformBasePointer;

   /** Set the transformation to an Identity
    *
    * This sets the matrix to identity and the Offset to null. */
   virtual void SetIdentity( void );

   /** Set matrix of an ESMMatrixOffsetTransformBase
    *
    * This method sets the matrix of an ESMMatrixOffsetTransformBase to a
    * value specified by the user.
    *
    * This updates the Offset wrt to current translation
    * and center.  See the warning regarding offset-versus-translation
    * in the documentation for SetCenter.
    *
    * To define an affine transform, you must set the matrix,
    * center, and translation OR the matrix and offset */
   virtual void SetMatrix(const MatrixType &matrix)
   {
      m_Matrix = matrix;
      this->ComputeOffset();
      this->ComputeMatrixParameters();
      m_MatrixMTime.Modified();
      this->Modified();
      return;
   }

   /** Get matrix of an ESMMatrixOffsetTransformBase
    *
    * This method returns the value of the matrix of the
    * ESMMatrixOffsetTransformBase.
    * To define an affine transform, you must set the matrix,
    * center, and translation OR the matrix and offset */
   const MatrixType & GetMatrix() const
   {
      return m_Matrix;
   }

   /** Set offset (origin) of an MatrixOffset TransformBase.
    *
    * This method sets the offset of an ESMMatrixOffsetTransformBase to a
    * value specified by the user.
    * This updates Translation wrt current center.  See the warning regarding
    * offset-versus-translation in the documentation for SetCenter.
    * To define an affine transform, you must set the matrix,
    * center, and translation OR the matrix and offset */
   void SetOffset(const VectorType &offset)
   {
      m_Offset = offset;
      this->ComputeTranslation();
      this->Modified();
      return;
   }

   /** Get offset of an ESMMatrixOffsetTransformBase
    *
    * This method returns the offset value of the ESMMatrixOffsetTransformBase.
    * To define an affine transform, you must set the matrix,
    * center, and translation OR the matrix and offset */
   const VectorType & GetOffset(void) const
   {
      return m_Offset;
   }

   /** Set center of rotation of an ESMMatrixOffsetTransformBase
    *
    * This method sets the center of rotation of an ESMMatrixOffsetTransformBase
    * to a fixed point - for most transforms derived from this class,
    * this point is not a "parameter" of the transform - the exception is that
    * "centered" transforms have center as a parameter during optimization.
    *
    * This method updates offset wrt to current translation and matrix.
    * That is, changing the center changes the transform!
    *
    * WARNING: When using the Center, we strongly recommend only changing the
    * matrix and translation to define a transform.   Changing a transform's
    * center, changes the mapping between spaces - specifically, translation is
    * not changed with respect to that new center, and so the offset is updated
    * to * maintain the consistency with translation.   If a center is not used,
    * or is set before the matrix and the offset, then it is safe to change the
    * offset directly.
    *        As a rule of thumb, if you wish to set the center explicitly, set
    * before Offset computations are done.
    *
    * To define an affine transform, you must set the matrix,
    * center, and translation OR the matrix and offset */
   void SetCenter(const PointType & center)
   {
      m_Center = center;
      this->ComputeOffset();
      this->Modified();
      return;
   }

   /** Get center of rotation of the ESMMatrixOffsetTransformBase
    *
    * This method returns the point used as the fixed
    * center of rotation for the ESMMatrixOffsetTransformBase.
    * To define an affine transform, you must set the matrix,
    * center, and translation OR the matrix and offset */
   const PointType & GetCenter() const
   {
      return m_Center;
   }

   /** Set translation of an ESMMatrixOffsetTransformBase
    *
    * This method sets the translation of an ESMMatrixOffsetTransformBase.
    * This updates Offset to reflect current translation.
    * To define an affine transform, you must set the matrix,
    * center, and translation OR the matrix and offset */
   void SetTranslation(const VectorType & translation)
   {
      m_Translation = translation;
      this->ComputeOffset();
      this->Modified();
      return;
   }

   /** Get translation component of the ESMMatrixOffsetTransformBase
    *
    * This method returns the translation used after rotation
    * about the center point.
    * To define an affine transform, you must set the matrix,
    * center, and translation OR the matrix and offset */
   const VectorType & GetTranslation(void) const
   {
      return m_Translation;
   }


   /** Set the transformation from a container of parameters.
    * The first (NOutputDimension x NInputDimension) parameters define the
    * matrix and the last NOutputDimension parameters the translation.
    * Offset is updated based on current center. */
   void SetParameters( const ParametersType & parameters );

   /** Get the Transformation Parameters. */
   const ParametersType& GetParameters(void) const;

   /** Set the fixed parameters and update internal transformation. */
   virtual void SetFixedParameters( const ParametersType & );

   /** Get the Fixed Parameters. */
   virtual const ParametersType& GetFixedParameters(void) const;


   /** Compose with another ESMMatrixOffsetTransformBase
    *
    * This method composes self with another ESMMatrixOffsetTransformBase of the
    * same dimension, modifying self to be the composition of self
    * and other.  If the argument pre is true, then other is
    * precomposed with self; that is, the resulting transformation
    * consists of first applying other to the source, followed by
    * self.  If pre is false or omitted, then other is post-composed
    * with self; that is the resulting transformation consists of
    * first applying self to the source, followed by other.
    * This updates the Translation based on current center. */
   void Compose(const Self * other, bool pre=0);

   /** Transform by an affine transformation
    *
    * This method applies the affine transform given by self to a
    * given point or vector, returning the transformed point or
    * vector.  The TransformPoint method transforms its argument as
    * an affine point, whereas the TransformVector method transforms
    * its argument as a vector. */
   PointType     TransformPoint(const PointType & point) const;
   VectorType    TransformVector(const VectorType & vector) const;
   VnlVectorType TransformVector(const VnlVectorType & vector) const;
   CovariantVectorType TransformCovariantVector(
      const CovariantVectorType &vector) const;

   /** Compute the Jacobian of the transformation
    *
    * This method computes the Jacobian matrix of the transformation.
    * given point or vector, returning the transformed point or
    * vector. The rank of the Jacobian will also indicate if the transform
    * is invertible at this point. */
   const JacobianType & GetJacobian(const PointType & point ) const;

   /** Return an inverse of this transform. */
   virtual InverseESMTransformBasePointer GetInverseESMTransform() const;


   /** Indicates that this transform is linear. That is, given two
    * points P and Q, and scalar coefficients a and b, then
    *
    *           T( a*P + b*Q ) = a * T(P) + b * T(Q)
    */
   virtual bool IsLinear() const
   {
      return true;
   }

   /** Use the given update vector to update the current transform
    *  Typically this will be done with    this <- this o exp ( update )
    */
   virtual void IncrementalUpdate( const ParametersType & updatevector );

   /** Get the derivative of exp ( update ) ( point ) w.r.t the update at zero.
    *  This is a SpaceDimension x ParametersDimension matrix.
    */
   virtual const JacobianType & GetIncrementalUpdateJacobian(const PointType  & point ) const;

   /** Get the derivative of s ( point ) w.r.t the point coordinates.
    *  This is a SpaceDimension x SpaceDimension matrix.
    */
   virtual const JacobianType & GetSpatialJacobian(const PointType  & ) const;

protected:
   /** Construct an ESMMatrixOffsetTransformBase object
    *
    * This method constructs a new ESMMatrixOffsetTransformBase object and
    * initializes the matrix and offset parts of the transformation
    * to values specified by the caller.  If the arguments are
    * omitted, then the ESMMatrixOffsetTransformBase is initialized to an identity
    * transformation in the appropriate number of dimensions.   **/
   ESMMatrixOffsetTransformBase(const MatrixType &matrix,
                                const VectorType &offset);
   ESMMatrixOffsetTransformBase(unsigned int outputDims,
                                unsigned int paramDims);
   ESMMatrixOffsetTransformBase();

   /** Destroy an ESMMatrixOffsetTransformBase object   **/
   virtual ~ESMMatrixOffsetTransformBase();

   /** Print contents of an ESMMatrixOffsetTransformBase */
   void PrintSelf(std::ostream &s, Indent indent) const;

   const InverseMatrixType & GetInverseMatrix( void ) const;

   const InverseMatrixType & GetVarInverseMatrix( void ) const
   {
      return m_InverseMatrix;
   }

   void SetVarInverseMatrix(const InverseMatrixType & matrix) const
   {
      m_InverseMatrix = matrix;
      m_InverseMatrixMTime.Modified();
   }

   bool InverseMatrixIsOld(void) const
   {
      if(m_MatrixMTime != m_InverseMatrixMTime)
      {
         return true;
      }
      else
      {
         return false;
      }
   }

   virtual void ComputeMatrixParameters(void);

   virtual void ComputeMatrix(void);

   virtual void ComputeTranslation(void);

   virtual void ComputeOffset(void);

   MatrixType                  m_Matrix;         // Matrix of the transformation
   VectorType                  m_Offset;         // Offset of the transformation
   PointType                   m_Center;
   VectorType                  m_Translation;

   VnlIncrementMatrixType      m_LogIncMatrix;   // Matrix for the incremental update

   /** To avoid recomputation of the inverse if not needed */
   TimeStamp                   m_MatrixMTime;

private:
   ESMMatrixOffsetTransformBase(const Self & other);
   const Self & operator=( const Self & );

   void FillJacobian( JacobianType & jac, const PointType & point ) const;

   mutable InverseMatrixType   m_InverseMatrix;  // Inverse of the matrix
   mutable bool                m_Singular;       // Is m_Inverse singular?

   /** To avoid recomputation of the inverse if not needed */
   mutable TimeStamp           m_InverseMatrixMTime;

}; //class ESMMatrixOffsetTransformBase

}  // namespace itk

#include "itkESMMatrixOffsetTransformBase.txx"

#endif /* MZ_ESMMatrixOffsetTransformBase_H_ */


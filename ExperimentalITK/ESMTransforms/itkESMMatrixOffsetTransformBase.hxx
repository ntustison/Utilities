#ifndef MZ_itkESMMatrixOffsetTransformBase_TXX_
#define MZ_itkESMMatrixOffsetTransformBase_TXX_

#include "vnl_sd_matrix_tools.h"

#include <itkNumericTraits.h>
#include "itkESMMatrixOffsetTransformBase.h"

namespace itk
{

// Constructor with default arguments
template<class TScalarType, unsigned int NDimensions>
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::ESMMatrixOffsetTransformBase()
   : Superclass(SpaceDimension, ParametersDimension)
{
   m_Matrix.SetIdentity();
   m_MatrixMTime.Modified();
   m_Offset.Fill( 0 );
   m_Center.Fill( 0 );
   m_Translation.Fill( 0 );
   m_Singular = false;
   m_InverseMatrix.SetIdentity();
   m_InverseMatrixMTime = m_MatrixMTime;
   m_LogIncMatrix.fill( 0.0 );
   this->m_FixedParameters.SetSize ( NDimensions );
   this->m_FixedParameters.Fill ( 0.0 );
}


// Constructor with default arguments
template<class TScalarType, unsigned int NDimensions>
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::ESMMatrixOffsetTransformBase( unsigned int dimension,
                                unsigned int paramDims   )
   : Superclass(dimension, paramDims)
{
   m_Matrix.SetIdentity();
   m_MatrixMTime.Modified();
   m_Offset.Fill( 0 );
   m_Center.Fill( 0 );
   m_Translation.Fill( 0 );
   m_Singular = false;
   m_InverseMatrix.SetIdentity();
   m_InverseMatrixMTime = m_MatrixMTime;
   m_LogIncMatrix.fill( 0.0 );
}



// Constructor with explicit arguments
template<class TScalarType, unsigned int NDimensions>
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::ESMMatrixOffsetTransformBase(const MatrixType &matrix,
                               const VectorType &offset)
{
   m_Matrix = matrix;
   m_MatrixMTime.Modified();
   m_Offset = offset;
   m_Center.Fill( 0 );
   m_Translation.Fill(0);
   for(unsigned int i=0; i<NDimensions; i++)
   {
      m_Translation[i] = offset[i];
   }
   this->ComputeMatrixParameters();
   m_LogIncMatrix.set_identity();
}



// Destructor
template<class TScalarType, unsigned int NDimensions>
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::~ESMMatrixOffsetTransformBase()
{
   return;
}



// Print self
template<class TScalarType, unsigned int NDimensions>
void
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::PrintSelf(std::ostream &os, Indent indent) const
{
   Superclass::PrintSelf(os,indent);

   unsigned int i, j;

   os << indent << "Matrix: " << std::endl;
   for (i = 0; i < NDimensions; i++)
   {
      os << indent.GetNextIndent();
      for (j = 0; j < NDimensions; j++)
      {
         os << m_Matrix[i][j] << " ";
      }
      os << std::endl;
   }

   os << indent << "Offset: " << m_Offset << std::endl;
   os << indent << "Center: " << m_Center << std::endl;
   os << indent << "Translation: " << m_Translation << std::endl;

   os << indent << "Inverse: " << std::endl;
   const InverseMatrixType & invmat = this->GetInverseMatrix();
   for (i = 0; i < NDimensions; i++)
   {
      os << indent.GetNextIndent();
      for (j = 0; j < NDimensions; j++)
      {
         os << invmat[i][j] << " ";
      }
      os << std::endl;
   }
   os << indent << "Singular: " << m_Singular << std::endl;
}

// Constructor with explicit arguments
template<class TScalarType, unsigned int NDimensions>
void
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::SetIdentity( void )
{
   m_Matrix.SetIdentity();
   m_MatrixMTime.Modified();
   m_Offset.Fill( 0.0 );
   m_Translation.Fill( 0.0 );
   m_Center.Fill( 0.0 );
   m_Singular = false;
   m_InverseMatrix.SetIdentity();
   m_InverseMatrixMTime = m_MatrixMTime;
   this->Modified();
}


// Compose with another affine transformation
template<class TScalarType, unsigned int NDimensions>
void
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::Compose(const Self * other, bool pre)
{
   if (pre)
   {
      m_Offset = m_Matrix * other->m_Offset + m_Offset;
      m_Matrix = m_Matrix * other->m_Matrix;
   }
   else
   {
      m_Offset = other->m_Matrix * m_Offset + other->m_Offset;
      m_Matrix = other->m_Matrix * m_Matrix;
   }

   this->ComputeTranslation();
   this->ComputeMatrixParameters();

   m_MatrixMTime.Modified();
   this->Modified();

   return;
}



// Transform a point
template<class TScalarType, unsigned int NDimensions>
typename ESMMatrixOffsetTransformBase<TScalarType,
                                      NDimensions>::PointType
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::TransformPoint(const PointType &point) const
{
   return m_Matrix * point + m_Offset;
}


// Transform a vector
template<class TScalarType, unsigned int NDimensions>
typename ESMMatrixOffsetTransformBase<TScalarType,
                                      NDimensions>::VectorType
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::TransformVector(const VectorType &vect) const
{
   return m_Matrix * vect;
}


// Transform a vnl_vector_fixed
template<class TScalarType, unsigned int NDimensions>
typename ESMMatrixOffsetTransformBase<TScalarType,
                                      NDimensions>::VnlVectorType
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::TransformVector(const VnlVectorType &vect) const
{
   return m_Matrix * vect;
}


// Transform a CovariantVector
template<class TScalarType, unsigned int NDimensions>
typename ESMMatrixOffsetTransformBase<TScalarType,
                                      NDimensions>::CovariantVectorType
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::TransformCovariantVector(const CovariantVectorType &vec) const
{
   CovariantVectorType  result;    // Converted vector

   const InverseMatrixType & invmat = this->GetInverseMatrix();

   for (unsigned int i = 0; i < NDimensions; i++)
   {
      result[i] = NumericTraits<ScalarType>::Zero;
      for (unsigned int j = 0; j < NDimensions; j++)
      {
         result[i] += invmat[j][i]*vec[j]; // Inverse transposed
      }
   }
   return result;
}

// Recompute the inverse matrix (internal)
template<class TScalarType, unsigned int NDimensions>
const typename ESMMatrixOffsetTransformBase<TScalarType,
                                            NDimensions>::InverseMatrixType &
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::GetInverseMatrix( void ) const
{
   // If the transform has been modified we recompute the inverse
   if(m_InverseMatrixMTime != m_MatrixMTime)
   {
      m_Singular = false;
      try
      {
         m_InverseMatrix  = m_Matrix.GetInverse();
      }
      catch(...)
      {
         m_Singular = true;
      }
      m_InverseMatrixMTime = m_MatrixMTime;
   }

   return m_InverseMatrix;
}

template<class TScalarType, unsigned int NDimensions>
typename ESMMatrixOffsetTransformBase<
   TScalarType, NDimensions>::InverseESMTransformBasePointer
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::GetInverseESMTransform() const
{
  Pointer inverse = New();

  this->GetInverseMatrix();
  if(m_Singular)
  {
     return NULL;
  }

  inverse->m_Matrix         = this->GetInverseMatrix();
  inverse->m_InverseMatrix  = m_Matrix;
  inverse->m_Offset         = -(this->GetInverseMatrix() * m_Offset);
  inverse->m_Center         = m_Center;
  inverse->ComputeTranslation();
  inverse->ComputeMatrixParameters();

  return inverse.GetPointer();
}


// Get fixed parameters
template<class TScalarType, unsigned int NDimensions>
void
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::SetFixedParameters( const ParametersType & fp )
{
   this->m_FixedParameters = fp;
   PointType c;
   for ( unsigned int i = 0; i < NDimensions; i++ )
   {
      c[i] = this->m_FixedParameters[i];
   }
   this->SetCenter ( c );
}

/** Get the Fixed Parameters. */
template<class TScalarType, unsigned int NDimensions>
const typename ESMMatrixOffsetTransformBase<TScalarType,
                                            NDimensions>::ParametersType &
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::GetFixedParameters(void) const
{
   this->m_FixedParameters.SetSize ( NDimensions );
   for ( unsigned int i = 0; i < NDimensions; i++ )
   {
      this->m_FixedParameters[i] = this->m_Center[i];
   }
   return this->m_FixedParameters;
}



// Get parameters
template<class TScalarType, unsigned int NDimensions>
const typename ESMMatrixOffsetTransformBase<TScalarType,
                                            NDimensions>::ParametersType &
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::GetParameters( void ) const
{
   // Transfer the linear part
   unsigned int par = 0;

   for(unsigned int row=0; row<NDimensions; row++)
   {
      for(unsigned int col=0; col<NDimensions; col++)
      {
         this->m_Parameters[par] = m_Matrix[row][col];
         ++par;
      }
   }

   // Transfer the constant part
   for(unsigned int i=0; i<NDimensions; i++)
   {
      this->m_Parameters[par] = m_Translation[i];
      ++par;
   }

   return this->m_Parameters;
}


// Set parameters
template<class TScalarType, unsigned int NDimensions>
void
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::SetParameters( const ParametersType & parameters )
{
   if (parameters.Size() <
       (NDimensions * NDimensions + NDimensions))
   {
      itkExceptionMacro
         (<< "Error setting parameters: parameters array size ("
          << parameters.Size() << ") is less than expected "
          << " (NDimensions * NDimensions + NDimensions) "
          << " (" << NDimensions << " * " << NDimensions
          << " + " << NDimensions
          << " = " << NDimensions*NDimensions+NDimensions << ")"
            );
   }

   unsigned int par = 0;

   this->m_Parameters = parameters;

   for(unsigned int row=0; row<NDimensions; row++)
   {
      for(unsigned int col=0; col<NDimensions; col++)
      {
         m_Matrix[row][col] = this->m_Parameters[par];
         ++par;
      }
   }

   // Transfer the constant part
   for(unsigned int i=0; i<NDimensions; i++)
   {
      m_Translation[i] = this->m_Parameters[par];
      ++par;
   }

   m_MatrixMTime.Modified();

   this->ComputeMatrix();  // Not necessary since parameters explicitly define
   //    the matrix
   this->ComputeOffset();

   // Modified is always called since we just have a pointer to the
   // parameters and cannot know if the parameters have changed.
   this->Modified();

}


// Compute the Jacobian in one position
template<class TScalarType, unsigned int NDimensions>
const typename ESMMatrixOffsetTransformBase<TScalarType, NDimensions>::JacobianType &
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::GetJacobian( const PointType & p ) const
{
   this->FillJacobian( this->m_Jacobian, p );

   return this->m_Jacobian;
}

template<class TScalarType, unsigned int NDimensions>
void ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::FillJacobian( JacobianType & jac, const PointType & p ) const
{
   // The Jacobian of the affine transform is composed of
   // subblocks of diagonal matrices, each one of them having
   // a constant value in the diagonal.

   jac.Fill( 0.0 );

   const VectorType v = p - this->GetCenter();

   unsigned int blockOffset = 0;

   for(unsigned int block=0; block < NDimensions; block++)
   {
      for(unsigned int dim=0; dim < NDimensions; dim++ )
      {
         jac( block , blockOffset + dim ) = v[dim];
      }
      blockOffset += NDimensions;
   }

   for(unsigned int dim=0; dim < NDimensions; dim++ )
   {
      jac( dim , blockOffset + dim ) = 1.0;
   }
}

// Computes offset based on center, matrix, and translation variables
template<class TScalarType, unsigned int NDimensions>
void
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::ComputeOffset( void )
{
   const MatrixType & matrix = this->GetMatrix();

   OffsetType offset;
   for(unsigned int i=0; i<NDimensions; i++)
   {
      offset[i] = m_Translation[i] + m_Center[i];
      for(unsigned int j=0; j<NDimensions; j++)
      {
         offset[i] -= matrix[i][j] * m_Center[j];
      }
   }

   m_Offset = offset;
}

// Computes translation based on offset, matrix, and center
template<class TScalarType, unsigned int NDimensions>
void
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::ComputeTranslation( void )
{
   const MatrixType & matrix = this->GetMatrix();

   OffsetType translation;
   for(unsigned int i=0; i<NDimensions; i++)
   {
      translation[i] = m_Offset[i] - m_Center[i];
      for(unsigned int j=0; j<NDimensions; j++)
      {
         translation[i] += matrix[i][j] * m_Center[j];
      }
   }

   m_Translation = translation;
}


// Computes matrix - base class does nothing.  In derived classes is
//    used to convert, for example, versor into a matrix
template<class TScalarType, unsigned int NDimensions>
void
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::ComputeMatrix( void )
{
   // Since parameters explicitly define the matrix in this base class, this
   // function does nothing.  Normally used to compute a matrix when
   // its parameterization (e.g., the class' versor) is modified.
}


// Computes parameters - base class does nothing.  In derived classes is
//    used to convert, for example, matrix into a versor
template<class TScalarType, unsigned int NDimensions>
void
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::ComputeMatrixParameters( void )
{
   // Since parameters explicitly define the matrix in this base class, this
   // function does nothing.  Normally used to update the parameterization
   // of the matrix (e.g., the class' versor) when the matrix is explicitly
   // set.
}

template<class TScalarType, unsigned int NDimensions>
void
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::IncrementalUpdate( const ParametersType & updatevector )
{
   // Store update parameters in homogeneous matrix
   unsigned int par = 0;
   const unsigned int NDim2 = NDimensions*NDimensions;

   for(unsigned int row=0; row<NDimensions; row++)
   {
      this->m_LogIncMatrix(row,NDimensions) = updatevector[NDim2+row];
      for(unsigned int col=0; col<NDimensions; col++)
      {
         this->m_LogIncMatrix(row,col) = updatevector[par];
         ++par;
      }
   }

   // Compute matrix exponential of the homogeneous matrix
   const vnl_matrix<double> incexp = sdtools::GetExponential(m_LogIncMatrix.as_ref());

   // Compose with current transfomation
   // First the translation part: Trans_new = Trans_old + Mat_old * Trans_expup
   for(unsigned int i=0; i<NDimensions; i++)
   {
      for(unsigned int j=0; j<NDimensions; j++)
      {
         this->m_Translation[i] += this->m_Matrix[i][j]*incexp(j,NDimensions);
      }
   }
   // Then the matrix part: Mat_new = Mat_old * Mat_expup
   this->m_Matrix *= incexp.extract(NDimensions, NDimensions, 0, 0);

   // Mae sure the offset is up-to-date
   this->ComputeOffset();

   this->ComputeMatrixParameters();

   this->m_MatrixMTime.Modified();
   this->Modified();
}

template<class TScalarType, unsigned int NDimensions>
const typename ESMMatrixOffsetTransformBase<TScalarType, NDimensions>::JacobianType &
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::GetIncrementalUpdateJacobian(const PointType  & p) const
{
   this->FillJacobian( this->m_IncrementalUpdateJacobian, p );

   return this->m_IncrementalUpdateJacobian;
}

template<class TScalarType, unsigned int NDimensions>
const typename ESMMatrixOffsetTransformBase<TScalarType, NDimensions>::JacobianType &
ESMMatrixOffsetTransformBase<TScalarType, NDimensions>
::GetSpatialJacobian(const PointType  &) const
{
   // For linear transformation, the derivative of
   // M*p+o w.r.t. p is simply M
   this->m_SpatialJacobian = this->GetMatrix().GetVnlMatrix();
   return this->m_SpatialJacobian;
}

} // namespace

#endif


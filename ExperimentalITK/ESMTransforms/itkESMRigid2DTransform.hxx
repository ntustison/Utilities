#ifndef MZ_ESMRigid2DTransform_TXX_
#define MZ_ESMRigid2DTransform_TXX_

#include "itkESMRigid2DTransform.h"

#include <vnl/vnl_math.h>


namespace itk
{

// Constructor with default arguments
template<class TScalarType>
ESMRigid2DTransform<TScalarType>::
ESMRigid2DTransform():
   Superclass(SpaceDimension, ParametersDimension)
{
   m_Angle = NumericTraits< TScalarType >::Zero;

   this->m_IncrementalUpdateJacobian(0,1) = 1.0;
   this->m_IncrementalUpdateJacobian(0,2) = 0.0;
   this->m_IncrementalUpdateJacobian(1,1) = 0.0;
   this->m_IncrementalUpdateJacobian(1,2) = 1.0;
}


// Constructor with arguments
template<class TScalarType>
ESMRigid2DTransform<TScalarType>::
ESMRigid2DTransform( unsigned int spaceDimension,
                     unsigned int parametersDimension):
   Superclass(spaceDimension,parametersDimension)
{
   m_Angle = NumericTraits< TScalarType >::Zero;

   this->m_IncrementalUpdateJacobian(0,1) = 1.0;
   this->m_IncrementalUpdateJacobian(0,2) = 0.0;
   this->m_IncrementalUpdateJacobian(1,1) = 0.0;
   this->m_IncrementalUpdateJacobian(1,2) = 1.0;
}


// Destructor
template<class TScalarType>
ESMRigid2DTransform<TScalarType>::
~ESMRigid2DTransform()
{
}


// Print self
template<class TScalarType>
void
ESMRigid2DTransform<TScalarType>::
PrintSelf(std::ostream &os, Indent indent) const
{
   Superclass::PrintSelf(os,indent);
   os << indent << "Angle       = " << m_Angle        << std::endl;
}


// Set the rotation matrix
template<class TScalarType>
void
ESMRigid2DTransform<TScalarType>::
SetMatrix(const MatrixType & matrix )
{
   itkDebugMacro("setting  m_Matrix  to " << matrix );
   // The matrix must be orthogonal otherwise it is not
   // representing a valid rotation in 2D space
   typename MatrixType::InternalMatrixType test =
      matrix.GetVnlMatrix() * matrix.GetTranspose();

   const double tolerance = 1e-10;
   if( !test.is_identity( tolerance ) )
   {
      itk::ExceptionObject ex(__FILE__,__LINE__,"Attempt to set a Non-Orthogonal matrix",ITK_LOCATION);
      throw ex;
   }

   this->m_Matrix = matrix;
   this->ComputeOffset();
   this->ComputeMatrixParameters();
   this->Modified();

}


/** Compute the Angle from the Rotation Matrix */
template <class TScalarType>
void
ESMRigid2DTransform<TScalarType>
::ComputeMatrixParameters( void )
{
   const MatrixType & mat = this->GetMatrix();

   if( vcl_abs(mat[0][0]-mat[1][1])>0.000001 ||
       vcl_abs(mat[0][1]+mat[1][0])>0.000001 ||
       vcl_abs(vnl_math_sqr(mat[0][0])+vnl_math_sqr(mat[1][0])-1.0)>0.000001 )
   {
      itkWarningMacro("Bad Rotation Matrix " << this->GetMatrix() );
   }

   m_Angle = atan2(mat[1][0],mat[0][0]);
}


// Compose with a translation
template<class TScalarType>
void
ESMRigid2DTransform<TScalarType>::
Translate(const OffsetType &offset, bool)
{
   VectorType newOffset = this->GetOffset();
   newOffset += offset;
   this->SetOffset(newOffset);
   this->ComputeTranslation();
}

// Return an inverse of this transform
template<class TScalarType>
typename ESMRigid2DTransform<TScalarType>::InverseESMTransformBasePointer
ESMRigid2DTransform<TScalarType>
::GetInverseESMTransform() const
{
   Pointer inverse = New();

   inverse->SetCenter( this->GetCenter() );  // inverse have the same center
   inverse->SetAngle( -this->GetAngle() );
   inverse->SetTranslation( -( this->GetInverseMatrix() * this->GetTranslation() ) );

   return inverse.GetPointer();
}

// Create and return a clone of the transformation
template<class TScalarType>
void
ESMRigid2DTransform<TScalarType>::
CloneTo( Pointer & result ) const
{
   result = New();
   result->SetCenter( this->GetCenter() );
   result->SetAngle( this->GetAngle() );
   result->SetTranslation( this->GetTranslation() );
}


// Reset the transform to an identity transform
template<class TScalarType >
void
ESMRigid2DTransform< TScalarType >::
SetIdentity( void )
{
   this->Superclass::SetIdentity();
   m_Angle = NumericTraits< TScalarType >::Zero;
}

// Set the angle of rotation
template <class TScalarType>
void
ESMRigid2DTransform<TScalarType>
::SetAngle(TScalarType angle)
{
   m_Angle = angle;
   this->ComputeMatrix();
   this->ComputeOffset();
   this->Modified();
}


// Set the angle of rotation
template <class TScalarType>
void
ESMRigid2DTransform<TScalarType>
::SetAngleInDegrees(TScalarType angle)
{
   const TScalarType angleInRadians = angle * vcl_atan(1.0) / 45.0;
   this->SetAngle( angleInRadians );
}

// Compute the matrix from the angle
template <class TScalarType>
void
ESMRigid2DTransform<TScalarType>
::ComputeMatrix( void )
{
   const double ca = vcl_cos(m_Angle );
   const double sa = vcl_sin(m_Angle );

   this->m_Matrix[0][0]= ca;
   this->m_Matrix[0][1]=-sa;
   this->m_Matrix[1][0]= sa;
   this->m_Matrix[1][1]= ca;
}

// Set Parameters
template <class TScalarType>
void
ESMRigid2DTransform<TScalarType>::
SetParameters( const ParametersType & parameters )
{
   itkDebugMacro( << "Setting parameters " << parameters );

   // Set angle
   this->SetVarAngle( parameters[0] );

   // Set translation
   for(unsigned int i=0; i < SpaceDimension; i++)
   {
      this->m_Translation[i] = parameters[i+1];
   }

   // Update matrix and offset
   this->ComputeMatrix();
   this->ComputeOffset();

   // Modified is always called since we just have a pointer to the
   // parameters and cannot know if the parameters have changed.
   this->Modified();

   itkDebugMacro(<<"After setting parameters ");
}

// Get Parameters
template <class TScalarType>
const typename ESMRigid2DTransform<TScalarType>::ParametersType &
ESMRigid2DTransform<TScalarType>::
GetParameters( void ) const
{
   itkDebugMacro( << "Getting parameters ");

   // Get the angle
   this->m_Parameters[0] = this->GetAngle();

   // Get the translation
   for(unsigned int i=0; i < SpaceDimension; i++)
   {
      this->m_Parameters[i+1] = this->GetTranslation()[i];
   }

   itkDebugMacro(<<"After getting parameters " << this->m_Parameters );

   return this->m_Parameters;
}

// Compute transformation Jacobian
template<class TScalarType>
const typename ESMRigid2DTransform<TScalarType>::JacobianType &
ESMRigid2DTransform<TScalarType>::
GetJacobian( const PointType & p ) const
{

   const double ca = vcl_cos(this->GetAngle() );
   const double sa = vcl_sin(this->GetAngle() );

   this->m_Jacobian.Fill(0.0);

   const double cx = this->GetCenter()[0];
   const double cy = this->GetCenter()[1];

   // derivatives with respect to the angle
   this->m_Jacobian[0][0] = -sa * ( p[0] - cx ) - ca * ( p[1] - cy );
   this->m_Jacobian[1][0] =  ca * ( p[0] - cx ) - sa * ( p[1] - cy );

   // compute derivatives for the translation part
   unsigned int blockOffset = 1;
   for(unsigned int dim=0; dim < SpaceDimension; dim++ )
   {
      this->m_Jacobian[ dim ][ blockOffset + dim ] = 1.0;
   }

   return this->m_Jacobian;

}

template<class TScalarType>
void
ESMRigid2DTransform<TScalarType>::
IncrementalUpdate( const ParametersType & updatevector )
{
   // We have a closed-form solution for the matrix exponential in this case

   // Compose with current transfomation
   // First the translation part: Trans_new = Trans_old + Mat_old * Trans_expup
   VectorType inc_trans;

   const TScalarType absalpha = vcl_abs( updatevector[0] );
   if ( absalpha > 1e-4 )
   {
      // Compute sinc function directly
      const TScalarType sincalpha = vcl_sin( updatevector[0] ) / updatevector[0];
      const TScalarType cosalpham1overalpha = (vcl_cos( updatevector[0] )-1.0) / updatevector[0];

      inc_trans[0] = sincalpha*updatevector[1] + cosalpham1overalpha*updatevector[2];
      inc_trans[1] = -cosalpham1overalpha*updatevector[1] + sincalpha*updatevector[2];
   }
   else
   {
      // Use Taylor expansion up to order 2
      const TScalarType sincalpha = 1.0 - vnl_math_sqr(updatevector[0])/6.0;
      const TScalarType cosalpham1overalpha = -updatevector[0]/2.0;

      inc_trans[0] = sincalpha*updatevector[1] + cosalpham1overalpha*updatevector[2];
      inc_trans[1] = -cosalpham1overalpha*updatevector[1] + sincalpha*updatevector[2];
   }

   this->m_Translation += this->m_Matrix*inc_trans;

   // Set Angle takes care of updating the matrix and offset
   // and calls modified
   this->SetAngle( this->m_Angle + updatevector[0] );
}

template<class TScalarType>
const typename ESMRigid2DTransform<TScalarType>::JacobianType &
ESMRigid2DTransform<TScalarType>::
GetIncrementalUpdateJacobian( const PointType & p ) const
{
   const PointType & center = this->GetCenter();

   this->m_IncrementalUpdateJacobian(0,0) = center[1] - p[1];
   this->m_IncrementalUpdateJacobian(1,0) = p[0] - center[0];

   return this->m_IncrementalUpdateJacobian;
}

} // namespace

#endif

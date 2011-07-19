/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDecomposeTensorImageFilter.hxx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:16:52 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkDecomposeTensorImageFilter_hxx
#define _itkDecomposeTensorImageFilter_hxx

#include "itkDecomposeTensorImageFilter.h"

#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"

#include "vnl/algo/vnl_qr.h"
#include "vnl/algo/vnl_real_eigensystem.h"
#include "vnl/algo/vnl_svd.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "vnl/vnl_matrix.h"
#include "vxl/vcl/vcl_complex.h"

#include <algorithm>
#include <vector>


namespace itk
{
template <typename TInputImage, typename TRealType, typename TOutputImage>
DecomposeTensorImageFilter<TInputImage, TRealType, TOutputImage>
::DecomposeTensorImageFilter()
{
  this->m_SymmetricTensors = false;
  this->SetNumberOfRequiredInputs( 1 );

  this->m_WhichDecomposition = 0;
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
void
DecomposeTensorImageFilter<TInputImage, TRealType, TOutputImage>
::GenerateData()
{
  /**
   * Set up outputs and allocate
   */
  switch ( this->m_WhichDecomposition )
    {
    case 0: case 1: case 2: case 3:
      this->SetNumberOfRequiredOutputs( 2 );
      break;      
    case 4:
      this->SetNumberOfRequiredOutputs( 3 );
      break;
    }   
  for ( unsigned int i = 0; i < this->GetNumberOfRequiredOutputs(); i++ )
    {
    this->SetNthOutput( i, ( OutputImageType::New() ).GetPointer() );
    this->GetOutput( i )->SetOrigin( this->GetInput()->GetOrigin() );
    this->GetOutput( i )->SetSpacing( this->GetInput()->GetSpacing() );
    this->GetOutput( i )->SetRegions( this->GetInput()->GetRequestedRegion() );  
    this->GetOutput( i )->Allocate();
    } 

  ImageRegionConstIteratorWithIndex<InputImageType> It
    ( this->GetInput(), this->GetInput()->GetRequestedRegion() );

  OutputMatrixType D;
  OutputMatrixType Q;
  OutputMatrixType R;
  OutputMatrixType S;
  OutputMatrixType U;
  OutputMatrixType W;
  OutputMatrixType V;

  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    InputMatrixType M = It.Get();
    switch ( this->m_WhichDecomposition )
      {
      case 0:  // EigenDecomposition
        this->EigenDecomposition( M, D, V );
        this->GetOutput( 0 )->SetPixel( It.GetIndex(), D );  // Eigenvalues
        this->GetOutput( 1 )->SetPixel( It.GetIndex(), V );  // Eigenvectors (columnwise)
        break;      
      case 1:  // Right Polar Decomposition
        this->RightPolarDecomposition( M, R, S );
        this->GetOutput( 0 )->SetPixel( It.GetIndex(), R );  // Rotation
        this->GetOutput( 1 )->SetPixel( It.GetIndex(), S );  // Stretch
        break;
      case 2:  // Left Polar Decomposition
        this->LeftPolarDecomposition( M, S, R );
        this->GetOutput( 0 )->SetPixel( It.GetIndex(), S );  // Stretch
        this->GetOutput( 1 )->SetPixel( It.GetIndex(), R );  // Rotation
        break;
      case 3:  // QR Decomposition
        this->QRDecomposition( M, Q, R );
        this->GetOutput( 0 )->SetPixel( It.GetIndex(), Q );  // Orthogonal matrix
        this->GetOutput( 1 )->SetPixel( It.GetIndex(), R );  // Upper triangular matrix
        break;      
      case 4:  // SVD Decomposition
        this->SVDDecomposition( M, U, W, V );
        this->GetOutput( 0 )->SetPixel( It.GetIndex(), U );  // Unitary matrix - U*U^* = U^*U = I
        this->GetOutput( 1 )->SetPixel( It.GetIndex(), W );  // Non-negative diagonal matrix
        this->GetOutput( 2 )->SetPixel( It.GetIndex(), V );  // Unitary matrix
        break;      
      }   
    } 
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
void
DecomposeTensorImageFilter<TInputImage, TRealType, TOutputImage>
::EigenDecomposition( InputMatrixType &M, OutputMatrixType &D, OutputMatrixType &V )
{
  D.SetSize( RowDimensions, RowDimensions );
  V.SetSize( RowDimensions, RowDimensions );

  D.Fill( 0.0 );
  if ( this->m_SymmetricTensors )
    {
    vnl_symmetric_eigensystem<RealType> eig( M.GetVnlMatrix() );

    for ( unsigned int j = 0; j < ColumnDimensions; j++ )
      {
      for ( unsigned int i = 0; i < RowDimensions; i++ )
        {
        V[i][j] = (eig.V).get( i, j );
        if ( i == j )
          {
          D[i][j] = eig.D(j);
          }  
        }  
      }
    }
  else
    {
    vnl_matrix<double> v( RowDimensions, ColumnDimensions );
    for ( unsigned int j = 0; j < ColumnDimensions; j++ )
      {
      for ( unsigned int i = 0; i < RowDimensions; i++ )
        {
        v.put( i, j, M[i][j] );
        }  
      }
    vnl_real_eigensystem eig( v );
    for ( unsigned int j = 0; j < ColumnDimensions; j++ )
      {
      for ( unsigned int i = 0; i < RowDimensions; i++ )
        {
        V[i][j] = static_cast<RealType>( (eig.Vreal).get( i, j ) );
        if ( i == j )
          {
          D[i][j] = static_cast<RealType>( vcl_real( eig.D(j) ) );
          }  
        }  
      }
    }
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
void
DecomposeTensorImageFilter<TInputImage, TRealType, TOutputImage>
::RightPolarDecomposition( InputMatrixType &M, OutputMatrixType &R, OutputMatrixType &S )
{
  OutputMatrixType U;
  OutputMatrixType W;
  OutputMatrixType V;
  this->SVDDecomposition( M, U, W, V );

  R = U * V.GetTranspose();
  S = V * W * V.GetTranspose();
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
void
DecomposeTensorImageFilter<TInputImage, TRealType, TOutputImage>
::LeftPolarDecomposition( InputMatrixType &M, OutputMatrixType &S, OutputMatrixType &R )
{
  OutputMatrixType U;
  OutputMatrixType W;
  OutputMatrixType V;
  this->SVDDecomposition( M, U, W, V );

  R = U * V.GetTranspose();
  S = U * W * U.GetTranspose();
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
void
DecomposeTensorImageFilter<TInputImage, TRealType, TOutputImage>
::QRDecomposition( InputMatrixType &M, OutputMatrixType &Q, OutputMatrixType &R )
{
  Q.SetSize( RowDimensions, ColumnDimensions );
  R.SetSize( ColumnDimensions, ColumnDimensions );

  vnl_qr<RealType> qr( M.GetVnlMatrix() );
  
  for ( unsigned int i = 0; i < RowDimensions; i++ )
    {
    for ( unsigned int j = 0; j < ColumnDimensions; j++ )
      {
      Q[i][j] = qr.Q()(i, j);  
      }
    }
  for ( unsigned int i = 0; i < ColumnDimensions; i++ )
    {
    for ( unsigned int j = 0; j < ColumnDimensions; j++ )
      {
      R[i][j] = qr.R()(i, j);  
      }
    }
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
void
DecomposeTensorImageFilter<TInputImage, TRealType, TOutputImage>
::SVDDecomposition( InputMatrixType &M, OutputMatrixType &U, OutputMatrixType &W, OutputMatrixType &V )
{
  U.SetSize( RowDimensions, RowDimensions );
  V.SetSize( ColumnDimensions, ColumnDimensions );
  W.SetSize( RowDimensions, ColumnDimensions );

  vnl_svd<RealType> svd( M.GetVnlMatrix() );
  
  for ( unsigned int i = 0; i < RowDimensions; i++ )
    {
    for ( unsigned int j = 0; j < RowDimensions; j++ )
      {
      U[i][j] = svd.U(i, j);  
      }
    }
  for ( unsigned int i = 0; i < ColumnDimensions; i++ )
    {
    for ( unsigned int j = 0; j < ColumnDimensions; j++ )
      {
      V[i][j] = svd.V(i, j);  
      }
    }

  W.Fill( 0.0 );
  unsigned int minDimensions = ColumnDimensions;
  if ( static_cast<unsigned int>( RowDimensions ) < 
       static_cast<unsigned int>( ColumnDimensions ) )
    {
    minDimensions = RowDimensions;
    }
  for ( unsigned int i = 0; i < minDimensions; i++ )
    {
    W[i][i] = svd.W(i, i);
    }  
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename DecomposeTensorImageFilter<TInputImage, TRealType, TOutputImage>::RealImageType *
DecomposeTensorImageFilter<TInputImage, TRealType, TOutputImage>
::GetEigenValueImage( unsigned int i )
{
  /**
   * Eigenvalues are sorted in increasing order ( i = 0 
   * is smallest eigenvalue ) --- ordering only holds for
   * symmetric matrices.
   */

  if ( this->m_WhichDecomposition != 0 || i >= ColumnDimensions )
    {
    return NULL;
    }

  this->m_EigenValueImage = RealImageType::New();
  this->m_EigenValueImage->SetOrigin( this->GetInput()->GetOrigin() );
  this->m_EigenValueImage->SetSpacing( this->GetInput()->GetSpacing() );
  this->m_EigenValueImage->SetRegions( this->GetInput()->GetRequestedRegion() );  
  this->m_EigenValueImage->Allocate();

  ImageRegionIterator<RealImageType> ItE( this->m_EigenValueImage, 
    this->m_EigenValueImage->GetLargestPossibleRegion() );
  ImageRegionIterator<OutputImageType> ItD( this->GetOutput( 0 ), 
    this->GetOutput( 0 )->GetLargestPossibleRegion() );
 
  ItE.GoToBegin();
  ItD.GoToBegin();
  while ( !ItE.IsAtEnd() && !ItD.IsAtEnd() )
    {
    ItE.Set( ItD.Get()[i][i] );
    ++ItE;
    ++ItD;
    } 

  return this->m_EigenValueImage.GetPointer();
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
typename DecomposeTensorImageFilter<TInputImage, TRealType, TOutputImage>::VectorImageType *
DecomposeTensorImageFilter<TInputImage, TRealType, TOutputImage>
::GetEigenVectorImage( unsigned int i )
{
  /**
   * Eigenvectors are sorted in increasing order of corresponding
   * eigenvalue ( i = 0 is smallest eigenvalue ) --- ordering only 
   * holds for symmetric matrices.
   */

  if ( this->m_WhichDecomposition != 0 || i >= ColumnDimensions )
    {
    return NULL;
    }

  this->m_EigenVectorImage = VectorImageType::New();
  this->m_EigenVectorImage->SetOrigin( this->GetInput()->GetOrigin() );
  this->m_EigenVectorImage->SetSpacing( this->GetInput()->GetSpacing() );
  this->m_EigenVectorImage->SetRegions( this->GetInput()->GetRequestedRegion() );  
  this->m_EigenVectorImage->Allocate();

  ImageRegionIterator<RealImageType> ItE( this->m_EigenVectorImage, 
    this->m_EigenVectorImage->GetLargestPossibleRegion() );
  ImageRegionIterator<OutputImageType> ItV( this->GetOutput( 1 ), 
    this->GetOutput( 1 )->GetLargestPossibleRegion() );
 
  ItE.GoToBegin();
  ItV.GoToBegin();
  while ( !ItE.IsAtEnd() && !ItV.IsAtEnd() )
    {
    VectorType V;
    for ( unsigned int j = 0; j < RowDimensions; j++ )
      {
      V[j] = ItV.Get()[j][i];
      }  
    ItE.Set( V );
    ++ItE;
    ++ItV;
    } 

  return this->m_EigenVectorImage.GetPointer();
}

template <typename TInputImage, typename TRealType, typename TOutputImage>
void
DecomposeTensorImageFilter<TInputImage, TRealType, TOutputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "m_SymmetricTensors = " << this->m_SymmetricTensors << std::endl;
  os << indent << "m_WhichDecomposition = " << this->m_WhichDecomposition << std::endl;
}
  
} // end namespace itk

#endif

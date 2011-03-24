/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSparseMatrixEigenAnalysis.txx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:20:04 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkSparseMatrixEigenAnalysis_txx
#define __itkSparseMatrixEigenAnalysis_txx

#include "itkSparseMatrixEigenAnalysis.h"

#include "arpackf.h"

namespace itk
{

template <class TRealType, class TSparseMatrix, class TEigenVector>
SparseMatrixEigenAnalysis<TRealType, TSparseMatrix, TEigenVector >
::SparseMatrixEigenAnalysis()
{
  this->m_IsSymmetric = false;
  this->m_SolveForSmallestEigenValues = true;
  this->m_NumberOfEigenPairs = 2;

  this->m_Tolerance = vnl_math::eps;
  this->m_MaximumNumberOfIterations = 300;
  this->m_NumberOfLanczosVectors = 0;

  this->m_Matrix1 = NULL;
  this->m_Matrix2 = NULL;
}  

/**
 * Destructor
 */
template <class TRealType, class TSparseMatrix, class TEigenVector>
SparseMatrixEigenAnalysis<TRealType, TSparseMatrix, TEigenVector >
::~SparseMatrixEigenAnalysis()
{
}


template <class TRealType, class TSparseMatrix, class TEigenVector>
void
SparseMatrixEigenAnalysis<TRealType, TSparseMatrix, TEigenVector >
::GenerateData()
{
  if ( typeid( RealType ) == typeid( float ) )
    {
    this->CalculateEigenPairsFloat();
    }
  else
    {
    this->CalculateEigenPairsDouble();
    }
}

template <class TRealType, class TSparseMatrix, class TEigenVector>
void
SparseMatrixEigenAnalysis<TRealType, TSparseMatrix, TEigenVector >
::CalculateEigenPairsFloat()
{
  SparseMatrixType M = *this->m_Matrix1;
 
  if ( M.rows() != M.columns() )
    {
    itkExceptionMacro( "The number of rows must equal the number of columns." );
    }

  /**
   * ARPACK work variables
   */
  int N = M.rows();
  int ldv = N;
  char which[2];
  if ( this->m_SolveForSmallestEigenValues )
    {
    which[0] = 'S';
    which[1] = 'M'; 
    }
  else
    {
    which[0] = 'L';
    which[1] = 'M'; 
    }
  char bmat[] = "I";  

  if ( this->m_NumberOfEigenPairs <= 0 )
    {
    this->SetNumberOfEigenPairs( 1 );
    itkWarningMacro( "Setting this->m_NumberOfEigenPairs to " << this->m_NumberOfEigenPairs );
    }
  else if ( this->m_NumberOfEigenPairs >= N )
    {
    this->SetNumberOfEigenPairs( N-1 ); 
    itkWarningMacro( "Setting this->m_NumberOfEigenPairs to " << this->m_NumberOfEigenPairs );
    }  
  int nev = this->m_NumberOfEigenPairs;

  if ( this->m_NumberOfLanczosVectors == 0 )
    {
    this->m_NumberOfLanczosVectors = 2 * this->m_NumberOfEigenPairs;
    }
  if ( this->m_NumberOfLanczosVectors <= this->m_NumberOfEigenPairs+1 )
    {
    this->SetNumberOfLanczosVectors( this->m_NumberOfEigenPairs + 2 );
    itkWarningMacro( "Setting this->m_NumberOfLanczosVectors to " << this->m_NumberOfLanczosVectors );
    }
  else if ( this->m_NumberOfLanczosVectors >= N )
    {
    this->SetNumberOfLanczosVectors( N-1 ); 
    itkWarningMacro( "Setting this->m_NumberOfLanczosVectors to " << this->m_NumberOfLanczosVectors );
    }  
  int ncv = this->m_NumberOfLanczosVectors;

  int ido = 0;
  int info = 0;
  int iparam[11];
  iparam[0] = 1;
  iparam[1] = 0;
  iparam[2] = this->m_MaximumNumberOfIterations;
  iparam[3] = 0;
  iparam[4] = 0;
  iparam[5] = 0;
  iparam[6] = 1;
  iparam[7] = 0;
  iparam[8] = 0;
  iparam[9] = 0;
  iparam[10] = 0;

  int ipntr[11];
  for ( unsigned int i = 0; i < 11; i++ )
    {
    ipntr[i] = 0;
    } 

  int rvec = true;               
  int *select = new int[ncv];    
  char howmny[] = "A";
  int ldz = nev;

  float *resid = new float[N];
  float *v = new float[ncv*N];
  float *workd = new float[3*N];
  for ( unsigned int i = 0; i < 3*N; i++ )
    {
    workd[i] = 0;
    }
  int lworkl;
  if ( this->m_IsSymmetric )
    {
    lworkl = ncv * ( ncv + 8 );
    }
  else
    {
    lworkl = 3 * ncv * ( ncv + 2 );  
    }  
  float *workl = new float[lworkl];
  for ( unsigned int i = 0; i < lworkl; i++ )
    {
    workl[i] = 0;
    }
  float tol = static_cast<float>( this->m_Tolerance );

  unsigned int iter = 0;
  while ( iter++ < this->m_MaximumNumberOfIterations )
    {
    if ( this->m_IsSymmetric )
      { 
      ssaupd_( &ido, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, 
               iparam, ipntr, workd, workl, &lworkl, &info );
      }
    else
      {
      snaupd_( &ido, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, 
               iparam, ipntr, workd, workl, &lworkl, &info );
      }   
    if ( ido == 1 || ido == -1 ) 
      {
      VectorType result( M.rows() );
      VectorType rhs( M.columns() );

      for ( unsigned int i = 0; i < M.columns(); i++ )
        {
        rhs[i] = workd[ipntr[0]-1+i];
        } 
      M.mult( rhs, result ); 

      for ( unsigned int i = 0; i < M.rows(); i++ )
        {
        workd[ipntr[1]-1+i] = result[i];
        } 
      }
    else
      {
      break;
      } 
    }
  if ( info < 0 )
    {
    if ( this->m_IsSymmetric )
      {
      itkWarningMacro( "Error with dsaupd_ routine (info = " << info << ")." );
      }
    else
      {
      itkWarningMacro( "Error with dnaupd_ routine (info = " << info << ")." );
      }   
    }

  float *dr = new float[nev+1];
  float *di = new float[nev+1];
  float *z = new float[N*(nev+1)];
  float sigma = 0;
  float workev[3*ncv];

  if ( this->m_IsSymmetric )
    { 
    sseupd_( &rvec, howmny, select, dr, z, &ldv, &sigma, bmat, &N, which, 
             &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, 
             &lworkl, &info );
    }
  else
    {
    sneupd_( &rvec, howmny, select, dr, di, z, &ldz, &sigma, &sigma, workev, bmat, &N, which, 
             &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, 
             &lworkl, &info );
    }   
  if ( info < 0 )
    {
    if ( this->m_IsSymmetric )
      {
      itkWarningMacro( "Error with dseupd_ routine (info = " << info << ")." );
      }
    else
      {
      itkWarningMacro( "Error with dneupd_ routine (info = " << info << ")." );
      }   
    }
  if( iparam[4] < nev )
    {
    itkWarningMacro( "Only " << iparam[4] << " of the " << nev << " requested eigenvalues converged." );
    }

  /**
   * Now store the resulting eigenvalues/eigenvectors for subsequent retrieval
   */
  this->m_EigenValues.set_size( nev );
  this->m_EigenVectors.clear();

  for ( unsigned int i = 0; i < nev; i++ )
    {
    this->m_EigenValues[i] = dr[i];   
    VectorType eigenVector( N );
    for ( unsigned int j = 0; j < N; j++ )
      { 
      eigenVector[j] = z[i*N + j];
      }
    this->m_EigenVectors.push_back( eigenVector );    
    } 

  delete [] dr;
  delete [] di;
  delete [] z;
  delete [] v;
  delete [] resid;
  delete [] workd;
  delete [] workl;
  delete [] select;
}

template <class TRealType, class TSparseMatrix, class TEigenVector>
void
SparseMatrixEigenAnalysis<TRealType, TSparseMatrix, TEigenVector >
::CalculateEigenPairsDouble()
{

  SparseMatrixType M = *this->m_Matrix1;
 
  if ( M.rows() != M.columns() )
    {
    itkExceptionMacro( "The number of rows must equal the number of columns." );
    }

  /**
   * ARPACK work variables
   */
  int N = M.rows();
  int ldv = N;
  char which[2];
  if ( this->m_SolveForSmallestEigenValues )
    {
    which[0] = 'S';
    which[1] = 'M'; 
    }
  else
    {
    which[0] = 'L';
    which[1] = 'M'; 
    }
  char bmat[] = "I";  

  if ( this->m_NumberOfEigenPairs <= 0 )
    {
    this->SetNumberOfEigenPairs( 1 );
    itkWarningMacro( "Setting this->m_NumberOfEigenPairs to " << this->m_NumberOfEigenPairs );
    }
  else if ( this->m_NumberOfEigenPairs >= N )
    {
    this->SetNumberOfEigenPairs( N-1 ); 
    itkWarningMacro( "Setting this->m_NumberOfEigenPairs to " << this->m_NumberOfEigenPairs );
    }  
  int nev = this->m_NumberOfEigenPairs;

  if ( this->m_NumberOfLanczosVectors == 0 )
    {
    this->m_NumberOfLanczosVectors = 2 * this->m_NumberOfEigenPairs;
    }
  if ( this->m_NumberOfLanczosVectors <= this->m_NumberOfEigenPairs+1 )
    {
    this->SetNumberOfLanczosVectors( this->m_NumberOfEigenPairs + 2 );
    itkWarningMacro( "Setting this->m_NumberOfLanczosVectors to " << this->m_NumberOfLanczosVectors );
    }
  else if ( this->m_NumberOfLanczosVectors >= N )
    {
    this->SetNumberOfLanczosVectors( N-1 ); 
    itkWarningMacro( "Setting this->m_NumberOfLanczosVectors to " << this->m_NumberOfLanczosVectors );
    }  
  int ncv = this->m_NumberOfLanczosVectors;

  int ido = 0;
  int info = 0;
  int iparam[11];
  iparam[0] = 1;
  iparam[1] = 0;
  iparam[2] = this->m_MaximumNumberOfIterations;
  iparam[3] = 0;
  iparam[4] = 0;
  iparam[5] = 0;
  iparam[6] = 1;
  iparam[7] = 0;
  iparam[8] = 0;
  iparam[9] = 0;
  iparam[10] = 0;

  int ipntr[11];
  for ( unsigned int i = 0; i < 11; i++ )
    {
    ipntr[i] = 0;
    } 

  int rvec = true;               
  int *select = new int[ncv];    
  char howmny[] = "A";
  int ldz = nev;

  double *resid = new double[N];
  double *v = new double[ncv*N];
  double *workd = new double[3*N];
  for ( unsigned int i = 0; i < 3*N; i++ )
    {
    workd[i] = 0;
    }
  int lworkl;
  if ( this->m_IsSymmetric )
    {
    lworkl = ncv * ( ncv + 8 );
    }
  else
    {
    lworkl = 3 * ncv * ( ncv + 2 );  
    }  
  double *workl = new double[lworkl];
  for ( unsigned int i = 0; i < lworkl; i++ )
    {
    workl[i] = 0;
    }
  double tol = static_cast<double>( this->m_Tolerance );

  unsigned int iter = 0;
  while ( iter++ < this->m_MaximumNumberOfIterations )
    {
    if ( this->m_IsSymmetric )
      { 
      dsaupd_( &ido, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, 
               iparam, ipntr, workd, workl, &lworkl, &info );
      }
    else
      {
      dnaupd_( &ido, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, 
               iparam, ipntr, workd, workl, &lworkl, &info );
      }   
    if ( ido == 1 || ido == -1 ) 
      {
      VectorType result( M.rows() );
      VectorType rhs( M.columns() );

      for ( unsigned int i = 0; i < M.columns(); i++ )
        {
        rhs[i] = workd[ipntr[0]-1+i];
        } 
      M.mult( rhs, result ); 

      for ( unsigned int i = 0; i < M.rows(); i++ )
        {
        workd[ipntr[1]-1+i] = result[i];
        } 
      }
    else
      {
      break;
      } 
    }
  if ( info < 0 )
    {
    if ( this->m_IsSymmetric )
      {
      itkWarningMacro( "Error with dsaupd_ routine (info = " << info << ")." );
      }
    else
      {
      itkWarningMacro( "Error with dnaupd_ routine (info = " << info << ")." );
      }   
    }

  double *dr = new double[nev+1];
  double *di = new double[nev+1];
  double *z = new double[N*(nev+1)];
  double sigma = 0;
  double workev[3*ncv];

  if ( this->m_IsSymmetric )
    { 
    dseupd_( &rvec, howmny, select, dr, z, &ldv, &sigma, bmat, &N, which, 
             &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, 
             &lworkl, &info );
    }
  else
    {
    dneupd_( &rvec, howmny, select, dr, di, z, &ldz, &sigma, &sigma, workev, bmat, &N, which, 
             &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, 
             &lworkl, &info );
    }   
  if ( info < 0 )
    {
    if ( this->m_IsSymmetric )
      {
      itkWarningMacro( "Error with dseupd_ routine (info = " << info << ")." );
      }
    else
      {
      itkWarningMacro( "Error with dneupd_ routine (info = " << info << ")." );
      }   
    }
  if( iparam[4] < nev )
    {
    itkWarningMacro( "Only " << iparam[4] << " of the " << nev << " requested eigenvalues converged." );
    }

  /**
   * Now store the resulting eigenvalues/eigenvectors for subsequent retrieval
   */
  this->m_EigenValues.set_size( nev );
  this->m_EigenVectors.clear();

  for ( unsigned int i = 0; i < nev; i++ )
    {
    this->m_EigenValues[i] = dr[i];   
    VectorType eigenVector( N );
    for ( unsigned int j = 0; j < N; j++ )
      { 
      eigenVector[j] = z[i*N + j];
      }
    this->m_EigenVectors.push_back( eigenVector );    
    } 

  delete [] dr;
  delete [] di;
  delete [] z;
  delete [] v;
  delete [] resid;
  delete [] workd;
  delete [] workl;
  delete [] select;
}

template <class TRealType, class TSparseMatrix, class TEigenVector>
void
SparseMatrixEigenAnalysis<TRealType, TSparseMatrix, TEigenVector >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "m_IsSymmetric = " << this->m_IsSymmetric << std::endl;
  os << indent << "m_SolveForSmallestEigenValues = " << this->m_IsSymmetric << std::endl;
  os << indent << "m_Tolerance = " << this->m_Tolerance << std::endl;
  os << indent << "m_MaximumNumberOfIterations = " << this->m_MaximumNumberOfIterations << std::endl;
  os << indent << "m_NumberOfEigenPairs = " << this->m_NumberOfEigenPairs << std::endl;
  os << indent << "m_NumberOfLanczosVectors = " << this->m_NumberOfLanczosVectors << std::endl;

}



}  // end namespace itk
 
#endif
 

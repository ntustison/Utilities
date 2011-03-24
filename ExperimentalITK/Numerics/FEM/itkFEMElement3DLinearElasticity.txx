/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFEMElement3DLinearElasticity.txx,v $
  Language:  C++
  Date:      $Date: 2008/10/18 00:22:48 $
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkFEMElement3DLinearElasticity_txx
#define __itkFEMElement3DLinearElasticity_txx

#include "itkFEMElement3DLinearElasticity.h"

namespace itk {
namespace fem {




template<class TBaseClass>
Element3DLinearElasticity<TBaseClass>
::Element3DLinearElasticity() : Superclass(), m_mat(0) {}




//////////////////////////////////////////////////////////////////////////
/*
 * Methods related to the physics of the problem.
 */

template<class TBaseClass>
void Element3DLinearElasticity<TBaseClass>
::GetStrainDisplacementMatrix(MatrixType& B, const MatrixType& shapeDgl) const
{
  unsigned int p;
  unsigned int Nn=3*this->GetNumberOfNodes();
  B.set_size(6,Nn);

  // Initialize the B matrix to zero - subsequently, only the nonzero
  // terms will be filled in
  B.fill(0.0);

  // Copy the shape function derivatives wrt global coordinates
  // in right position in B matrix.

  for (unsigned int i=0; i<Nn; i++) {  
    p = i / 3;
    
    switch(i % 3) 
    {
    
      case 0:  /** Columns 1, 4, 7, ..., 22 */
        B[0][i] = shapeDgl[0][p];
        B[3][i] = shapeDgl[1][p];
        B[5][i] = shapeDgl[2][p];
        break;
        
        case 1:  /** Columns 2, 5, 8, ..., 23 */
        B[1][i] = shapeDgl[1][p];
        B[3][i] = shapeDgl[0][p];
        B[4][i] = shapeDgl[2][p];
        break;

      case 2:  /** Columns 3, 6, 9, ..., 24 */
        B[2][i] = shapeDgl[2][p];
        B[4][i] = shapeDgl[1][p];
        B[5][i] = shapeDgl[0][p];
        break;    
    }
  }

}




template<class TBaseClass>
void
Element3DLinearElasticity<TBaseClass>
::GetMassMatrix(MatrixType& Me) const
{
  // Call the parent's get matrix function
  Superclass::GetMassMatrix(Me);

  // Since parent class doesn't have the material properties,
  // we need to adjust Me matrix here for the density of the element.
  Me=Me*m_mat->RhoC;
}




template<class TBaseClass>
void
Element3DLinearElasticity<TBaseClass>
::GetMaterialMatrix(MatrixType& D) const
{
  D.set_size(6,6);
  D.fill(0.0);
 
  Float lambda =  m_mat->nu;
  Float mu =  m_mat->E;
    

  // Set the elements in the top left quadrant 
  for (int j=0; j < 3; j++) {
    for (int k=0; k < 3; k++) {
      D[j][k] = lambda;
    }
  }

  // Set the diagonal elements 
  for (int k=0; k < 3; k++) {
    D[k][k] = lambda + 2.0*mu;
  }
  for (int k=3; k < 6; k++) {
    D[k][k] = mu;
  }

/*
  D.fill(0.);

  D[0][0] = 1.;
  D[0][1] = (-1.0)*mu;
  D[0][2] = (-1.0)*mu;

  D[1][0] = (-1.0)*mu;
  D[1][1] = 1.;
  D[1][2] = (-1.0)*mu;

  D[2][0] = (-1.0)*mu;
  D[2][1] = (-1.0)*mu;
  D[2][2] = 1.;
  
  D[3][3] = 1.+mu;
  D[4][4] = 1.+mu;
  D[5][5] = 1.+mu;

  D=D*1.0/lambda;

*/
}




template<class TBaseClass>
void
Element3DLinearElasticity<TBaseClass>
::Read( std::istream& f, void* info )
{
  int n;
  /*
   * Convert the info pointer to a usable objects
   */
  ReadInfoType::MaterialArrayPointer mats=static_cast<ReadInfoType*>(info)->m_mat;


  /* first call the parent's read function */
  Superclass::Read(f,info);

  try
  {
    /*
     * Read and set the material pointer
     */
    SkipWhiteSpace(f); f>>n; if(!f) goto out;
    m_mat=dynamic_cast<const MaterialLinearElasticity*>( &*mats->Find(n));

  }
  catch ( FEMExceptionObjectNotFound e )
  {
    throw FEMExceptionObjectNotFound(__FILE__,__LINE__,"Element3DLinearElasticity::Read()",e.m_baseClassName,e.m_GN);
  }

  // Check if the material object was of correct class
  if(!m_mat)
  {
    throw FEMExceptionWrongClass(__FILE__,__LINE__,"Element3DStress::Read()");
  }


out:

  if( !f )
  { 
    throw FEMExceptionIO(__FILE__,__LINE__,"Element3DLinearElasticity::Read()","Error reading FEM element!");
  }

}



/*
 * Write the element to the output stream.
 */
template<class TBaseClass>
void
Element3DLinearElasticity<TBaseClass>
::Write( std::ostream& f ) const {

  // First call the parent's write function
  Superclass::Write(f);

  /*
   * then write the actual data (material number)
   * We also add some comments in the output file
   */
  f<<"\t"<<m_mat->GN<<"\t% MaterialLinearElasticity ID\n";

  // check for errors
  if (!f)
  { 
    throw FEMExceptionIO(__FILE__,__LINE__,"Element3DLinearElasticity::Write()","Error writing FEM element!");
  }

}




#ifdef _MSC_VER
// Declare a static dummy function to prevent a MSVC 6.0 SP5 from crashing.
// I have no idea why things don't work when this is not declared, but it
// looks like this declaration makes compiler forget about some of the
// troubles it has with templates.
static void Dummy( void );
#endif // #ifdef _MSC_VER

}} // end namespace itk::fem

#endif // #ifndef __itkFEMElement3DLinearElasticity_txx

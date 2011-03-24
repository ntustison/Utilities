/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkNormalizedCutsSegmentationImageFilter.txx,v $
  Language:  C++
  Date:      
  Version:   $Revision: 1.1.1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkNormalizedCutsSegmentationImageFilter_txx_
#define __itkNormalizedCutsSegmentationImageFilter_txx_

#include "itkNormalizedCutsSegmentationImageFilter.h"

#include "itkConstantScalarOperator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkNeighborhoodInnerProduct.h"

#include "vnl/vnl_math.h"

namespace itk
{

template <class TInputImage, class TLabelImage>
NormalizedCutsSegmentationImageFilter<TInputImage, TLabelImage>
::NormalizedCutsSegmentationImageFilter()
{
  this->m_Radius.Fill( 10 );
  this->m_SpatialSigma = 10.0;
  this->m_DataSigma = 50.0;

  this->m_NumberOfClasses = 2;
  this->m_MaskImage = NULL;
  this->m_ForegroundValue = 1;
  this->m_NumberOfSplittingPoints = 20;

  this->m_PartitionStrategy = Zero;
}

/** Generate the data */
template <class TInputImage, class TLabelImage>
void
NormalizedCutsSegmentationImageFilter<TInputImage, TLabelImage>
::GenerateData()
{
  this->GenerateEigenSystemFromInputImage();
  this->SolveEigensystem();
//  this->MulticlassSpectralClustering();
  this->EigenVectorClustering();
}

template <class TInputImage, class TLabelImage>
void
NormalizedCutsSegmentationImageFilter<TInputImage, TLabelImage>
::GenerateEigenSystemFromInputImage()
{
  unsigned long numberOfVoxels = 1;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    numberOfVoxels *= this->GetInput()->GetLargestPossibleRegion().GetSize()[i];
    }

  /** Fill in the sparse matrices **/  
  this->m_W.resize( numberOfVoxels, numberOfVoxels );
  this->m_D.resize( numberOfVoxels, numberOfVoxels );
  this->m_E.resize( numberOfVoxels, numberOfVoxels );

  SizeType size = this->GetInput()->GetLargestPossibleRegion().GetSize();

  unsigned long numberOfNeighbors = 1;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    numberOfNeighbors *= ( 2*this->m_Radius[i] + 1 );
    }

  typedef typename NeighborhoodAlgorithm
    ::ImageBoundaryFacesCalculator<InputImageType> FaceCalculatorType;
  
  FaceCalculatorType faceCalculator;
  typename FaceCalculatorType::FaceListType faceList;
  typename FaceCalculatorType::FaceListType::iterator fit;

  faceList = faceCalculator( this->GetInput(), this->GetInput()->GetLargestPossibleRegion(),
                             this->m_Radius );
 
  for ( fit = faceList.begin(); fit != faceList.end(); ++fit )
    {
    ConstNeighborhoodIteratorType It( this->m_Radius, this->GetInput(), *fit );
    for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      if ( this->m_MaskImage && 
           ( this->m_MaskImage->GetPixel( It.GetIndex() ) != this->m_ForegroundValue ) )
        {
        continue;
        }
      unsigned long iLabel = this->IndexToNumber( It.GetIndex(), size ); 
      RealType d = 0.0;
      for ( unsigned int j = 0; j < numberOfNeighbors; j++ )
        {
        unsigned long jLabel = this->IndexToNumber( It.GetIndex( j ), size );
        if ( iLabel == jLabel || 
             !this->GetInput()->GetLargestPossibleRegion().IsInside( It.GetIndex( j ) ) 
             || ( this->m_MaskImage && 
             ( this->m_MaskImage->GetPixel( It.GetIndex( j ) ) != this->m_ForegroundValue ) ) )
          {
          continue;
          }
        RealType weight = this->EvaluateWeightFunction( It.GetIndex(), It.GetIndex( j ) );
        this->m_W( iLabel, jLabel ) = weight;  
        d += weight;
        }
      
      this->m_D( iLabel, iLabel ) = d;
      this->m_E( iLabel, iLabel ) = 1.0 / vcl_sqrt( d + vnl_math::eps );
      }
    }  
}

template <class TInputImage, class TLabelImage>
typename NormalizedCutsSegmentationImageFilter<TInputImage, TLabelImage>:: RealType
NormalizedCutsSegmentationImageFilter<TInputImage, TLabelImage>
::EvaluateWeightFunction( IndexType index1, IndexType index2 )
{
  RealType iValue = static_cast<RealType>( this->GetInput()->GetPixel( index1 ) );
  InputPointType iPoint;
  this->GetInput()->TransformIndexToPhysicalPoint( index1, iPoint ); 

  RealType jValue = static_cast<RealType>( this->GetInput()->GetPixel( index2 ) );
  InputPointType jPoint;
  this->GetInput()->TransformIndexToPhysicalPoint( index2, jPoint ); 

  RealType weight = vcl_exp( -0.5*vnl_math_sqr( ( iValue - jValue ) / this->m_DataSigma ) ) 
    * vcl_exp( -0.5*( iPoint - jPoint ).GetSquaredNorm() / vnl_math_sqr( this->m_SpatialSigma ) );
  
  return weight;
}

template <class TInputImage, class TLabelImage>
void
NormalizedCutsSegmentationImageFilter<TInputImage, TLabelImage>
::SolveEigensystem()
{
  SparseMatrixType A;
  SparseMatrixType B; 
  SparseMatrixType C;

  this->m_D.subtract( this->m_W, A );
  this->m_E.mult( A, B );
  B.mult( this->m_E, C );

  this->m_EigenSystem = EigenSystemType::New();
  this->m_EigenSystem->SetNumberOfEigenPairs( this->m_NumberOfClasses );
  this->m_EigenSystem->SetNumberOfLanczosVectors( 35 );
  this->m_EigenSystem->SetSolveForSmallestEigenValues( true );
  this->m_EigenSystem->SetTolerance( 1e-6 );
  this->m_EigenSystem->SetMatrix( &C );
  this->m_EigenSystem->Update();
}

template <class TInputImage, class TLabelImage>
void
NormalizedCutsSegmentationImageFilter<TInputImage, TLabelImage>
::EigenVectorClustering()
{
// Use second smallest eigenvector to partition the graph
  
  typename LabelImageType::Pointer output = LabelImageType::New();
  output->SetRegions( this->GetInput()->GetLargestPossibleRegion() );
  output->SetOrigin( this->GetInput()->GetOrigin() );
  output->SetSpacing( this->GetInput()->GetSpacing() );
  output->Allocate();

  typename EigenSystemType::VectorType minCutVector = this->m_EigenSystem->GetEigenVector( 1 );

  switch( this->m_PartitionStrategy )
    {
    case Zero: default: 
      for ( unsigned int i = 0; i < minCutVector.size(); i++ )
        {
        if ( minCutVector[i] <= 0 )
          {
          minCutVector[i] = -1;
          }
        else
          {
          minCutVector[i] = 1;
          }      
        }
      break;

    case Mean:
      RealType mean = minCutVector.mean();     
      
      for ( unsigned int i = 0; i < minCutVector.size(); i++ )
        {
        if ( minCutVector[i] <= mean )
          {
          minCutVector[i] = -1;
          }
        else
          {
          minCutVector[i] = 1;
          }      
        }
      break;

    case SplittingPoints: 
      SparseMatrixType A;
      this->m_D.subtract( this->m_W, A );
   
      RealType minValue = this->m_EigenSystem->GetEigenVector( 1 ).min_value();
      RealType maxValue = this->m_EigenSystem->GetEigenVector( 1 ).max_value();
      RealType dx = vnl_math_abs( maxValue - minValue ) 
                    / static_cast<RealType>( this->m_NumberOfSplittingPoints );
    
      RealType minCutValue = NumericTraits<RealType>::max();
      vnl_vector<RealType> minCutVector(
        this->m_EigenSystem->GetEigenVector( 1 ).size() );
      minCutVector.fill( 0.0 );

      if ( minValue == maxValue )
        {
        break;
        }  

      for ( RealType x = minValue+dx; x <= maxValue-dx; x+=dx )
        {
        vnl_vector<RealType> binaryEigenVector 
          = this->m_EigenSystem->GetEigenVector( 1 );
    
        RealType bNumerator = 0;
        RealType bDenominator = 0;
        for ( unsigned int i = 0; i < binaryEigenVector.size(); i++ )
          {
          if ( binaryEigenVector[i] < x )
            {
            bNumerator += this->m_D( i, i );
            }
          else
            {
            bDenominator += this->m_D( i, i );
            }      
          }
        RealType b = bNumerator / bDenominator;
    
        for ( unsigned int i = 0; i < binaryEigenVector.size(); i++ )
          {
          if ( binaryEigenVector[i] < x )
            {
            binaryEigenVector[i] = 1;
            }
          else
            {
            binaryEigenVector[i] = -b;
            }      
          }
    
        RealType denominator = binaryEigenVector.squared_magnitude();
    
        if ( denominator == 0 )
          {
          continue;
          }
     
        vnl_vector<RealType> result;
        A.mult( binaryEigenVector, result );
    
        RealType numerator = 0.0;
        for ( unsigned int i = 0; i < binaryEigenVector.size(); i++ )
          {
          numerator += ( result[i] * binaryEigenVector[i] );
          }   
     
        if ( minCutValue >= ( numerator / denominator ) )
          {
          minCutValue = ( numerator / denominator );
          minCutVector = binaryEigenVector;
          }
        }
      break;
    }  

  SizeType size = this->GetInput()->GetLargestPossibleRegion().GetSize();

  for ( unsigned long i = 0; i < minCutVector.size(); i++ )
    {
    if ( minCutVector[i] < 0 )
      { 
      output->SetPixel( this->NumberToIndex( i, size ), 0 );    
      } 
    else
      {
      output->SetPixel( this->NumberToIndex( i, size ), 1 );    
      } 
    }

  this->SetNthOutput( 0, output );
}

template <class TInputImage, class TLabelImage>
typename NormalizedCutsSegmentationImageFilter<TInputImage, TLabelImage>::RealImageType::Pointer
NormalizedCutsSegmentationImageFilter<TInputImage, TLabelImage>
::GetEigenImage( unsigned int n )
{
  typename RealImageType::Pointer output = RealImageType::New();
  output->SetRegions( this->GetInput()->GetLargestPossibleRegion() );
  output->SetOrigin( this->GetInput()->GetOrigin() );
  output->SetSpacing( this->GetInput()->GetSpacing() );
  output->Allocate();
  output->FillBuffer( 0 );

  if ( n < this->m_NumberOfClasses )
    {
     vnl_vector<RealType> eigenVector 
       = this->m_EigenSystem->GetEigenVector( n );

    typename RealImageType::Pointer output = RealImageType::New();
    output->SetRegions( this->GetInput()->GetLargestPossibleRegion() );
    output->SetOrigin( this->GetInput()->GetOrigin() );
    output->SetSpacing( this->GetInput()->GetSpacing() );
    output->Allocate();
  
    SizeType size = this->GetInput()->GetLargestPossibleRegion().GetSize();

    for ( unsigned long i = 0; i < eigenVector.size(); i++ )
      {
      output->SetPixel( this->NumberToIndex( i, size ), eigenVector[i] );    
      }
    }

  return output;  
}

template <class TInputImage, class TLabelImage>
void
NormalizedCutsSegmentationImageFilter<TInputImage, TLabelImage>
::MulticlassSpectralClustering()
{
/*
  SparseMatrixType A;

  this->m_D.subtract( this->m_W, A );

// Use second smallest eigenvector to partition the graph
  
  typename LabelImageType::Pointer output = LabelImageType::New();
  output->SetRegions( this->GetInput()->GetLargestPossibleRegion() );
  output->SetOrigin( this->GetInput()->GetOrigin() );
  output->SetSpacing( this->GetInput()->GetSpacing() );
  output->Allocate();

  SizeType size = this->GetInput()->GetLargestPossibleRegion().GetSize();

  RealType minValue = this->m_EigenSystem.GetEigenVector( 1 ).min_value();
  RealType maxValue = this->m_EigenSystem.GetEigenVector( 1 ).max_value();
  RealType dx = vnl_math_abs( maxValue - minValue ) 
                / static_cast<RealType>( this->m_NumberOfSplittingPoints );

  RealType minCutValue = NumericTraits<RealType>::max();
  vnl_vector<RealType> minCutVector(
    this->m_EigenSystem.GetEigenVector( 1 ).size() );
  minCutVector.fill( 0.0 );

  for ( RealType x = minValue+dx; x <= maxValue-dx; x+=dx )
    {
    vnl_vector<RealType> binaryEigenVector 
      = this->m_EigenSystem.GetEigenVector( 1 );

    std::cout << x << " out of " << maxValue-dx << " "<< std::endl;


    RealType bNumerator = 0;
    RealType bDenominator = 0;
    for ( unsigned int i = 0; i < binaryEigenVector.size(); i++ )
      {
      if ( binaryEigenVector[i] < x )
        {
        bNumerator += this->m_D( i, i );
        }
      else
        {
        bDenominator += this->m_D( i, i );
        }      
      }
    RealType b = bNumerator / bDenominator;

    for ( unsigned int i = 0; i < binaryEigenVector.size(); i++ )
      {
      if ( binaryEigenVector[i] < x )
        {
        binaryEigenVector[i] = 1;
        }
      else
        {
        binaryEigenVector[i] = -b;
        }      
      }

    RealType denominator = binaryEigenVector.squared_magnitude();

    if ( denominator == 0 )
      {
      continue;
      }
 
    vnl_vector<RealType> result;
    A.mult( binaryEigenVector, result );

    RealType numerator = 0.0;
    for ( unsigned int i = 0; i < binaryEigenVector.size(); i++ )
      {
      numerator += ( result[i] * binaryEigenVector[i] );
      }   
 
    std::cout << x << " out of " << maxValue-dx << ": " << numerator / denominator << std::endl;


    if ( minCutValue >= ( numerator / denominator ) )
      {
      minCutValue = ( numerator / denominator );
      minCutVector = binaryEigenVector;
      }
    }   

  for ( unsigned long i = 0; i < minCutVector.size(); i++ )
    {
    if ( minCutVector[i] < 0 )
      { 
      output->SetPixel( this->NumberToIndex( i, size ), 0 );    
      } 
    else
      {
      output->SetPixel( this->NumberToIndex( i, size ), 1 );    
      } 
    }

  this->SetNthOutput( 0, output );
*/
}

template <class TInputGraph, class TLabelImage>
void
NormalizedCutsSegmentationImageFilter<TInputGraph, TLabelImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
}

} // end namespace itk

#endif
    
    

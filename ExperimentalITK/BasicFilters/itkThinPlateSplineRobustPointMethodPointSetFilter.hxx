#ifndef __itkThinPlateSplineRobustPointMethodPointSetFilter_hxx
#define __itkThinPlateSplineRobustPointMethodPointSetFilter_hxx

#include "itkThinPlateSplineRobustPointMethodPointSetFilter.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkThinPlateSplineKernelTransform.h"
#include "itkThinPlateR2LogRSplineKernelTransform.h"

#include "vnl/vnl_math.h"

#include "fstream.h"

namespace itk
{

template <class TPointSet, class TOutputImage>
ThinPlateSplineRobustPointMethodPointSetFilter<TPointSet, TOutputImage>
::ThinPlateSplineRobustPointMethodPointSetFilter()
{
  this->SetNumberOfRequiredInputs( 2 );

  this->m_FinalTemperature = NumericTraits<RealType>::max();
  this->m_InitialTemperature = NumericTraits<RealType>::max();
  this->m_InitialLambda = 1.0;
  this->m_AnnealingRate = 0.93;
  this->m_NumberOfIterationsPerTemperature = 5;

  this->m_Origin.Fill( 0 );
  this->m_Size.Fill( 64 );
  this->m_Spacing.Fill( 1 );

  this->m_UseBoundingBox = true;

  this->m_SolveSimplerLeastSquaresProblem = true;
}

template <class TPointSet, class TOutputImage>
ThinPlateSplineRobustPointMethodPointSetFilter<TPointSet, TOutputImage>
::~ThinPlateSplineRobustPointMethodPointSetFilter()
{  
}

template <class TPointSet, class TOutputImage>
void
ThinPlateSplineRobustPointMethodPointSetFilter<TPointSet, TOutputImage>
::GenerateData()
{
  this->Initialize();

  this->m_CurrentTemperature = this->m_InitialTemperature;

  this->m_CurrentLambda = this->m_InitialLambda * this->m_CurrentTemperature * 
    static_cast<RealType>( this->GetInput( 0 )->GetNumberOfPoints() ) / this->m_AnnealingRate;

  while ( this->m_CurrentTemperature >= this->m_FinalTemperature )
    {

    for ( unsigned int i = 0; i < this->m_NumberOfIterationsPerTemperature; i++ )
      {
      /**
       * Step 1:  Update the correspondence matrix
       */
      this->UpdateCorrespondenceMatrix();
  
      /**
       * Step 2:  Update the transformation
       */
      this->UpdateTransformation();  
      }

    this->m_CurrentTemperature *= this->m_AnnealingRate;
    this->m_CurrentLambda *= this->m_AnnealingRate;

    std::cout << this->m_CurrentTemperature << ", " << this->m_CurrentLambda << std::endl;
    if ( this->m_CurrentTemperature < 0.0942 )
      {  
      this->VisualizeCurrentState();
      exit( 0 );
      }

    }    



  /**
   * Generate output
   */   

  typename PointSetType::Pointer targetLandmarks = PointSetType::New();
  targetLandmarks->Initialize();
  typename PointSetType::Pointer sourceLandmarks = PointSetType::New();
  sourceLandmarks->Initialize();
  
  unsigned int count = 0;
  for ( unsigned int i = 0; i < this->m_VPoints->GetNumberOfPoints(); i++ )
    {
    typename InputPointSetType::PointType V;
    this->GetInput( 1 )->GetPoint( i, &V );
  
    if ( this->m_SolveSimplerLeastSquaresProblem )
      { 
      typename PointSetType::PointType Y;
      Y.Fill( 0 );
      for ( unsigned int j = 0; j < this->GetInput( 0 )->GetNumberOfPoints(); j++ ) 
        {
        RealType m = this->m_CorrespondenceMatrix( i, j );
        if ( m > 0 )
          {  
          typename InputPointSetType::PointType X;
          this->GetInput( 0 )->GetPoint( j, &X );
          for ( unsigned int d = 0; d < Dimension; d++ )
            {
            Y[d] += ( X[d] * m ); 
            }
          }     
        }
      targetLandmarks->SetPoint( count, Y );  
      sourceLandmarks->SetPoint( count, V ); 
      count++;
      }
    else
      {
      std::cerr << "Not implemented." << std::endl;
      exit( 0 );
      }     
    } 

  typename TransformType::Pointer tps;
  if ( Dimension == 2 )
    {
    tps = ThinPlateR2LogRSplineKernelTransform<RealType, Dimension>::New();
    }
  else if ( Dimension == 3 )
    {
    tps = ThinPlateSplineKernelTransform<RealType, Dimension>::New();
    }
  tps->SetSourceLandmarks( sourceLandmarks );
  tps->SetTargetLandmarks( targetLandmarks );
  tps->SetStiffness( this->m_CurrentLambda / this->m_AnnealingRate );
  tps->ComputeWMatrix();
 
  typename OutputImageType::Pointer output = OutputImageType::New();
  output->SetRegions( this->m_Size );
  output->SetSpacing( this->m_Spacing );
  output->SetOrigin( this->m_Origin );
  output->Allocate();

  ImageRegionIteratorWithIndex<OutputImageType> It( output,
    output->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    typename OutputImageType::PointType point;
    output->TransformIndexToPhysicalPoint( It.GetIndex(), point );
    
    typename PointSetType::PointType X;
    for ( unsigned int d = 0; d < Dimension; d++ )
      {
      X[d] = point[d];
      } 
    typename PointSetType::PointType Y = tps->TransformPoint( X );
    
    typename OutputImageType::PixelType V;
    for ( unsigned int d = 0; d < Dimension; d++ )
      {
      V[d] = Y[d] - X[d]; 
      }
    It.Set( V ); 
    }  

  this->GraftOutput( output );
}      
         
template <class TPointSet, class TOutputImage>
void
ThinPlateSplineRobustPointMethodPointSetFilter<TPointSet, TOutputImage>
::Initialize()
{
  /**
   * Create a copy of the moving point set since it will be updated.
   */

  this->m_VPoints = InputPointSetType::New();
  this->m_VPoints->Initialize();

  for ( unsigned int i = 0; i < this->GetInput( 1 )->GetNumberOfPoints(); i++ )
    {
    typename InputPointSetType::PointType point;      
    typename InputPointSetType::PixelType pixel;      
    this->GetInput( 1 )->GetPoint( i, &point );   
    this->GetInput( 1 )->GetPointData( i, &pixel );   
    this->m_VPoints->SetPoint( i, point ); 
    this->m_VPoints->SetPointData( i, pixel );
    } 

  /**
   * If UseBoundingBox, calculate the size, origin of the region of transformation
   */
  if ( this->m_UseBoundingBox )
    {
    for ( unsigned int i = 0; i < Dimension; i++ )
      {
      RealType minimum = 
        vnl_math_min( this->GetInput( 0 )->GetBoundingBox()->GetMinimum()[i], 
                      this->GetInput( 1 )->GetBoundingBox()->GetMinimum()[i] );
      RealType maximum = 
        vnl_math_max( this->GetInput( 0 )->GetBoundingBox()->GetMaximum()[i], 
                      this->GetInput( 1 )->GetBoundingBox()->GetMaximum()[i] );

      /** 
       * Expand the bounding box to ensure coverage  
       */
      RealType expandedMinimum = minimum - 0.5 * ( maximum - minimum );
      RealType expandedMaximum = maximum + 0.5 * ( maximum - minimum );
      this->m_Origin[i] = expandedMinimum;
      this->m_Size[i] = static_cast<unsigned int>( ( expandedMaximum - expandedMinimum ) / this->m_Spacing[i] ) + 1;       
      }  
    }   

  /**
   * Initialize correspondence matrix and outlier row/column
   */
  RealType K = static_cast<RealType>( this->GetInput( 1 )->GetNumberOfPoints() );
  RealType N = static_cast<RealType>( this->GetInput( 0 )->GetNumberOfPoints() );

  this->m_CorrespondenceMatrix.SetSize( K, N );
  this->m_OutlierColumn.SetSize( K );
  this->m_OutlierRow.SetSize( N );

  this->m_OutlierColumn.Fill( 0.1 / K );
  this->m_OutlierRow.Fill( 0.1 / K );

  if ( this->m_FinalTemperature == NumericTraits<RealType>::max() || 
       this->m_InitialTemperature == NumericTraits<RealType>::max() )
    {
    this->CalculateInitialAndFinalTemperatures();
    }

}

template <class TPointSet, class TOutputImage>
void
ThinPlateSplineRobustPointMethodPointSetFilter<TPointSet, TOutputImage>
::CalculateInitialAndFinalTemperatures()
{
  RealType maxSquaredDistance = 0;
  RealType sumMinSquaredDistance = 0;

  for ( unsigned int i = 0; i < this->GetInput( 1 )->GetNumberOfPoints(); i++ )
    {
    typename InputPointSetType::PointType V;
    this->GetInput( 1 )->GetPoint( i, &V );    

    RealType minSquaredDistance = NumericTraits<RealType>::max();  
  
    for ( unsigned int j = 0; j < this->GetInput( 0 )->GetNumberOfPoints(); j++ )
      {
      typename InputPointSetType::PointType X;      
      this->GetInput( 0 )->GetPoint( j, &X );      
      
      RealType squaredDistance = ( X - V ).GetSquaredNorm();

      if ( squaredDistance > maxSquaredDistance )
        {
        maxSquaredDistance = squaredDistance;
        } 
      if ( squaredDistance < minSquaredDistance )
        {
        minSquaredDistance = squaredDistance;
        } 
      }  
    if ( minSquaredDistance != NumericTraits<RealType>::max() )
      {
      sumMinSquaredDistance += minSquaredDistance;
      }
    } 

  this->m_InitialTemperature = maxSquaredDistance;
  this->m_FinalTemperature = sumMinSquaredDistance 
    / static_cast<RealType>( this->GetInput( 1 )->GetNumberOfPoints() );
}

template <class TPointSet, class TOutputImage>
void
ThinPlateSplineRobustPointMethodPointSetFilter<TPointSet, TOutputImage>
::UpdateCorrespondenceMatrix()
{
  for ( unsigned int i = 0; i < this->m_VPoints->GetNumberOfPoints(); i++ )
    {
    typename InputPointSetType::PointType V;
    this->m_VPoints->GetPoint( i, &V );    

    for ( unsigned int j = 0; j < this->GetInput( 0 )->GetNumberOfPoints(); j++ )
      {
      typename InputPointSetType::PointType X;      
      this->GetInput( 0 )->GetPoint( j, &X );      
 
      this->m_CorrespondenceMatrix( i, j )  
        = vcl_exp( -( X - V ).GetSquaredNorm() / this->m_CurrentTemperature );
      }  
    } 

  RealType K = static_cast<RealType>( this->GetInput( 1 )->GetNumberOfPoints() );
  if ( this->m_CurrentTemperature == this->m_InitialTemperature )
    {
    this->m_OutlierColumn.Fill( 0.0017 /*0.1 / K*/ );
    this->m_OutlierRow.Fill( 0.0017 /*0.1 / K*/ );
    }
  else
    {  
    this->m_OutlierColumn.Fill( 0.3078 /*0.1 / K*/ );
    this->m_OutlierRow.Fill( 0.3078 /*0.1 / K*/ );
    }

  /**
   * Normalize correspondence matrix
   */
 
  RealType epsilon = 0.05;
  RealType deviation = NumericTraits<RealType>::max(); 

  OutlierVectorType rowSum;
  OutlierVectorType columnSum;

  unsigned int maximumNumberOfIterations = 10;
  unsigned int iterations = 0;

  while ( deviation > epsilon*epsilon && iterations++ < maximumNumberOfIterations )
    {
    /**
     * do column normalization
     */
    columnSum = this->m_OutlierRow;

    for ( unsigned int j = 0; j < this->m_CorrespondenceMatrix.Cols(); j++ )
      {
      for ( unsigned int i = 0; i < this->m_CorrespondenceMatrix.Rows(); i++ )
        {
        columnSum[j] += this->m_CorrespondenceMatrix( i, j );
        }
      } 
    for ( unsigned int j = 0; j < this->m_CorrespondenceMatrix.Cols(); j++ )
      {
      for ( unsigned int i = 0; i < this->m_CorrespondenceMatrix.Rows(); i++ )
        {
        this->m_CorrespondenceMatrix( i, j ) /= columnSum[j];
        }
      this->m_OutlierRow[j] /= columnSum[j];
      } 

    /**
     * do row normalization
     */
    rowSum = this->m_OutlierColumn;

    for ( unsigned int i = 0; i < this->m_CorrespondenceMatrix.Rows(); i++ )
      {
      for ( unsigned int j = 0; j < this->m_CorrespondenceMatrix.Cols(); j++ )
        {
        rowSum[i] += this->m_CorrespondenceMatrix( i, j );
        }
      } 
    for ( unsigned int i = 0; i < this->m_CorrespondenceMatrix.Rows(); i++ )
      {
      for ( unsigned int j = 0; j < this->m_CorrespondenceMatrix.Cols(); j++ )
        {
        this->m_CorrespondenceMatrix( i, j ) /= rowSum[i];
        }
      this->m_OutlierColumn[i] /= rowSum[i];
      } 
    

    /**
     * Calculate current deviation from 1
     */
    rowSum -= 1.0;
    columnSum -= 1.0;

    deviation = ( rowSum.GetSquaredNorm() + columnSum.GetSquaredNorm() )
      / ( static_cast<RealType>( rowSum.GetSize() + columnSum.GetSize() ) );
    } 

}

template <class TPointSet, class TOutputImage>
void
ThinPlateSplineRobustPointMethodPointSetFilter<TPointSet, TOutputImage>
::UpdateTransformation()
{
  typename PointSetType::Pointer targetLandmarks = PointSetType::New();
  targetLandmarks->Initialize();
  typename PointSetType::Pointer sourceLandmarks = PointSetType::New();
  sourceLandmarks->Initialize();
  
  unsigned int count = 0;
  for ( unsigned int i = 0; i < this->GetInput( 1 )->GetNumberOfPoints(); i++ )
    {
    typename InputPointSetType::PointType V;
    this->GetInput( 1 )->GetPoint( i, &V );
  
    if ( this->m_SolveSimplerLeastSquaresProblem )
      { 
      typename PointSetType::PointType Y;
      Y.Fill( 0 );
      for ( unsigned int j = 0; j < this->GetInput( 0 )->GetNumberOfPoints(); j++ ) 
        {
        RealType m = this->m_CorrespondenceMatrix( i, j );
        if ( m > 0 )
          {  
          typename InputPointSetType::PointType X;
          this->GetInput( 0 )->GetPoint( j, &X );
          for ( unsigned int d = 0; d < Dimension; d++ )
            {
            Y[d] += ( X[d] * m ); 
            }
          }     
        }
      targetLandmarks->SetPoint( count, Y );  
      sourceLandmarks->SetPoint( count, V ); 
      count++;
      }
    else
      {
      std::cerr << "Not implemented." << std::endl;
      exit( 0 );
      }     
    } 

  typename TransformType::Pointer tps;
  if ( Dimension == 2 )
    {
    tps = ThinPlateR2LogRSplineKernelTransform<RealType, Dimension>::New();
    }
  else if ( Dimension == 3 )
    {
    tps = ThinPlateSplineKernelTransform<RealType, Dimension>::New();
    }
  tps->SetSourceLandmarks( sourceLandmarks );
  tps->SetTargetLandmarks( targetLandmarks );
  tps->SetStiffness( this->m_CurrentLambda );
  tps->ComputeWMatrix();
 
  /**
   * Update the V points
   */
  for ( unsigned int i = 0; i < this->m_VPoints->GetNumberOfPoints(); i++ )
    {
    typename InputPointSetType::PointType V;
    this->GetInput( 1 )->GetPoint( i, &V );

    typename PointSetType::PointType fV;
    for ( unsigned d = 0; d < Dimension; d++ )
      {
      fV[d] = V[d];
      }
    
    fV = tps->TransformPoint( fV );
    for ( unsigned d = 0; d < Dimension; d++ )
      {
      V[d] = fV[d];
      }
    this->m_VPoints->SetPoint( i, V );
    }      
}

template <class TPointSet, class TOutputImage>
void
ThinPlateSplineRobustPointMethodPointSetFilter<TPointSet, TOutputImage>
::VisualizeCurrentState()
{
  std::ofstream str( "CurrentPoints.txt" );
  str << "0 0 0 0" << std::endl;

  for ( unsigned int i = 0; i < this->m_VPoints->GetNumberOfPoints(); i++ )
    {
    typename InputPointSetType::PointType V;
    this->m_VPoints->GetPoint( i, &V );
    str << V[0] << " " << V[1] << " 0 " << i+1 << std::endl;
    } 

  str << "0 0 0 0" << std::endl;
  str.close();

  RealType K = static_cast<RealType>( this->GetInput( 1 )->GetNumberOfPoints() );
  RealType N = static_cast<RealType>( this->GetInput( 0 )->GetNumberOfPoints() );

  std::ofstream str2( "CurrentInfluence.txt" );
  str2 << "0 0 0 0" << std::endl;

  for ( unsigned int i = 0; i < this->m_VPoints->GetNumberOfPoints(); i++ )
    {
    typename InputPointSetType::PointType V;
    this->m_VPoints->GetPoint( i, &V );

    for ( RealType theta = 0.0; theta <= 2*vnl_math::pi; theta+= 0.01 )
      {
      str2 << V[0] + vcl_sqrt( this->m_CurrentTemperature ) * vcl_cos( theta )<< " "
           << V[1] + vcl_sqrt( this->m_CurrentTemperature ) * vcl_sin( theta )<< " "
           << "0 " << i+1 << std::endl;
      } 
    } 
  str2 << "0 0 0 0" << std::endl;
  str2.close();


  std::ofstream str3( "CurrentNeighbors.txt" );
  str3 << "0 0 0 0" << std::endl;
  unsigned int count = 0;

  for ( unsigned int i = 0; i < this->m_VPoints->GetNumberOfPoints(); i++ )
    {
    typename InputPointSetType::PointType V;
    this->m_VPoints->GetPoint( i, &V );
    for ( unsigned int j = 0; j < this->GetInput( 0 )->GetNumberOfPoints(); j++ )
      {
      if ( this->m_CorrespondenceMatrix( i, j ) > 1.0 / static_cast<RealType>( K ) )
        {
        typename InputPointSetType::PointType X;
        this->GetInput( 0 )->GetPoint( j, &X );
        str3 << X[0] << " " << X[1] << " 0 " << i+1 << std::endl;      
        str3 << V[0] << " " << V[1] << " 0 " << i+1 << std::endl;   
        } 
      } 
    } 
  str3 << "0 0 0 0" << std::endl;
  str3.close();


  std::ofstream str4( "CorrespondenceMatrix.txt" );

  for ( unsigned int i = 0; i < this->m_CorrespondenceMatrix.Rows(); i++ )
    {
    for ( unsigned int j = 0; j < this->m_CorrespondenceMatrix.Cols(); j++ )
      {
      str4 << this->m_CorrespondenceMatrix( i, j ) << " ";
      } 
    str4 << std::endl;   
    } 
  str4.close();

}



/**
 * Standard "PrintSelf" method
 */
template <class TPointSet, class TOutputImage>
void
ThinPlateSplineRobustPointMethodPointSetFilter<TPointSet, TOutputImage>
::PrintSelf(
  std::ostream& os, 
  Indent indent) const
{
}



}  //end namespace itk

#endif

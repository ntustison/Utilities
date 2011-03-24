#ifndef __itkBSplineRobustPointMethodPointSetFilter_txx
#define __itkBSplineRobustPointMethodPointSetFilter_txx

#include "itkBSplineRobustPointMethodPointSetFilter.h"

#include "itkBSplineControlPointImageFilter.h"

#include "vnl/vnl_math.h"

#include "fstream.h"

namespace itk
{

template <class TPointSet, class TOutputImage>
BSplineRobustPointMethodPointSetFilter<TPointSet, TOutputImage>
::BSplineRobustPointMethodPointSetFilter()
{
  this->SetNumberOfRequiredInputs( 2 );

  this->m_SplineOrder = 3;
  this->m_NumberOfLevels = 1;
  this->m_NumberOfControlPoints.Fill( this->m_SplineOrder+1 );

  this->m_ControlPointLattice = ControlPointLatticeType::New();
  this->m_ControlPointLattice = NULL;

  this->m_FinalTemperature = NumericTraits<RealType>::max();
  this->m_InitialTemperature = NumericTraits<RealType>::max();
  this->m_AnnealingRate = 0.9;
  this->m_NumberOfIterationsPerTemperature = 1;

  this->m_Origin.Fill( 0 );
  this->m_Size.Fill( 64 );
  this->m_Spacing.Fill( 1 );

  this->m_UseBoundingBox = true;

  this->m_SolveSimplerLeastSquaresProblem = false;
}

template <class TPointSet, class TOutputImage>
BSplineRobustPointMethodPointSetFilter<TPointSet, TOutputImage>
::~BSplineRobustPointMethodPointSetFilter()
{  
}

template <class TPointSet, class TOutputImage>
void
BSplineRobustPointMethodPointSetFilter<TPointSet, TOutputImage>
::GenerateData()
{
  this->m_CurrentTemperature = NumericTraits<RealType>::max();

  this->Initialize();

  this->m_CurrentTemperature = 0.01 * this->m_InitialTemperature;

  this->VisualizeCurrentState();

  while ( this->m_CurrentTemperature >= 0.01 * this->m_FinalTemperature )
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
    }    

  /**
   * Generate output
   */   
  typedef BSplineControlPointImageFilter
    <ControlPointLatticeType, VectorFieldType> BSplineControlPointsFilterType;
  typename BSplineControlPointsFilterType::Pointer bspliner 
    = BSplineControlPointsFilterType::New();

  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetInput( this->m_ControlPointLattice );
  bspliner->SetOrigin( this->m_Origin );
  bspliner->SetSize( this->m_Size );
  bspliner->SetSpacing( this->m_Spacing );
  bspliner->Update();

  this->GraftOutput( bspliner->GetOutput() );
}      
         
template <class TPointSet, class TOutputImage>
void
BSplineRobustPointMethodPointSetFilter<TPointSet, TOutputImage>
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
  this->m_CorrespondenceMatrix.SetSize( 
    this->m_VPoints->GetNumberOfPoints(), this->GetInput( 0 )->GetNumberOfPoints() );
  this->m_OutlierColumn.SetSize( this->m_VPoints->GetNumberOfPoints() );
  this->m_OutlierRow.SetSize( this->GetInput( 0 )->GetNumberOfPoints() );

  RealType K = static_cast<RealType>( this->m_VPoints->GetNumberOfPoints() );
  RealType N = static_cast<RealType>( this->GetInput( 0 )->GetNumberOfPoints() );

  this->m_CorrespondenceMatrix.Fill( 1.0 / ( N * K ) );

  this->m_OutlierColumn.Fill( 1.0 / ( 1000 * N * K ) );
  this->m_OutlierRow.Fill( 1.0 / ( 1000 * N * K ) );


  /**
   * Initialize the outlier center of mass
   */ 
/*
  this->m_OutlierPointV.Fill( 0 );
  this->m_OutlierPointX.Fill( 0 );

  for ( unsigned int i = 0; i < this->m_VPoints->GetNumberOfPoints(); i++ )
    {
    typename InputPointSetType::PointType pointV;
    this->m_VPoints->GetPoint( i, &pointV );    
    for ( unsigned int d = 0; d < Dimension; d++ )
      {
      this->m_OutlierPointV[d] += ( pointV[d] / K );
      }
    }  
  for ( unsigned int d = 0; d < Dimension; d++ )
    {
    this->m_OutlierPointV[d] /= K;
    }

  for ( unsigned int j = 0; j < this->GetInput( 0 )->GetNumberOfPoints(); j++ )
    {
    typename InputPointSetType::PointType pointX;
    this->GetInput( 0 )->GetPoint( j, &pointX );    
    for ( unsigned int d = 0; d < Dimension; d++ )
      {
      this->m_OutlierPointX[d] += ( pointX[d] );
      }
    }  
  for ( unsigned int d = 0; d < Dimension; d++ )
    {
    this->m_OutlierPointX[d] /= N;
    }
*/
  this->NormalizeCorrespondenceMatrix();

  this->UpdateTransformation();

  if ( this->m_FinalTemperature == NumericTraits<RealType>::max() && 
       this->m_InitialTemperature == NumericTraits<RealType>::max() )
    {
    this->CalculateInitialAndFinalTemperatures();
    }
}

template <class TPointSet, class TOutputImage>
void
BSplineRobustPointMethodPointSetFilter<TPointSet, TOutputImage>
::CalculateInitialAndFinalTemperatures()
{
  RealType maxSquaredDistance = 0;
  RealType sumMinSquaredDistance = 0;

  for ( unsigned int i = 0; i < this->m_VPoints->GetNumberOfPoints(); i++ )
    {
    typename InputPointSetType::PointType pointV;
    this->m_VPoints->GetPoint( i, &pointV );    

    RealType minSquaredDistance = NumericTraits<RealType>::max();  
  
    for ( unsigned int j = 0; j < this->GetInput( 0 )->GetNumberOfPoints(); j++ )
      {
      typename InputPointSetType::PointType pointX;      
      this->GetInput( 0 )->GetPoint( j, &pointX );      
      
      RealType squaredDistance = ( pointX - pointV ).GetSquaredNorm();

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
    / static_cast<RealType>( this->m_VPoints->GetNumberOfPoints() );
}

template <class TPointSet, class TOutputImage>
void
BSplineRobustPointMethodPointSetFilter<TPointSet, TOutputImage>
::UpdateCorrespondenceMatrix()
{
  for ( unsigned int i = 0; i < this->m_VPoints->GetNumberOfPoints(); i++ )
    {
    typename InputPointSetType::PointType pointV;
    this->m_VPoints->GetPoint( i, &pointV );    

    for ( unsigned int j = 0; j < this->GetInput( 0 )->GetNumberOfPoints(); j++ )
      {
      typename InputPointSetType::PointType pointX;      
      this->GetInput( 0 )->GetPoint( j, &pointX );      
 
      this->m_CorrespondenceMatrix( i, j )  
        = vcl_exp( -0.5 * ( pointX - pointV ).GetSquaredNorm() / this->m_CurrentTemperature )
          / this->m_CurrentTemperature;

      }  
    } 

/*
  this->m_OutlierColumn.Fill( 0 );
  this->m_OutlierRow.Fill( 0 );


  for ( unsigned int i = 0; i < this->m_VPoints->GetNumberOfPoints(); i++ )
    {
    typename InputPointSetType::PointType pointV;
    this->m_VPoints->GetPoint( i, &pointV );    

    this->m_OutlierColumn[i] 
      = vcl_exp( -0.5 * ( this->m_OutlierPointX - pointV ).GetSquaredNorm() 
          / this->m_InitialTemperature )/ this->m_InitialTemperature;
    } 

  for ( unsigned int j = 0; j < this->GetInput( 0 )->GetNumberOfPoints(); j++ )
    {
    typename InputPointSetType::PointType pointX;
    this->GetInput( 0 )->GetPoint( j, &pointX );    

    this->m_OutlierRow[j] 
      = vcl_exp( -0.5 * ( this->m_OutlierPointV - pointX ).GetSquaredNorm() 
          / this->m_InitialTemperature )/ this->m_InitialTemperature;
    } 
*/
  this->NormalizeCorrespondenceMatrix();
}

template <class TPointSet, class TOutputImage>
void
BSplineRobustPointMethodPointSetFilter<TPointSet, TOutputImage>
::NormalizeCorrespondenceMatrix()
{
  RealType epsilon = 1e-4;
  RealType deviation = NumericTraits<RealType>::max(); 

  OutlierVectorType rowSum;
  OutlierVectorType columnSum;

  unsigned int maximumNumberOfIterations = 100;
  unsigned int iterations = 0;

  while ( vnl_math_abs( deviation ) > epsilon && iterations++ < maximumNumberOfIterations )
    {
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
BSplineRobustPointMethodPointSetFilter<TPointSet, TOutputImage>
::UpdateTransformation()
{
  typename PointSetType::Pointer points = PointSetType::New();
  points->Initialize();
  typename WeightsContainerType::Pointer weights = WeightsContainerType::New();
 
  unsigned int count = 0;
  for ( unsigned int i = 0; i < this->m_VPoints->GetNumberOfPoints(); i++ )
    {
    typename InputPointSetType::PointType V;
    this->m_VPoints->GetPoint( i, &V );

    if ( this->m_SolveSimplerLeastSquaresProblem )
      { 
      typename InputPointSetType::PointType Y;
      Y.Fill( 0 );
      RealType weight = 0;
      for ( unsigned int j = 0; j < this->GetInput( 0 )->GetNumberOfPoints(); j++ ) 
        {
        typename InputPointSetType::PointType X;
        this->GetInput( 0 )->GetPoint( j, &X );
        for ( unsigned int d = 0; d < Dimension; d++ )
          {
          Y[d] += ( X[d] * this->m_CorrespondenceMatrix( i, j ) ); 
          }
        
        weight += this->m_CorrespondenceMatrix( i, j );
        }
      VectorType vector = Y - V;

      points->SetPoint( count, V );  
      points->SetPointData( count, vector );
      weights->InsertElement( count, 1.0 );
      count++;
      }
    else
      {
      for ( unsigned int j = 0; j < this->GetInput( 0 )->GetNumberOfPoints(); j++ ) 
        {
        if ( this->m_CorrespondenceMatrix( i, j ) <= 0 )
          {
          continue;
          } 
        typename InputPointSetType::PointType X;
        this->GetInput( 0 )->GetPoint( j, &X );

        VectorType vector = X - V;

        points->SetPoint( count, V );  
        points->SetPointData( count, vector );
        weights->InsertElement( count, this->m_CorrespondenceMatrix( i, j ) );
        count++;
        }
      }     
    } 

  typename BSplineFitterType::Pointer fitter = BSplineFitterType::New();

  unsigned int nLevels = this->m_NumberOfLevels;
  if ( this->m_CurrentTemperature == NumericTraits<RealType>::max() )
    {
    nLevels = 1;
    }

//  fitter->DebugOn();
  fitter->SetInput( points );
  fitter->SetPointWeights( weights );
  fitter->SetOrigin( this->m_Origin );
  fitter->SetSpacing( this->m_Spacing );
  fitter->SetSize( this->m_Size );
  fitter->SetNumberOfLevels( nLevels );              
  fitter->SetSplineOrder( this->m_SplineOrder );
  fitter->SetNumberOfControlPoints( this->m_NumberOfControlPoints );       
  fitter->SetGenerateOutputImage( false );
  fitter->Update();

  /**
   * Update total transformation defined by the control point lattice.
   */
  if ( !this->m_ControlPointLattice )
    {
    typedef BSplineControlPointImageFilter
      <ControlPointLatticeType, VectorFieldType> BSplineControlPointsFilterType;
    typename BSplineControlPointsFilterType::Pointer bspliner = BSplineControlPointsFilterType::New();
  
    bspliner->SetInput( fitter->GetPhiLattice() );
    bspliner->SetSplineOrder( this->m_SplineOrder );
    bspliner->SetOrigin( this->m_Origin );
    bspliner->SetSpacing( this->m_Spacing );
    bspliner->SetSize( this->m_Size );

    typename BSplineControlPointsFilterType::ArrayType nLevels;
    nLevels.Fill( this->m_NumberOfLevels );  
    this->m_ControlPointLattice = bspliner->RefineControlLattice( nLevels );
    }
  else
    {
    ImageRegionIterator<ControlPointLatticeType> ItC( 
      this->m_ControlPointLattice, 
      this->m_ControlPointLattice->GetLargestPossibleRegion() );
    ImageRegionIterator<ControlPointLatticeType> ItP( 
      fitter->GetPhiLattice(), fitter->GetPhiLattice()->GetLargestPossibleRegion() );
    for ( ItC.GoToBegin(), ItP.GoToBegin(); !ItC.IsAtEnd(); ++ItC, ++ItP )
      {
      ItC.Set( ItC.Get() + ItP.Get() );
      }
    }  

  /**
   *  As a check, ensure that the total distance between points is decreasing
   */

  RealType error = 0.0;

  for ( unsigned int i = 0; i < this->m_VPoints->GetNumberOfPoints(); i++ )
    {
    typename InputPointSetType::PointType V;
    this->m_VPoints->GetPoint( i, &V );

    for ( unsigned int j = 0; j < this->GetInput( 0 )->GetNumberOfPoints(); j++ ) 
      {
      typename InputPointSetType::PointType X;
      this->GetInput( 0 )->GetPoint( j, &X );

      error += ( this->m_CorrespondenceMatrix( i, j ) * ( X - V ).GetNorm() );  
      }
    }     

  error /= static_cast<RealType>( this->m_VPoints->GetNumberOfPoints() * this->GetInput( 0 )->GetNumberOfPoints() );

  std::cout << "Before error = " << error << std::endl;  

  /**
   * Update the V points
   */
  for ( unsigned int i = 0; i < this->m_VPoints->GetNumberOfPoints(); i++ )
    {
    typename InputPointSetType::PointType V;
    this->m_VPoints->GetPoint( i, &V );

    typename VectorFieldType::PixelType vector;
    typename PointSetType::PointType point;
    for ( unsigned int d = 0; d < Dimension; d++ )
      {
      point[d] = V[d];
      }       
    fitter->EvaluateAtPoint( point, vector );
    V += vector;
    this->m_VPoints->SetPoint( i, V );
    }      

  /**
   *  As a check, ensure that the total distance between points is decreasing
   */

  for ( unsigned int i = 0; i < this->m_VPoints->GetNumberOfPoints(); i++ )
    {
    typename InputPointSetType::PointType V;
    this->m_VPoints->GetPoint( i, &V );

    for ( unsigned int j = 0; j < this->GetInput( 0 )->GetNumberOfPoints(); j++ ) 
      {
      typename InputPointSetType::PointType X;
      this->GetInput( 0 )->GetPoint( j, &X );

      error += ( this->m_CorrespondenceMatrix(i, j) * ( X - V ).GetNorm() );  
      }
    }     

  error /= static_cast<RealType>( this->m_VPoints->GetNumberOfPoints() * this->GetInput( 0 )->GetNumberOfPoints() );

  std::cout << "After error = " << error << std::endl;  


}

template <class TPointSet, class TOutputImage>
void
BSplineRobustPointMethodPointSetFilter<TPointSet, TOutputImage>
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


}



/**
 * Standard "PrintSelf" method
 */
template <class TPointSet, class TOutputImage>
void
BSplineRobustPointMethodPointSetFilter<TPointSet, TOutputImage>
::PrintSelf(
  std::ostream& os, 
  Indent indent) const
{
}



}  //end namespace itk

#endif

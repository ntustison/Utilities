#ifndef __itkBSplineGridTaggingImageFilter_hxx
#define __itkBSplineGridTaggingImageFilter_hxx

#include "itkBSplineGridTaggingImageFilter.h"

#include "itkBSplineControlPointImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"

#include "vnl/vnl_math.h"

namespace itk
{

template <class TControlPointLattice, class TCandidatePointImage, class TOutputImage>
BSplineGridTaggingImageFilter<TControlPointLattice, TCandidatePointImage, TOutputImage>
::BSplineGridTaggingImageFilter()
{
  this->SetNumberOfRequiredInputs( 2 );

  this->m_SplineOrder = 3;
  this->m_NumberOfLevels = 1;
  this->m_NumberOfControlPoints.Fill( 4 );

  this->m_FinalTemperature = 0.0;
  this->m_AnnealingRate = 0.9;
  this->m_NumberOfIterationsPerTemperature = 2;

  this->m_SolveSimplerLeastSquaresProblem = true;
}

template <class TControlPointLattice, class TCandidatePointImage, class TOutputImage>
BSplineGridTaggingImageFilter<TControlPointLattice, TCandidatePointImage, TOutputImage>
::~BSplineGridTaggingImageFilter()
{  
}

template <class TControlPointLattice, class TCandidatePointImage, class TOutputImage>
void
BSplineGridTaggingImageFilter<TControlPointLattice, TCandidatePointImage, TOutputImage>
::SetInput1( const TControlPointLattice * image1 ) 
{
  // Process object is not const-correct so the const casting is required.
  SetNthInput(0, const_cast<TControlPointLattice *>( image1 ));
}

template <class TControlPointLattice, class TCandidatePointImage, class TOutputImage>
void
BSplineGridTaggingImageFilter<TControlPointLattice, TCandidatePointImage, TOutputImage>
::SetInput2( const TCandidatePointImage * image2 ) 
{
  // Process object is not const-correct so the const casting is required.
  SetNthInput(1, const_cast<TCandidatePointImage *>( image2 ));
}

template <class TControlPointLattice, class TCandidatePointImage, class TOutputImage>
void
BSplineGridTaggingImageFilter<TControlPointLattice, TCandidatePointImage, TOutputImage>
::GenerateData()
{

  this->InitializeLabelPointSets();

  this->m_CurrentTemperature = NumericTraits<RealType>::max();

  while ( 0 )
//  while ( this->m_CurrentTemperature > this->m_FinalTemperature )
    {
    for ( unsigned int i = 0; i < this->m_NumberOfIterationsPerTemperature; i++ )
      {
      /**
       * Step 1:  Update the correspondence matrix
       */
      this->UpdateCorrespondenceMatrix();
  
      if ( this->m_CurrentTemperature == NumericTraits<RealType>::max() )
        {
        this->m_CurrentTemperature = this->m_InitialTemperature;
        } 
      /**
       * Step 2:  Update the transformation
       */
      this->UpdateTransformation();   
      } 
    this->m_CurrentTemperature *= this->m_AnnealingRate;
    }    


  typedef ImageDuplicator<OutputImageType> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage( this->GetInput( 0 ) );
  duplicator->Update();
  
  ImageRegionIterator<ScalarControlPointLatticeType> ItS( this->m_ScalarControlPointLattice,
    this->m_ScalarControlPointLattice->GetLargestPossibleRegion() );
  ImageRegionIterator<OutputImageType> ItD( duplicator->GetOutput(),
    duplicator->GetOutput()->GetLargestPossibleRegion() );
  for ( ItS.GoToBegin(), ItD.GoToBegin(); !ItS.IsAtEnd(); ++ItS, ++ItD )
    {
    typename OutputImageType::PixelType pixel = ItD.Get(); 
    pixel[this->m_WhichParametricDimension] = ItS.Get()[0];
    ItD.Set( pixel ); 
    }  

  this->GraftOutput( duplicator->GetOutput() );
}      
         
template <class TControlPointLattice, class TCandidatePointImage, class TOutputImage>
void
BSplineGridTaggingImageFilter<TControlPointLattice, TCandidatePointImage, TOutputImage>
::InitializeLabelPointSets()
{
  typename ControlPointLatticeType::ConstPointer initialControlPointLattice
    = dynamic_cast<const TControlPointLattice*>( ProcessObject::GetInput( 0 ) );
  typename CandidatePointImageType::ConstPointer candidatePointImage
    = dynamic_cast<const TCandidatePointImage*>( ProcessObject::GetInput( 1 ) );

  /**
   * Initialize control point lattice
   */

  this->m_ScalarControlPointLattice = ScalarControlPointLatticeType::New();
  this->m_ScalarControlPointLattice->SetOrigin( initialControlPointLattice->GetOrigin() );  
  this->m_ScalarControlPointLattice->SetSpacing( initialControlPointLattice->GetSpacing() );  
  this->m_ScalarControlPointLattice->SetRegions( initialControlPointLattice->GetLargestPossibleRegion() );  
  this->m_ScalarControlPointLattice->Allocate();
  ScalarType S;
  S.Fill( 0 );
  this->m_ScalarControlPointLattice->FillBuffer( S );

  ImageRegionIterator<ScalarControlPointLatticeType> ItS( this->m_ScalarControlPointLattice,
    this->m_ScalarControlPointLattice->GetLargestPossibleRegion() );
  ImageRegionConstIterator<ControlPointLatticeType> ItC( initialControlPointLattice,
    initialControlPointLattice->GetLargestPossibleRegion() );
  for ( ItS.GoToBegin(), ItC.GoToBegin(); !ItS.IsAtEnd(); ++ItS, ++ItC )
    {
    ItS.Set( ItC.Get()[this->m_WhichParametricDimension] ); 
    }  

  /**
   * Create the V label image from the input control point lattice
   */
  typename CandidatePointImageType::Pointer VLabelImage = CandidatePointImageType::New();
  VLabelImage->SetOrigin( candidatePointImage->GetOrigin() );
  VLabelImage->SetSpacing( candidatePointImage->GetSpacing() );
  VLabelImage->SetRegions( candidatePointImage->GetLargestPossibleRegion() );
  VLabelImage->Allocate();
  VLabelImage->FillBuffer( 0 ); 

  typename ScalarControlPointLatticeType::PointType origin;
  typename ScalarControlPointLatticeType::SizeType size;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    origin[i] = this->m_BoundingBoxMinimum[i];
    size[i] = static_cast<unsigned int>(
      ( candidatePointImage->GetSpacing()[i] +
        this->m_BoundingBoxMaximum[i] - this->m_BoundingBoxMinimum[i] )
      / candidatePointImage->GetSpacing()[i] );  
    }
  typename CandidatePointImageType::IndexType originIndex;
  candidatePointImage->TransformPhysicalPointToIndex( origin, originIndex ); 

  typedef BSplineControlPointImageFilter
    <ScalarControlPointLatticeType, ScalarFieldType> BSplineControlPointsFilterType;
  typename BSplineControlPointsFilterType::Pointer bspliner 
    = BSplineControlPointsFilterType::New();

  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetInput( this->m_ScalarControlPointLattice );
  bspliner->SetOrigin( origin );
  bspliner->SetSize( size );
  bspliner->SetSpacing( candidatePointImage->GetSpacing() );
 // bspliner->Update();

  for ( unsigned int n = 0; n < this->m_NumberOfTagLabels; n++ )
    {
    typename BSplineControlPointsFilterType::PointType parametricPoint;
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      if ( i != this->m_WhichParametricDimension )
        { 
        parametricPoint[i] = 1e10;
        }
      else
        {
        parametricPoint[i] = this->m_FirstTagLocation + 
          this->m_InitialTagSpacing * static_cast<RealType>( n ) ;
        }    
      } 
    typename ScalarFieldType::Pointer tagPlane = bspliner->GenerateOutputImageAt( parametricPoint ); 

    ImageRegionIteratorWithIndex<ScalarFieldType> It( tagPlane, 
      tagPlane->GetLargestPossibleRegion() );
    for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      typename ScalarFieldType::IndexType index = It.GetIndex();
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        if ( i != this->m_WhichParametricDimension )
          {
          index[i] += originIndex[i];
          }
        }
      typename CandidatePointImageType::PointType pt;
      candidatePointImage->TransformIndexToPhysicalPoint( index, pt );
      pt[this->m_WhichParametricDimension] = It.Get()[0] 
        + ( this->m_FirstTagLocation + this->m_InitialTagSpacing * ( n ) ); 
      candidatePointImage->TransformPhysicalPointToIndex( pt, index );
      
      if ( VLabelImage->GetLargestPossibleRegion().IsInside( index ) )
        {  
        VLabelImage->SetPixel( index, static_cast<typename LabelImageType::PixelType>( n+1 ) );   
        }
      }
    }  

  {
  /**
   * Label each V point that's in the bounding box.
   */

  this->m_VPoints = LabelPointSetType::New();
  this->m_VPoints->Initialize();

  ImageRegionIteratorWithIndex<LabelImageType> ItV( VLabelImage, 
    VLabelImage->GetLargestPossibleRegion() );

  unsigned int K = 0;
  for ( ItV.GoToBegin(); !ItV.IsAtEnd(); ++ItV )
    {
    if ( ItV.Get() > 0 )
      {
      typename LabelImageType::PointType pt;
      VLabelImage->TransformIndexToPhysicalPoint( ItV.GetIndex(), pt );
      typename LabelPointSetType::PointType point;
      LabelType label;

      bool isOutside = false;
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        point[i] = pt[i];
        if ( pt[i] <= this->m_BoundingBoxMinimum[i] || 
             pt[i] >= this->m_BoundingBoxMaximum[i] )
          {
          isOutside = true;
          break;
          }
        } 
      if ( !isOutside ) 
        { 
        this->m_VPoints->SetPoint( K, point ); 
        label[0] = ItV.Get();
        this->m_VPoints->SetPointData( K, label ); 
        ItV.Set( K+1 );
        K++;
        } 
      }
    }     
  }   

  /**
   * Label each X point that's in the bounding box.
   */

  typename ImageDuplicator<CandidatePointImageType>::Pointer XLabelImage = 
    ImageDuplicator<CandidatePointImageType>::New();
  XLabelImage->SetInputImage( candidatePointImage );
  XLabelImage->Update(); 

  {
  this->m_XPoints = LabelPointSetType::New();
  this->m_XPoints->Initialize();

  ImageRegionIteratorWithIndex<CandidatePointImageType> ItX( XLabelImage->GetOutput(), 
    XLabelImage->GetOutput()->GetLargestPossibleRegion() );

  unsigned int N = 0;
  for ( ItX.GoToBegin(); !ItX.IsAtEnd(); ++ItX )
    {
    if ( ItX.Get() > 0 )
      {
      typename CandidatePointImageType::PointType pt;
      XLabelImage->GetOutput()->TransformIndexToPhysicalPoint( ItX.GetIndex(), pt );
      typename LabelPointSetType::PointType point;

      bool isOutside = false;
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        point[i] = pt[i];
        if ( pt[i] <= this->m_BoundingBoxMinimum[i] || 
             pt[i] >= this->m_BoundingBoxMaximum[i] )
          {
          isOutside = true;
          break;
          }
        } 
      if ( !isOutside ) 
        {
        this->m_XPoints->SetPoint( N, point ); 
        ItX.Set( N+1 );
        N++;
        } 
      }
    }     
  }

  /**
   * Initialize correspondence matrix
   */

  this->m_CorrespondenceMatrix.set_size( 
    this->m_VPoints->GetNumberOfPoints(), this->m_XPoints->GetNumberOfPoints() );
  this->m_OutlierColumn.set_size( this->m_VPoints->GetNumberOfPoints() );
  this->m_OutlierRow.set_size( this->m_XPoints->GetNumberOfPoints() );

  ImageRegionIteratorWithIndex<LabelImageType> ItV( VLabelImage, 
    VLabelImage->GetLargestPossibleRegion() );
  ImageLinearIteratorWithIndex<LabelImageType> ItX( XLabelImage->GetOutput(), 
    XLabelImage->GetOutput()->GetLargestPossibleRegion() );

  RealType maxSquaredDistance = 0.0;
  RealType sumMinSquaredDistance = 0.0;

  ItX.SetDirection( this->m_WhichParametricDimension );
  for ( ItV.GoToBegin(); !ItV.IsAtEnd(); ++ItV )
    {
    unsigned long labelV = ItV.Get();
    if ( labelV > 0 )
      {
      typename CandidatePointImageType::PointType pointV;
      VLabelImage->TransformIndexToPhysicalPoint( ItV.GetIndex(), pointV );
 
      RealType minSquaredDistance = NumericTraits<RealType>::max();  

      vcl_vector<RealType> values;
      vcl_vector<int> columns;

      ItX.SetIndex( ItV.GetIndex() );
      for ( ItX.GoToBeginOfLine(); !ItX.IsAtEndOfLine(); ++ItX ) 
        {
        unsigned long labelX = ItX.Get();
        if ( labelX > 0 )
          {
          typename LabelImageType::PointType pointX;
          XLabelImage->GetOutput()->TransformIndexToPhysicalPoint( ItX.GetIndex(), pointX );
           
          RealType squaredDistance = ( pointX - pointV ).GetSquaredNorm();
          values.push_back( vcl_exp( -0.5 * squaredDistance ) );
          columns.push_back( labelX-1 );

          if ( squaredDistance > maxSquaredDistance )
            {
            maxSquaredDistance = squaredDistance; 
            } 
          if ( squaredDistance < minSquaredDistance )
            {
            minSquaredDistance = squaredDistance; 
            } 
          }  
        }
      this->m_CorrespondenceMatrix.set_row( labelV-1, columns, values ); 
      if ( minSquaredDistance != NumericTraits<RealType>::max() )
        {
        sumMinSquaredDistance += minSquaredDistance;
        }
      }
    } 

  RealType scaleFactor = 0;
  this->m_InitialTemperature = maxSquaredDistance;
  this->m_FinalTemperature = sumMinSquaredDistance 
    / static_cast<RealType>( this->m_CorrespondenceMatrix.rows() );
  scaleFactor = vcl_exp( -1.0 / this->m_InitialTemperature ) / this->m_InitialTemperature;

  for ( unsigned int i = 0; i < this->m_CorrespondenceMatrix.rows(); i++ )
    {
    this->m_CorrespondenceMatrix.scale_row( i, scaleFactor );  
    }  
}

template <class TControlPointLattice, class TCandidatePointImage, class TOutputImage>
void
BSplineGridTaggingImageFilter<TControlPointLattice, TCandidatePointImage, TOutputImage>
::UpdateCorrespondenceMatrix()
{
  typename LabelPointSetType::PointType pointX;
  typename LabelPointSetType::PointType pointV;

  RealType N = static_cast<RealType>( this->m_VPoints->GetNumberOfPoints() );
  RealType K = static_cast<RealType>( this->m_XPoints->GetNumberOfPoints() );

  typedef typename Statistics
     ::MersenneTwisterRandomVariateGenerator GeneratorType; 
  typename GeneratorType::Pointer generator = GeneratorType::New();
  generator->SetSeed();

  for ( unsigned int j = 0; j < this->m_CorrespondenceMatrix.rows(); j++ )
    {
    this->m_VPoints->GetPoint( j, &pointV ); 
    typename SparseMatrixType::row &rw = this->m_CorrespondenceMatrix.get_row( j );
    typename SparseMatrixType::row::iterator ri;
    for ( ri = rw.begin(); ri != rw.end(); ++ri )
      {
      this->m_XPoints->GetPoint( (*ri).first, &pointX );
      (*ri).second = vcl_exp( -( pointX - pointV ).GetSquaredNorm() 
        / this->m_CurrentTemperature ) + generator->GetNormalVariate( 0.0, 0.001 / N ); 
      } 
    }

  this->m_OutlierRow.fill( 0.1 / N );
  this->m_OutlierColumn.fill( 0.1 / K );

  this->NormalizeCorrespondenceMatrix();
}

template <class TControlPointLattice, class TCandidatePointImage, class TOutputImage>
void
BSplineGridTaggingImageFilter<TControlPointLattice, TCandidatePointImage, TOutputImage>
::NormalizeCorrespondenceMatrix()
{
  RealType epsilon = 0.0025;
  RealType deviation = NumericTraits<RealType>::max(); 

  unsigned int maximumNumberOfIterations = 100;
  unsigned int iterations = 0;

  while ( vnl_math_abs( deviation ) > epsilon && iterations++ < maximumNumberOfIterations )
    {
    /**
     * do row normalization
     */
    OutlierArrayType rowSum;
    rowSum = this->m_OutlierColumn;

    for ( unsigned int i = 0; i < this->m_CorrespondenceMatrix.rows(); i++ )
      {
      rowSum[i] += this->m_CorrespondenceMatrix.sum_row( i );
      } 
    for ( unsigned int i = 0; i < this->m_CorrespondenceMatrix.rows(); i++ )
      {
      this->m_CorrespondenceMatrix.scale_row( i, 1.0 / rowSum[i] );
      }
    for ( unsigned int i = 0; i < this->m_OutlierColumn.size(); i++ )
      {
      this->m_OutlierColumn[i] /= rowSum[i];     
      }

    /**
     * do column normalization
     */
    OutlierArrayType columnSum;
    columnSum = this->m_OutlierRow;

    for ( unsigned int i = 0; i < this->m_CorrespondenceMatrix.rows(); i++ )
      {
      typename SparseMatrixType::row &rw = this->m_CorrespondenceMatrix.get_row( i );
      typename SparseMatrixType::row::iterator ri;
      for ( ri = rw.begin(); ri != rw.end(); ++ri )
        { 
        columnSum[(*ri).first] += (*ri).second;
        }
      } 
    for ( unsigned int i = 0; i < this->m_CorrespondenceMatrix.rows(); i++ )
      {
      typename SparseMatrixType::row &rw = this->m_CorrespondenceMatrix.get_row( i );
      typename SparseMatrixType::row::iterator ri;
      for ( ri = rw.begin(); ri != rw.end(); ++ri )
        { 
        (*ri).second /= columnSum[(*ri).first];
        }
      }
    for ( unsigned int i = 0; i < this->m_OutlierRow.size(); i++ )
      {
      this->m_OutlierRow[i] /= columnSum[i];     
      }

    /**
     * Calculate current deviation from 1
     */
    for ( unsigned int i = 0; i < rowSum.size(); i++ )
      {
      rowSum[i] -= 1.0;     
      }
    for ( unsigned int i = 0; i < columnSum.size(); i++ )
      {
      columnSum[i] -= 1.0;     
      }

    deviation = ( rowSum.squared_magnitude() + columnSum.squared_magnitude() ) 
      / ( static_cast<RealType>( rowSum.size() + columnSum.size() ) );
    }  
}

template <class TControlPointLattice, class TCandidatePointImage, class TOutputImage>
void
BSplineGridTaggingImageFilter<TControlPointLattice, TCandidatePointImage, TOutputImage>
::UpdateTransformation()
{
  typename PointSetType::Pointer points = PointSetType::New();
  points->Initialize();
  typename WeightsContainerType::Pointer weights = WeightsContainerType::New();
 
  unsigned int count = 0;
  for ( unsigned int i = 0; i < this->m_CorrespondenceMatrix.rows(); i++ )
    {
    typename LabelPointSetType::PointType V;
    this->m_VPoints->GetPoint( i, &V );
    LabelType label;        
    this->m_VPoints->GetPointData( i, &label );
    V[this->m_WhichParametricDimension] 
      = this->m_FirstTagLocation + this->m_InitialTagSpacing * ( label[0] - 1 ); 

    typename SparseMatrixType::row &rw = this->m_CorrespondenceMatrix.get_row( i );
    typename SparseMatrixType::row::iterator ri;

    if ( this->m_SolveSimplerLeastSquaresProblem )
      { 
      typename LabelPointSetType::PointType Y;
      Y.Fill( 0 );
      for ( ri = rw.begin(); ri != rw.end(); ++ri ) 
        {
        typename LabelPointSetType::PointType X;
        this->m_XPoints->GetPoint( (*ri).first, &X );
        for ( unsigned int d = 0; d < ImageDimension; d++ )
          {
          Y[d] += ( X[d] * (*ri).second ); 
         }
        }

      ScalarType scalar;
      scalar[0] = Y[this->m_WhichParametricDimension] - V[this->m_WhichParametricDimension];

      points->SetPoint( count, V );  
      points->SetPointData( count, scalar );
      weights->InsertElement( count, 1.0 );
      count++;
      }
    else
      {
      for ( ri = rw.begin(); ri != rw.end(); ++ri ) 
        {
        typename LabelPointSetType::PointType X;
        this->m_XPoints->GetPoint( (*ri).first, &X );

        ScalarType scalar;
        scalar[0] = X[this->m_WhichParametricDimension] - V[this->m_WhichParametricDimension];

        points->SetPoint( count, V );  
        points->SetPointData( count, scalar );
        weights->InsertElement( count, (*ri).second );
        count++;
        }
      }     
    } 

  typename BSplineFitterType::Pointer fitter = BSplineFitterType::New();
  
  typename BSplineFitterType::ArrayType close;
  close.Fill( false );
  typename ScalarControlPointLatticeType::PointType origin;
  typename ScalarControlPointLatticeType::SizeType size;
  typename ScalarControlPointLatticeType::SpacingType spacing;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    origin[i] = this->m_BoundingBoxMinimum[i];
    spacing[i] = this->GetInput( 1 )->GetSpacing()[i]; 
    size[i] = static_cast<unsigned int>(
      ( spacing[i] + this->m_BoundingBoxMaximum[i] - this->m_BoundingBoxMinimum[i] )
      / spacing[i] );  
    spacing[i] = ( this->m_BoundingBoxMaximum[i] - this->m_BoundingBoxMinimum[i] )
      / static_cast<RealType>( size[i] - 1 );
    }

//   std::cout << " (number of points = " << points->GetNumberOfPoints() << ")  " << std::endl;

//  fitter->DebugOn();
  fitter->SetInput( points );
  fitter->SetPointWeights( weights );
  fitter->SetOrigin( origin );
  fitter->SetSpacing( spacing );
  fitter->SetSize( size );
  fitter->SetNumberOfLevels( this->m_NumberOfLevels );              
  fitter->SetSplineOrder( this->m_SplineOrder );
  fitter->SetNumberOfControlPoints( this->m_NumberOfControlPoints );       
  fitter->SetCloseDimension( close );
  fitter->SetGenerateOutputImage( false );
  fitter->Update();

  this->m_ScalarControlPointLattice = fitter->GetPhiLattice();

  ImageRegionIterator<ScalarControlPointLatticeType> ItS( 
    this->m_ScalarControlPointLattice, 
    this->m_ScalarControlPointLattice->GetLargestPossibleRegion() );
  ImageRegionIterator<ScalarControlPointLatticeType> ItP( 
    fitter->GetPhiLattice(), fitter->GetPhiLattice()->GetLargestPossibleRegion() );
  for ( ItS.GoToBegin(), ItP.GoToBegin(); !ItS.IsAtEnd(); ++ItS, ++ItP )
    {
    ItS.Set( ItS.Get() + ItP.Get() );
    }

  typedef BSplineControlPointImageFilter
    <ScalarControlPointLatticeType, ScalarFieldType> BSplineControlPointsFilterType;
  typename BSplineControlPointsFilterType::Pointer bspliner 
    = BSplineControlPointsFilterType::New();

  bspliner->SetSplineOrder( this->m_SplineOrder );
  bspliner->SetCloseDimension( close );
  bspliner->SetInput( this->m_ScalarControlPointLattice );
  bspliner->SetOrigin( origin );
  bspliner->SetSize( size );
  bspliner->SetSpacing( spacing );

  /**
   * Update the V points
   */
  for ( unsigned int i = 0; i < this->m_VPoints->GetNumberOfPoints(); i++ )
    {
    typename LabelPointSetType::PointType V;
    this->m_VPoints->GetPoint( i, &V );
    LabelType label;        
    this->m_VPoints->GetPointData( i, &label );
    V[this->m_WhichParametricDimension] 
      = ( this->m_FirstTagLocation + this->m_InitialTagSpacing * ( label[0] - 1 ) );
    typename ScalarFieldType::PixelType scalar;
    typename ScalarFieldType::PointType point;
    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      point[d] = V[d];
      }       
    bspliner->EvaluateAtPoint( point, scalar );
    V[this->m_WhichParametricDimension] += scalar[0];
    this->m_VPoints->SetPoint( i, V );
    }      
}

/**
 * Standard "PrintSelf" method
 */
template <class TControlPointLattice, class TCandidatePointImage, class TOutputImage>
void
BSplineGridTaggingImageFilter<TControlPointLattice, TCandidatePointImage, TOutputImage>
::PrintSelf(
  std::ostream& os, 
  Indent indent) const
{
}



}  //end namespace itk

#endif

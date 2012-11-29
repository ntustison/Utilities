#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryThinning3DImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkPointSet.h"

#include <fstream>

#include "vnl/vnl_cross.h"
#include "vcl_cmath.h"

template <unsigned int ImageDimension>
int Tortuosity( unsigned int argc, char *argv[] )
{

  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<unsigned int, ImageDimension> LabelImageType;

  typedef itk::ImageFileReader<LabelImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::Vector<float, ImageDimension> VectorType;
  typedef itk::Image<VectorType, 1> CurveImageType;

  typedef itk::PointSet<VectorType, 1> PointSetType;

  typename PointSetType::PixelType startPoint( 0.0 );
  typename PointSetType::PixelType endPoint( 0.0 );
  float startPointCount = 0;
  float endPointCount = 0;

  itk::ImageRegionConstIteratorWithIndex<LabelImageType> It( reader->GetOutput(),
    reader->GetOutput()->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if( It.Get() == 1 )
      {
      typename LabelImageType::PointType labelPoint;
      reader->GetOutput()->TransformIndexToPhysicalPoint( It.GetIndex(), labelPoint );

      startPoint += labelPoint.GetVectorFromOrigin();
      startPointCount++;
      }
    else if( It.Get() == 2 )
      {
      typename LabelImageType::PointType labelPoint;
      reader->GetOutput()->TransformIndexToPhysicalPoint( It.GetIndex(), labelPoint );

      endPoint += labelPoint.GetVectorFromOrigin();
      endPointCount++;
      }
    }
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    startPoint[d] /= startPointCount;
    endPoint[d] /= endPointCount;
    }

  // Extract the vessel and thin

  typedef itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType> ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( reader->GetOutput() );
  thresholder->SetLowerThreshold( 3 );
  thresholder->SetUpperThreshold( 3 );
  thresholder->SetInsideValue( 1 );
  thresholder->SetOutsideValue( 0 );

  typedef itk::BinaryThinning3DImageFilter<LabelImageType, LabelImageType> ThinnerType;
  typename ThinnerType::Pointer thinner = ThinnerType::New();
  thinner->SetInput( thresholder->GetOutput() );
  thinner->Update();

  // Extract the center points in order

  typename PointSetType::Pointer pointSet = PointSetType::New();
  pointSet->Initialize();

  typename PointSetType::PointType startParam;
  startParam[0] = 0;
  pointSet->SetPoint( 0, startParam );
  pointSet->SetPointData( 0, startPoint );

  typename PointSetType::PixelType currentPoint = startPoint;
  unsigned long pointCount = 1;

  bool morePoints = true;
  while( morePoints )
    {
    morePoints = false;
    float minDistance = itk::NumericTraits<float>::max();
    typename PointSetType::PixelType minPoint = currentPoint;
    typename LabelImageType::IndexType minIndex;
    minIndex = thinner->GetOutput()->GetLargestPossibleRegion().GetIndex();

    itk::ImageRegionIteratorWithIndex<LabelImageType> ItT( thinner->GetOutput(),
      thinner->GetOutput()->GetLargestPossibleRegion() );
    for( ItT.GoToBegin(); !ItT.IsAtEnd(); ++ItT )
      {
      if( ItT.Get() == 1 )
        {
        typename LabelImageType::PointType labelPoint;
        thinner->GetOutput()->TransformIndexToPhysicalPoint( ItT.GetIndex(), labelPoint );

        float distance = ( labelPoint.GetVectorFromOrigin() - currentPoint ).GetSquaredNorm();

        if( distance <= minDistance )
          {
          minDistance = distance;
          minPoint = labelPoint.GetVectorFromOrigin();
          minIndex = ItT.GetIndex();
          }
        morePoints = true;
        }
      }

    if( morePoints )
      {
      thinner->GetOutput()->SetPixel( minIndex, 2 );
      typename PointSetType::PointType currentParam;
      currentParam[0] = static_cast<float>( pointCount );
      currentPoint = minPoint;
      pointSet->SetPoint( pointCount, currentParam );
      pointSet->SetPointData( pointCount, currentPoint );
      pointCount++;
      }
    }

  typename PointSetType::PointType currentParam;
  currentParam[0] = static_cast<float>( pointCount );
  pointSet->SetPoint( pointCount, currentParam );
  pointSet->SetPointData( pointCount, endPoint );

  // Fit the points to a B-spline curve

  typedef itk::BSplineScatteredDataPointSetToImageFilter
     <PointSetType, CurveImageType>  FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  filter->SetInput( pointSet );
  filter->SetGenerateOutputImage( true );

  typename CurveImageType::PointType origin;
  origin.Fill( 0.0 );
  filter->SetOrigin( origin );

  // 1000 samples
  typename CurveImageType::SizeType size;
  size[0] = 1000;
  filter->SetSize( size );

  typename CurveImageType::SpacingType spacing;
  spacing[0] = pointCount / static_cast<float>( size[0] - 1 );
  filter->SetSpacing( spacing );

  filter->SetSplineOrder( 3 );

  typename FilterType::ArrayType ncps;
  ncps[0] = 5;
  if ( argc > 4 )
    {
    ncps[0] = 3 + atoi( argv[4] );
    }
  filter->SetNumberOfControlPoints( ncps );

  typename FilterType::ArrayType nlevels;
  nlevels[0] = 1;
  if ( argc > 5 )
    {
    nlevels[0] = atoi( argv[5] );
    }
  filter->SetNumberOfLevels( nlevels );

  filter->Update();

  // Spit out the result to a file
  {
  std::string filename = std::string( argv[3] );
  std::ofstream ostr( filename.c_str() );
  ostr << "0 0 0 0" << std::endl;

  itk::ImageRegionIterator<CurveImageType> It(
    filter->GetOutput(), filter->GetOutput()->GetLargestPossibleRegion() );
  float sample = 0;
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    ostr << It.Get()[0] << " " << It.Get()[1] << " " << It.Get()[2] << " " << 1 << std::endl;
    sample++;
    }
  ostr << "0 0 0 0" << std::endl;
  ostr.close();
  }

  // Calculate the DM metric
  float dm = 0.0;
  {
  typedef itk::ConstNeighborhoodIterator<CurveImageType> NeighborhoodIteratorType;
  typename NeighborhoodIteratorType::RadiusType radius;
  typename NeighborhoodIteratorType::OffsetType offset;
  radius.Fill( 1 );
  offset.Fill( 1 );

  NeighborhoodIteratorType ItN( radius, filter->GetOutput(),
    filter->GetOutput()->GetLargestPossibleRegion() );
  for( ItN.GoToBegin(); !ItN.IsAtEnd(); ++ItN )
    {
    if( filter->GetOutput()->GetLargestPossibleRegion().IsInside( ItN.GetIndex( offset ) ) )
      {
      dm += ( ItN.GetCenterPixel() - ItN.GetPixel( offset ) ).GetNorm();
      }
    }
  std::cout << "Distance Metric = " << dm / ( startPoint - endPoint ).GetNorm() << std::endl;
  }

  // Calculate the ICM metric
  {
  typedef itk::ConstNeighborhoodIterator<CurveImageType> NeighborhoodIteratorType;
  typename NeighborhoodIteratorType::RadiusType radius;
  typename NeighborhoodIteratorType::OffsetType poffset;
  typename NeighborhoodIteratorType::OffsetType noffset;
  radius.Fill( 1 );
  poffset.Fill( 1 );
  noffset.Fill( -1 );

  std::vector<VectorType> Ns;

  NeighborhoodIteratorType ItN( radius, filter->GetOutput(),
    filter->GetOutput()->GetLargestPossibleRegion() );
  for( ItN.GoToBegin(); !ItN.IsAtEnd(); ++ItN )
    {
    if( filter->GetOutput()->GetLargestPossibleRegion().IsInside( ItN.GetIndex( poffset ) ) &&
      filter->GetOutput()->GetLargestPossibleRegion().IsInside( ItN.GetIndex( noffset ) ) )
      {
      VectorType T1 = ItN.GetCenterPixel() - ItN.GetPixel( noffset );
      VectorType T2 = ItN.GetPixel( poffset ) - ItN.GetCenterPixel();
      VectorType V = ItN.GetPixel( poffset ) - ItN.GetPixel( noffset );
      VectorType A = T2 - T1;

      VectorType T = V / V.GetNorm();

      VectorType VxA;
      VxA.SetVnlVector( vnl_cross_3d<float>( V.GetVnlVector(), A.GetVnlVector() ) );

      VectorType N;
      N.SetVnlVector( vnl_cross_3d<float>( VxA.GetVnlVector(), V.GetVnlVector() ) );
      N /= N.GetNorm();

      VectorType B;
      B.SetVnlVector( vnl_cross_3d<float>( T.GetVnlVector(), N.GetVnlVector() ) );

      Ns.push_back( N );
      }
    }

  std::vector<float> deltaNMag;
  for( unsigned int i = 1; i < Ns.size(); ++i )
    {
    VectorType deltaN = Ns[i] - Ns[i-1];
    deltaNMag.push_back( deltaN.GetSquaredNorm() );
    }

  unsigned int icmCount = 1;
  for( unsigned int i = 1; i < deltaNMag.size()-1; ++i )
    {
    if( deltaNMag[i] > 1.0 && deltaNMag[i] > deltaNMag[i-1] && deltaNMag[i] >= deltaNMag[i-1] )
      {
      icmCount++;
      }
    }
  std::cout << "Inflection point count = " << icmCount << std::endl;
  }


  // Calculate the SOAM metric
  {
  typedef itk::ConstNeighborhoodIterator<CurveImageType> NeighborhoodIteratorType;
  typename NeighborhoodIteratorType::RadiusType radius;
  typename NeighborhoodIteratorType::OffsetType poffset;
  typename NeighborhoodIteratorType::OffsetType p2offset;
  typename NeighborhoodIteratorType::OffsetType noffset;
  radius.Fill( 3 );
  poffset.Fill( 1 );
  p2offset.Fill( 2 );
  noffset.Fill( -1 );

  float CP = 0.0;

  NeighborhoodIteratorType ItN( radius, filter->GetOutput(),
    filter->GetOutput()->GetLargestPossibleRegion() );
  for( ItN.GoToBegin(); !ItN.IsAtEnd(); ++ItN )
    {
    if( filter->GetOutput()->GetLargestPossibleRegion().IsInside( ItN.GetIndex( p2offset ) ) &&
      filter->GetOutput()->GetLargestPossibleRegion().IsInside( ItN.GetIndex( noffset ) ) )
      {
      VectorType T1 = ItN.GetCenterPixel() - ItN.GetPixel( noffset );
      VectorType T2 = ItN.GetPixel( poffset ) - ItN.GetCenterPixel();
      VectorType T3 = ItN.GetPixel( p2offset ) - ItN.GetPixel( poffset );

      VectorType T1xT2;
      T1xT2.SetVnlVector( vnl_cross_3d( T1.GetVnlVector(), T2.GetVnlVector() ) );
      VectorType T2xT3;
      T2xT3.SetVnlVector( vnl_cross_3d( T2.GetVnlVector(), T3.GetVnlVector() ) );

      VectorType T1norm = T1 / T1.GetNorm();
      VectorType T2norm = T2 / T2.GetNorm();
      VectorType T3norm = T3 / T3.GetNorm();
      VectorType T1xT2norm = T1xT2 / T1xT2.GetNorm();
      VectorType T2xT3norm = T2xT3 / T2xT3.GetNorm();

      float IP = vcl_acos( dot_product( T1norm.GetVnlVector(), T2norm.GetVnlVector() ) );
      float TP = vcl_acos( dot_product( T1xT2norm.GetVnlVector(), T2xT3norm.GetVnlVector() ) );

      CP += vcl_sqrt( vnl_math_sqr( IP ) + vnl_math_sqr( TP ) );
      }
    }

  std::cout << "SOAM = " << CP / dm << std::endl;
  }

  return 0;
}

int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension";
    std::cerr << " inputImage outputFile.txt [numberOfSpans] [numberOfLevels]" << std::endl;

    std::cerr << "Note: It is assumed that the input image has a label 1 "
     << "region at the beginning of the vessel, a vessel label of 3, and "
     << "a final region with label 2" << std::endl;
    return 1;
    }

  switch( atoi( argv[1] ) )
   {
   case 3:
     Tortuosity<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


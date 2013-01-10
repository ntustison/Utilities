#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkPointSet.h"

template <unsigned int ImageDimension>
int CreateSampledBSplineSurface( int argc, char *argv[] )
{
  typedef itk::Image<float, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  unsigned int whichAxis = atoi( argv[4] );

  typename ImageType::Pointer labelImage = reader->GetOutput();

  //
  // Smooth the results by creating a B-spline surface spanning
  // all the slices.
  //

  typedef itk::Vector<float, 1> ScalarPointType;
  typedef itk::PointSet<ScalarPointType, ImageDimension-1> PointSetType;
  typedef itk::Image<ScalarPointType, ImageDimension-1> ScalarSurfaceType;

  typename PointSetType::Pointer points = PointSetType::New();
  points->Initialize();

  unsigned long N = 0;

  typename ImageType::DirectionType direction = labelImage->GetDirection();
  typename ImageType::DirectionType identity;
  identity.SetIdentity();
  labelImage->SetDirection( identity );

  float averageHeight = 0;

  itk::ImageRegionIteratorWithIndex<ImageType> ItL( labelImage,
    labelImage->GetLargestPossibleRegion() );
  for( ItL.GoToBegin(); !ItL.IsAtEnd(); ++ItL )
    {
    if( ItL.Get() == 1 )
      {
      averageHeight += static_cast<float>( ItL.GetIndex()[whichAxis] );
      N++;
      }
    }

  averageHeight /= static_cast<float>( N );

  N = 0;

  for( ItL.GoToBegin(); !ItL.IsAtEnd(); ++ItL )
    {
    if( ItL.Get() == 1 )
      {
      typename PointSetType::PointType point;

      unsigned int count = 0;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        if( d != whichAxis )
          {
          point[count++] = ItL.GetIndex()[d];
          }
        }
      points->SetPoint( N, point );

      ScalarPointType surfacePoint;
      surfacePoint[0] = static_cast<float>( ItL.GetIndex()[whichAxis] ) - averageHeight;
      points->SetPointData( N, surfacePoint );
      N++;
      }
    }

  labelImage->SetDirection( direction );

  typedef itk::BSplineScatteredDataPointSetToImageFilter<PointSetType, ScalarSurfaceType> BSplinerType;

  typename BSplinerType::ArrayType ncps;
  ncps.Fill( 4 );

  typename ScalarSurfaceType::PointType   bsplineOrigin;
  typename ScalarSurfaceType::SizeType    bsplineSize;
  typename ScalarSurfaceType::SpacingType bsplineSpacing;

  unsigned int count = 0;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    if( d != whichAxis )
      {
      bsplineSize[count] = labelImage->GetLargestPossibleRegion().GetSize()[d];
      bsplineOrigin[count] = labelImage->GetLargestPossibleRegion().GetIndex()[d];
      bsplineSpacing[count] = 1;
      count++;
      }
    }

  typename BSplinerType::Pointer bspliner = BSplinerType::New();
  bspliner->SetInput( points );
  bspliner->SetGenerateOutputImage( true );
  bspliner->SetSplineOrder( 3 );
  bspliner->SetNumberOfLevels( atoi( argv[5] ) );
  bspliner->SetNumberOfControlPoints( ncps );
  bspliner->SetSize( bsplineSize );
  bspliner->SetOrigin( bsplineOrigin );
  bspliner->SetSpacing( bsplineSpacing );
  bspliner->Update();

  labelImage->FillBuffer( 0 );

  itk::ImageRegionIteratorWithIndex<ScalarSurfaceType> ItB( bspliner->GetOutput(), bspliner->GetOutput()->GetRequestedRegion() );
  for( ItB.GoToBegin(); !ItB.IsAtEnd(); ++ItB )
    {
    typename ImageType::IndexType index;
    int count = 0;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      if( d != whichAxis )
        {
        index[count++] = ItB.GetIndex()[d];
        }
      else
        {
        index[count++] = static_cast<long>( vcl_floor( ItB.Get()[0] + averageHeight + 0.5 ) );
        }
      }
    labelImage->SetPixel( index, 1 );
    }

  //
  // Fill in the left side of the image with all '1's.
  //

  itk::ImageLinearIteratorWithIndex<ImageType> ItR( labelImage, labelImage->GetLargestPossibleRegion() );
  ItR.SetDirection( whichAxis );

  ItR.GoToBegin();
  while( !ItR.IsAtEnd() )
    {
    bool penDown = true;
    while( !ItR.IsAtEndOfLine() )
      {
      if( ItR.Get() == 1 )
        {
        penDown = false;
        }
      if( penDown )
        {
        ItR.Set( 1 );
        }
      ++ItR;
      }
    ItR.NextLine();
    }


  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( labelImage );
  writer->SetFileName( argv[3] );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if( argc < 6 )
    {
    std::cerr << "Usage: "<< argv[0] << " imageDimension inputImage outputImage whichAxis numberOfLevels" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     CreateSampledBSplineSurface<2>( argc, argv );
     break;
   case 3:
     CreateSampledBSplineSurface<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


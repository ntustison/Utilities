#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkConfidenceConnectedImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkGridImageSource.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkMinimalPathImageFunction.h"
#include "itkMultiplyImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkPathIterator.h"
#include "itkPointSet.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkThresholdImageFilter.h"

#include "itkVector.h"
#include "itkVectorContainer.h"

template <unsigned int ImageDimension>
int SplitHeliumLungs2D( int argc, char *argv[] )
{
  std::cout << "Not implemented." << std::endl;
  return EXIT_FAILURE;
}

template <unsigned int ImageDimension>
int SplitHeliumLungs3D( int argc, char *argv[] )
{
  /**
   * Get the initial segmentation of both lungs.
   * Steps are:
   *  2. Binary Thresholding
   *  3. Check connected components
   *  4. Split lungs and relabel left lung as '2' and right lung as '3'
   */

  typedef itk::Image<unsigned int, ImageDimension> ImageType;
  typedef itk::Image<float, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<RealImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  /**
   * Step 2:  BinaryThresholding.
   */
  typedef itk::ConfidenceConnectedImageFilter<RealImageType, ImageType>
    ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();

  thresholder->SetInput( reader->GetOutput() );

  thresholder->SetNumberOfIterations( 5 );
  thresholder->SetReplaceValue( 1 );
  thresholder->SetInitialNeighborhoodRadius( 5 );

  typename ImageType::IndexType seed;
  seed.Fill( 10 );

  thresholder->SetSeed( seed );
  thresholder->SetMultiplier( 5 );
  thresholder->Update();


  itk::ImageRegionIterator<ImageType> It( thresholder->GetOutput(),
    thresholder->GetOutput()->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    It.Set( 1 - It.Get() );
    }

  /**
   * Step 3:  Get connected components.
   */

  typedef itk::ConnectedComponentImageFilter
    <ImageType, ImageType> ConnectedComponentType;
  typename ConnectedComponentType::Pointer connecter
    = ConnectedComponentType::New();
  connecter->SetInput( thresholder->GetOutput() );

  typedef itk::RelabelComponentImageFilter<ImageType, ImageType> RelabelerType;
  typename RelabelerType::Pointer relabeler = RelabelerType::New();
  relabeler->SetInput( connecter->GetOutput() );
  relabeler->Update();

// 		{
// 		typedef itk::ImageFileWriter<ImageType> WriterType;
// 		typename WriterType::Pointer writer = WriterType::New();
// 		writer->SetInput( relabeler->GetOutput() );
// 		writer->SetFileName( "relabeled.nii.gz" );
// 		writer->Update();
// 		}


  unsigned int numberOfObjects = relabeler->GetNumberOfObjects();
  if( numberOfObjects > 1 )
    {
    unsigned int start = 2;
    if( relabeler->GetSizeOfObjectInPhysicalUnits( 2 ) >
        0.5 * relabeler->GetSizeOfObjectInPhysicalUnits( 1 ) )
      {
      start = 3;
      }

    itk::ImageRegionIterator<ImageType> It( relabeler->GetOutput(),
      relabeler->GetOutput()->GetLargestPossibleRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      if( It.Get() >= start )
        {
        It.Set( 0 );
        }
      if( It.Get() < start && It.Get() > 0 )
        {
        It.Set( 1 );
        }
      }
    }

  /**
   * Step 4:  Split the lungs and relabel.
   */

  typename ImageType::Pointer labelImage = ImageType::New();
  labelImage->SetOrigin( reader->GetOutput()->GetOrigin() );
  labelImage->SetSpacing( reader->GetOutput()->GetSpacing() );
  labelImage->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  labelImage->SetDirection( reader->GetOutput()->GetDirection() );
  labelImage->Allocate();
  labelImage->FillBuffer( 0 );

  typedef itk::Image<unsigned int, 2> SliceType;
  typedef itk::Image<float, 2> RealSliceType;

  typename ImageType::RegionType region;
  typename ImageType::RegionType::SizeType imageSize
    = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
  typename ImageType::IndexType imageIndex
    = reader->GetOutput()->GetLargestPossibleRegion().GetIndex();

  typename ImageType::IndexType index;
  typename ImageType::RegionType::SizeType size;
  size[0] = static_cast<unsigned long>( vcl_floor( 0.5*imageSize[0] ) );
  size[1] = imageSize[1];
  size[2] = 0;
  region.SetSize( size );
  index[0] = static_cast<unsigned long>(
    vcl_floor( imageIndex[0] + 0.25*imageSize[0] ) );
  index[1] = imageIndex[1];

  long lastIndex = imageIndex[2] + imageSize[2]-1;

  for( long n = imageIndex[2]; n <= lastIndex; n++ )
    {
    index[2] = n;
    region.SetIndex( index );

    typedef itk::ExtractImageFilter<ImageType, SliceType> ExtracterType;
    typename ExtracterType::Pointer extracter = ExtracterType::New();
    extracter->SetInput( relabeler->GetOutput() );
    extracter->SetDirectionCollapseToIdentity();
    extracter->SetExtractionRegion( region );
    extracter->Update();


    typedef itk::MinimalPathImageFunction<RealSliceType> FunctionType;
    typename FunctionType::Pointer function = FunctionType::New();

    typename SliceType::IndexType anchor;
    typename SliceType::IndexType free;

    anchor[0] = static_cast<long>(
      vcl_floor( imageIndex[0] + 0.5*imageSize[0] + 0.5 ) );
    anchor[1] = imageIndex[1];

    free[0]   = static_cast<long>(
      vcl_floor( imageIndex[0] + 0.5*imageSize[0] + 0.5 ) );
    free[1]   = imageIndex[1] + imageSize[1]-1;



    typedef itk::BinaryThresholdImageFilter
      <SliceType, RealSliceType> ThresholderType;

    typename ThresholderType::Pointer thresholder
      = ThresholderType::New();

    thresholder->SetInput( extracter->GetOutput() );

    thresholder->SetLowerThreshold( 1 );

    thresholder->SetUpperThreshold( 1 );

    thresholder->SetInsideValue( 0 );

    thresholder->SetOutsideValue( 1 );
    thresholder->Update();

    typename RealSliceType::Pointer speedImage = RealSliceType::New();
    speedImage->SetRegions( thresholder->GetOutput()->GetLargestPossibleRegion() );
    speedImage->SetOrigin( thresholder->GetOutput()->GetOrigin() );
    speedImage->SetSpacing( thresholder->GetOutput()->GetSpacing() );
    speedImage->Allocate();

    float lambda = ( argc > 4 ) ? atof( argv[4] ) : 0.01;

    itk::ImageRegionIteratorWithIndex<RealSliceType> ItD( speedImage,
      speedImage->GetLargestPossibleRegion() );
    for( ItD.GoToBegin(); !ItD.IsAtEnd(); ++ItD )
      {
      typename ImageType::IndexType speedIndex;
      speedIndex[0] = ItD.GetIndex()[0];
      speedIndex[1] = ItD.GetIndex()[1];
      speedIndex[2] = n;
      float t = ItD.GetIndex()[0] - free[0];
      ItD.Set( lambda*t*t + reader->GetOutput()->GetPixel( speedIndex ) );
      }

    function->SetAnchorSeed( anchor );
    function->SetInputImage( speedImage );

    typename FunctionType::OutputType::Pointer path
      = function->EvaluateAtIndex( free );

    typedef itk::PathIterator<SliceType,
      typename FunctionType::OutputType> IteratorType;

    IteratorType It( extracter->GetOutput(), path );
    It.GoToBegin();
    while ( !It.IsAtEnd() )
      {
      typename ImageType::IndexType labelIndex;
      labelIndex[0] = It.GetIndex()[0];
      labelIndex[1] = It.GetIndex()[1];
      labelIndex[2] = n;
      labelImage->SetPixel( labelIndex, n );
      ++It;
      }
    }

//  {
//  typedef itk::ImageFileWriter<ImageType> WriterType;
//  typename WriterType::Pointer writer = WriterType::New();
//  writer->SetInput( labelImage );
//  writer->SetFileName( "slicePoints.nii.gz" );
//  writer->Update();
//  }

  //
  // Smooth the results by creating a B-spline surface spanning
  // all the slices.
  //

  typedef itk::Vector<float, 1> ScalarType;
  typedef itk::PointSet<ScalarType, 2> PointSetType;
  typedef itk::Image<ScalarType, 2> ScalarSliceType;

  typename PointSetType::Pointer points = PointSetType::New();
  points->Initialize();

  unsigned long count = 0;

  typename ImageType::DirectionType direction = labelImage->GetDirection();
  typename ImageType::DirectionType identity;
  identity.SetIdentity();
  labelImage->SetDirection( identity );

  itk::ImageRegionIteratorWithIndex<ImageType> ItL( labelImage,
    labelImage->GetLargestPossibleRegion() );
  for( ItL.GoToBegin(); !ItL.IsAtEnd(); ++ItL )
    {
    if( ItL.Get() != 0 )
      {
      typename PointSetType::PointType point;
      point[0] = ItL.GetIndex()[1];
      point[1] = ItL.GetIndex()[2];
      points->SetPoint( count, point );

      ScalarType scalar;
      scalar[0] = ItL.GetIndex()[0];
      points->SetPointData( count, scalar );

      count++;
      }
    }

  labelImage->SetDirection( direction );

  typedef itk::BSplineScatteredDataPointSetToImageFilter<PointSetType,
    ScalarSliceType> BSplinerType;

  typename BSplinerType::ArrayType ncps;
  ncps.Fill( 4 );

  typename ScalarSliceType::PointType   bsplineOrigin;
  typename ScalarSliceType::SizeType    bsplineSize;
  typename ScalarSliceType::SpacingType bsplineSpacing;

  for( unsigned int d = 1; d < ImageDimension; d++ )
    {
    bsplineSize[d-1]
      = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[d];
    bsplineOrigin[d-1]
      = reader->GetOutput()->GetLargestPossibleRegion().GetIndex()[d];
    bsplineSpacing[d-1] = 1;
    }

  typename BSplinerType::Pointer bspliner = BSplinerType::New();
  bspliner->SetInput( points );
  bspliner->SetGenerateOutputImage( true );
  bspliner->SetSplineOrder( 3 );
  bspliner->SetNumberOfLevels( 7 );
  bspliner->SetNumberOfControlPoints( ncps );
  bspliner->SetSize( bsplineSize );
  bspliner->SetOrigin( bsplineOrigin );
  bspliner->SetSpacing( bsplineSpacing );
  bspliner->Update();

  labelImage->FillBuffer( 0 );

  itk::ImageRegionIteratorWithIndex<ScalarSliceType> ItB( bspliner->GetOutput(),
    bspliner->GetOutput()->GetLargestPossibleRegion() );
  for( ItB.GoToBegin(); !ItB.IsAtEnd(); ++ItB )
    {
    typename ImageType::IndexType index;
    index[0] = static_cast<long>( vcl_floor( ItB.Get()[0] + 0.5 ) );
    index[1] = ItB.GetIndex()[0];
    index[2] = ItB.GetIndex()[1];

    labelImage->SetPixel( index, 1 );
    }

  //
  // Fill in the left side of the image with all '1's.
  //

  itk::ImageSliceIteratorWithIndex<ImageType> ItR( labelImage,
    labelImage->GetLargestPossibleRegion() );
  ItR.SetFirstDirection( 0 );
  ItR.SetSecondDirection( 1 );

  ItR.GoToBegin();
  while( !ItR.IsAtEnd() )
    {
    while( !ItR.IsAtEndOfSlice() )
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
    ItR.NextSlice();
    }

// 		{
// 		typedef itk::ImageFileWriter<ImageType> WriterType;
// 		typename WriterType::Pointer writer = WriterType::New();
// 		writer->SetInput( labelImage );
// 		writer->SetFileName( "bsplineSurface.nii.gz" );
// 		writer->Update();
// 		}

  /**
   * Now combine the two labeled images and write the output.
   */
  {
  itk::ImageRegionIteratorWithIndex<ImageType> ItL( labelImage,
    labelImage->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<ImageType> ItR( relabeler->GetOutput(),
    relabeler->GetOutput()->GetLargestPossibleRegion() );
  for( ItL.GoToBegin(), ItR.GoToBegin(); !ItL.IsAtEnd(); ++ItL, ++ItR  )
    {
    if( ItR.Get() == 1 && ItL.Get() == 1 )
      {
      ItL.Set( 2 );
      }
    else if( ItR.Get() == 1 && ItL.Get() == 0 )
      {
      ItL.Set( 3 );
      }
    else
      {
      ItL.Set( 0 );
      }
    }
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
  if ( argc < 4 )

    {

    std::cout << "Usage: "<< argv[0] << " imageDimension inputImage outputImage" << std::endl;

    exit( 1 );

    }


  switch( atoi( argv[1] ) )
   {
   case 2:
     SplitHeliumLungs2D<2>( argc, argv );
     break;
   case 3:
     SplitHeliumLungs3D<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


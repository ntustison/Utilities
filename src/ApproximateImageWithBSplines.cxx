#include "itkBSplineScatteredDataPointSetToImageFilter.h"
// #include "itkBSplineScatteredDataPointSetToImageFilter2.h"

#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkPointSet.h"
#include "itkVector.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

#include "itkTimeProbe.h"

#include "Common.h"

#include <string>
#include <vector>


template <unsigned int ImageDimension>
int ApproximateImageWithBSplines( int argc, char *argv[] )
{
  typedef float RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<RealImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typename RealImageType::DirectionType direction
    = reader->GetOutput()->GetDirection();

  typename RealImageType::DirectionType identity;
  identity.SetIdentity();
  reader->GetOutput()->SetDirection( identity );

  typedef itk::Vector<RealType, 1> ScalarType;
  typedef itk::Image<ScalarType, ImageDimension> ScalarImageType;

  typedef itk::PointSet<ScalarType, ImageDimension> PointSetType;
  typedef itk::BSplineScatteredDataPointSetToImageFilter
    <PointSetType, ScalarImageType>             BSplineFilterType;
  typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();

  unsigned int splineOrder = atoi( argv[4] );
  typename BSplineFilterType::ArrayType numberOfLevels;
  typename BSplineFilterType::ArrayType ncps;

  std::vector<unsigned int> nlevels
    = ConvertVector<unsigned int>( std::string( argv[5] ) );
  if ( nlevels.size() == 1 )
    {
    numberOfLevels.Fill( nlevels[0] );
    }
  else if ( nlevels.size() == ImageDimension )
    {
    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      numberOfLevels[d] = nlevels[d];
      }
    }
  else
    {
    std::cerr << "Invalid nlevels format." << std::endl;
    return EXIT_FAILURE;
    }

  std::vector<unsigned int> meshSize
    = ConvertVector<unsigned int>( std::string( argv[6] ) );
  if ( meshSize.size() == 1 )
    {
    ncps.Fill( meshSize[0] + splineOrder );
    }
  else if ( meshSize.size() == ImageDimension )
    {
    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      ncps[d] = meshSize[d] + splineOrder;
      }
    }
  else
    {
    std::cerr << "Invalid ncps format." << std::endl;
    return EXIT_FAILURE;
    }

  typename RealImageType::PointType origin = reader->GetOutput()->GetOrigin();
  typename RealImageType::SpacingType spacing
    = reader->GetOutput()->GetSpacing();
  typename RealImageType::SizeType size
    = reader->GetOutput()->GetLargestPossibleRegion().GetSize();


  typename RealImageType::Pointer maskImage;
  if( argc > 7 )
    {
    typename ReaderType::Pointer maskreader = ReaderType::New();
    maskreader->SetFileName( argv[7] );
    try
      {
      maskreader->Update();
      maskImage = maskreader->GetOutput();
      }
    catch(...)
      {
      maskImage = NULL;
      }
    }

  if( argc > 8 )
    {
    std::vector<RealType> org
      = ConvertVector<RealType>( std::string( argv[8] ) );
    if ( org.size() == 1 )
      {
      origin.Fill( org[0] );
      }
    else if ( org.size() == ImageDimension )
      {
      for ( unsigned int d = 0; d < ImageDimension; d++ )
        {
        origin[d] = org[d];
        }
      }
    else
      {
      std::cerr << "Invalid origin format." << std::endl;
      return EXIT_FAILURE;
      }
    }
  if( argc > 9 )
    {
    std::vector<RealType> sp
      = ConvertVector<RealType>( std::string( argv[9] ) );
    if ( sp.size() == 1 )
      {
      spacing.Fill( sp[0] );
      }
    else if ( sp.size() == ImageDimension )
      {
      for ( unsigned int d = 0; d < ImageDimension; d++ )
        {
        spacing[d] = sp[d];
        }
      }
    else
      {
      std::cerr << "Invalid spacing format." << std::endl;
      return EXIT_FAILURE;
      }
    }
  if( argc > 10 )
    {
    std::vector<unsigned int> sz
      = ConvertVector<unsigned int>( std::string( argv[10] ) );
    if ( sz.size() == 1 )
      {
      size.Fill( sz[0] );
      }
    else if ( sz.size() == ImageDimension )
      {
      for ( unsigned int d = 0; d < ImageDimension; d++ )
        {
        size[d] = sz[d];
        }
      }
    else
      {
      std::cerr << "Invalid size format." << std::endl;
      return EXIT_FAILURE;
      }
   }

  typename PointSetType::Pointer imagePoints = PointSetType::New();
  imagePoints->Initialize();

  typename BSplineFilterType::WeightsContainerType::Pointer weights =
    BSplineFilterType::WeightsContainerType::New();
  weights->Initialize();

  itk::ImageRegionConstIteratorWithIndex<RealImageType>
    It( reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );

  unsigned int N = 0;
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if( !maskImage || maskImage->GetPixel( It.GetIndex() ) )
      {
      typename PointSetType::PointType point;
      reader->GetOutput()->TransformIndexToPhysicalPoint( It.GetIndex(), point );

      ScalarType scalar;
      scalar[0] = It.Get();

      imagePoints->SetPointData( N, scalar );
      imagePoints->SetPoint( N, point );
      if( maskImage )
        {
        weights->InsertElement( N, maskImage->GetPixel( It.GetIndex() ) );
        }
      else
        {
        weights->InsertElement( N, 1 );
        }

      N++;
      }
    }

  if( argc > 11 )
    {
    typedef itk::ImageFileReader<RealImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[11] );
    reader->Update();
    reader->GetOutput()->SetDirection( identity );
//
    typedef itk::Vector<RealType, ImageDimension> VectorType;
    typedef itk::Image<VectorType, ImageDimension> VectorImageType;
    typename VectorImageType::Pointer vectorImage = VectorImageType::New();
//
    if( argc > 12 )
      {
      typedef itk::ImageFileReader<VectorImageType> VectorFileReaderType;
      typename VectorFileReaderType::Pointer vreader
        = VectorFileReaderType::New();
      vreader->SetFileName( argv[12] );
      vreader->Update();
      vectorImage = vreader->GetOutput();
      }
    else
      {
      vectorImage->SetOrigin( reader->GetOutput()->GetOrigin() );
      vectorImage->SetSpacing( reader->GetOutput()->GetSpacing() );
      vectorImage->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
      vectorImage->SetDirection( reader->GetOutput()->GetDirection() );
      vectorImage->Allocate();
//
      VectorType V;
      V.Fill( 0 );
      vectorImage->FillBuffer( V );
      }
//
    itk::ImageRegionConstIteratorWithIndex<RealImageType>
      ItI( reader->GetOutput(), reader->GetOutput()->GetRequestedRegion() );
    itk::ImageRegionConstIterator<VectorImageType>
      ItV( vectorImage, vectorImage->GetLargestPossibleRegion() );
    for( ItI.GoToBegin(), ItV.GoToBegin(); !ItI.IsAtEnd(); ++ItI, ++ItV )
      {
      if( !maskImage || maskImage->GetPixel( ItI.GetIndex() ) )
        {
        typename PointSetType::PointType point;
        reader->GetOutput()->TransformIndexToPhysicalPoint( ItI.GetIndex(), point );
  //
        point += ItV.Get();
  //
        bool isOutside = false;
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          if( point[d] < origin[d] || point[d] > origin[d] + (size[d]-1) * spacing[d] )
            {
            isOutside = true;
            break;
            }
          }
        if( isOutside )
          {
          continue;
          }
  //
        ScalarType scalar;
        scalar[0] = ItI.Get();
  //
        imagePoints->SetPointData( N, scalar );
        imagePoints->SetPoint( N, point );
        if( maskImage )
          {
          weights->InsertElement( N, maskImage->GetPixel( It.GetIndex() ) );
          }
        else
          {
          weights->InsertElement( N, 1 );
          }
        N++;
        }
      }
    }

  itk::TimeProbe timer;
  timer.Start();

//  bspliner->DebugOn();
  bspliner->SetOrigin( origin );
  bspliner->SetSpacing( spacing );
  bspliner->SetSize( size );
  bspliner->SetGenerateOutputImage( true );
  bspliner->SetNumberOfLevels( numberOfLevels );
  bspliner->SetSplineOrder( splineOrder );
  bspliner->SetNumberOfControlPoints( ncps );
  bspliner->SetInput( imagePoints );
  bspliner->SetPointWeights( weights );
  bspliner->Update();

  timer.Stop();

  std::cout << "Elapsed Time:  " << timer.GetMean() << std::endl;

  typename RealImageType::Pointer output = RealImageType::New();
  output->SetOrigin( origin );
  output->SetSpacing( spacing );
  output->SetRegions( size );
  output->SetDirection( direction );
  output->Allocate();

  itk::ImageRegionIterator<ScalarImageType> ItB( bspliner->GetOutput(),
    bspliner->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<RealImageType> ItO( output,
    output->GetLargestPossibleRegion() );
  for( ItB.GoToBegin(), ItO.GoToBegin(); !ItB.IsAtEnd(); ++ItB, ++ItO )
    {
    ItO.Set( ItB.Get()[0] );
    }

  typedef itk::ImageFileWriter<RealImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( output );
  writer->Update();

  std::string outputFilename( argv[3] );

  const int index = outputFilename.rfind( ".nii.gz", outputFilename.size() );
  if( index > 1 )
    {
    std::string cpFilename = outputFilename.substr( 0, index ) +
      std::string( "ControlPointLattice.nii.gz" );

    typedef itk::VectorIndexSelectionCastImageFilter<ScalarImageType, RealImageType>
    SelectorType;
    typename SelectorType::Pointer selector = SelectorType::New();
    selector->SetInput( bspliner->GetPhiLattice() );
    selector->SetIndex( 0 );
    selector->Update();

    typename WriterType::Pointer writer2 = WriterType::New();
    writer2->SetFileName( cpFilename.c_str() );
    writer2->SetInput( selector->GetOutput() );
    writer2->Update();
    }

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 7 )
    {
    std::cout << argv[0] << " imageDimension inputImage outputImage order "
              << "numberOfLevels meshSize [mask] [origin] [spacing] [size] "
               << "[inputImage2] [pushInputField2] " << std::endl;
    exit( 0 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     ApproximateImageWithBSplines<2>( argc, argv );
     break;
   case 3:
     ApproximateImageWithBSplines<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


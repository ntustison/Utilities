#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkPointSet.h"
#include "itkTimeProbe.h"
#include "itkVector.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

#include "Common.h"

template <unsigned int ImageDimension>
int SuperResolution( unsigned int argc, char *argv[] )
{

  typedef float                                   RealType;
  typedef itk::Image<RealType, ImageDimension>    ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[3] );

  typename ImageType::Pointer domainImage = reader->GetOutput();
  domainImage->Update();
  domainImage->DisconnectPipeline();

  typedef itk::Vector<RealType, 1>                         ScalarType;
  typedef itk::Image<ScalarType, ImageDimension>           ScalarImageType;
  typedef itk::PointSet<ScalarType, ImageDimension>        PointSetType;
  typedef itk::BSplineScatteredDataPointSetToImageFilter
    <PointSetType, ScalarImageType>                        BSplineFilterType;

  typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();

  typename PointSetType::Pointer bsplinePoints = PointSetType::New();
  bsplinePoints->Initialize();

  typename BSplineFilterType::WeightsContainerType::Pointer weights =
    BSplineFilterType::WeightsContainerType::New();
  weights->Initialize();

  unsigned int splineOrder = 3;
  typename BSplineFilterType::ArrayType numberOfLevels;
  typename BSplineFilterType::ArrayType ncps;

  std::vector<unsigned int> nlevels = ConvertVector<unsigned int>( std::string( argv[5] ) );
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

  std::vector<unsigned int> meshSize = ConvertVector<unsigned int>( std::string( argv[4] ) );
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

  unsigned int N = 0;
  for( unsigned int n = 6; n < argc; n++ )
    {
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[n] );

    typename ImageType::Pointer inputImage = reader->GetOutput();
    inputImage->Update();
    inputImage->DisconnectPipeline();

    itk::ImageRegionConstIteratorWithIndex<ImageType> It( inputImage, inputImage->GetRequestedRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      typename ImageType::PointType imagePoint;
      inputImage->TransformIndexToPhysicalPoint( It.GetIndex(), imagePoint );

      typename ImageType::IndexType index;
      bool isInside = domainImage->TransformPhysicalPointToIndex( imagePoint, index );

      if( !isInside )
        {
        continue;
        }

      ScalarType scalar;
      scalar[0] = It.Get();

      bsplinePoints->SetPointData( N, scalar );
      bsplinePoints->SetPoint( N, imagePoint );
      weights->InsertElement( N, 1 );

      N++;
      }
    }

  itk::TimeProbe timer;
  timer.Start();

  bspliner->SetOrigin( domainImage->GetOrigin() );
  bspliner->SetSpacing( domainImage->GetSpacing() );
  bspliner->SetSize( domainImage->GetRequestedRegion().GetSize() );
  bspliner->SetGenerateOutputImage( true );
  bspliner->SetNumberOfLevels( numberOfLevels );
  bspliner->SetSplineOrder( splineOrder );
  bspliner->SetNumberOfControlPoints( ncps );
  bspliner->SetInput( bsplinePoints );
  bspliner->SetPointWeights( weights );
  bspliner->Update();

  timer.Stop();

  std::cout << "Elapsed Time:  " << timer.GetMean() << std::endl;

  typedef itk::VectorIndexSelectionCastImageFilter<ScalarImageType, ImageType> SelectorType;
  typename SelectorType::Pointer selector = SelectorType::New();
  selector->SetInput( bspliner->GetOutput() );
  selector->SetIndex( 0 );
  selector->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( selector->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 7 )
    {
    std::cout << argv[0] << " imageDimension outputImage domainImage meshSize numberOfLevels inputImage1 ... inputImageN" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     return SuperResolution<2>( argc, argv );
     break;
   case 3:
     return SuperResolution<3>( argc, argv );
     break;
   case 4:
     return SuperResolution<4>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


#include "itkImage.h"
#include "itkFastMarchingImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLabelContourImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"


template <unsigned int ImageDimension>
int FastMarchingImageFilter( unsigned int argc, char *argv[] )
{
  typedef float InternalPixelType;
  typedef itk::Image<InternalPixelType, ImageDimension>  InternalImageType;

  typedef unsigned char OutputPixelType;
  typedef itk::Image<OutputPixelType, ImageDimension> OutputImageType;

  InternalPixelType stoppingValue = atof( argv[5] );

  typedef  itk::ImageFileReader< InternalImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef  itk::FastMarchingImageFilter
    <InternalImageType, InternalImageType> FastMarchingFilterType;
  typename FastMarchingFilterType::Pointer fastMarching
    = FastMarchingFilterType::New();
  fastMarching->SetInput( reader->GetOutput() );

  typedef typename FastMarchingFilterType::NodeContainer           NodeContainer;
  typedef typename FastMarchingFilterType::NodeType                NodeType;
  typedef typename FastMarchingFilterType::LabelImageType LabelImageType;

  typedef itk::ImageFileReader<LabelImageType> LabelImageReaderType;
  typename LabelImageReaderType::Pointer labelImageReader = LabelImageReaderType::New();
  labelImageReader->SetFileName( argv[4] );
  labelImageReader->Update();

  typedef itk::LabelContourImageFilter<LabelImageType, LabelImageType> ContourFilterType;
  typename ContourFilterType::Pointer contour = ContourFilterType::New();
  contour->SetInput( labelImageReader->GetOutput() );
  contour->FullyConnectedOff();
  contour->SetBackgroundValue( itk::NumericTraits<typename LabelImageType::PixelType>::Zero );
  contour->Update();

  typename NodeContainer::Pointer alivePoints = NodeContainer::New();
  alivePoints->Initialize();
  unsigned long aliveCount = 0;
  typename NodeContainer::Pointer trialPoints = NodeContainer::New();
  trialPoints->Initialize();
  unsigned long trialCount = 0;

  itk::ImageRegionIteratorWithIndex<LabelImageType> ItL( labelImageReader->GetOutput(),
    labelImageReader->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<LabelImageType> ItC( contour->GetOutput(),
    contour->GetOutput()->GetLargestPossibleRegion() );
  for( ItL.GoToBegin(), ItC.GoToBegin(); !ItL.IsAtEnd(); ++ItL, ++ItC )
//  for( ItC.GoToBegin(); !ItC.IsAtEnd(); ++ItC )
    {
    if( ItC.Get() != itk::NumericTraits<typename LabelImageType::PixelType>::Zero )
      {
      typename LabelImageType::IndexType position = ItC.GetIndex();

      NodeType node;
      const double value = 0.0;

      node.SetValue( value );
      node.SetIndex( position );
      trialPoints->InsertElement( trialCount++, node );
      }
    else if( ItL.Get() != itk::NumericTraits<typename LabelImageType::PixelType>::Zero )
      {
      typename LabelImageType::IndexType position = ItL.GetIndex();

      NodeType node;
      const double value = 0.0;

      node.SetValue( value );
      node.SetIndex( position );
      alivePoints->InsertElement( aliveCount++, node );
      }
    }
  fastMarching->SetTrialPoints(  trialPoints  );
  fastMarching->SetAlivePoints(  alivePoints  );

  fastMarching->SetStoppingValue( stoppingValue );
  fastMarching->SetTopologyCheck( FastMarchingFilterType::None );
  if( argc > 6 && atoi( argv[6] ) == 1 )
    {
    std::cout << "Strict." << std::endl;
    fastMarching->SetTopologyCheck( FastMarchingFilterType::Strict );
    }
  if( argc > 6 && atoi( argv[6] ) == 2 )
    {
    std::cout << "No handles." << std::endl;
    fastMarching->SetTopologyCheck( FastMarchingFilterType::NoHandles );
    }
  if( argc > 7 )
    {
    fastMarching->SetUseWellComposedness( static_cast<bool>( atoi( argv[7] ) ) );
    }
  if( argc > 8 )
    {
    fastMarching->SetSimplePointConnectivity(
      static_cast<unsigned int>( atoi( argv[8] ) ) );
    }

  try
    {
    fastMarching->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }

//    typedef itk::BinaryThresholdImageFilter
//      <InternalImageType, OutputImageType> ThresholdingFilterType;
//    typename ThresholdingFilterType::Pointer thresholder
//      = ThresholdingFilterType::New();
//
//    thresholder->SetLowerThreshold( 0.00001 );
//    thresholder->SetUpperThreshold( stoppingValue );
//    thresholder->SetOutsideValue( 0 );
//    thresholder->SetInsideValue( 1 );
//    thresholder->SetInput( fastMarching->GetOutput() );
//    thresholder->Update();
//
//   itk::ImageRegionIteratorWithIndex<OutputImageType> ItF( thresholder->GetOutput(),
//     thresholder->GetOutput()->GetLargestPossibleRegion() );
//   for( ItL.GoToBegin(), ItF.GoToBegin(); !ItL.IsAtEnd(); ++ItL, ++ItF )
//     {
//     if( ItL.Get() != itk::NumericTraits<typename LabelImageType::PixelType>::Zero )
//       {
//       ItF.Set( 1 );
//       }
//     }
//
//    typedef  itk::ImageFileWriter<OutputImageType> WriterType;
//    typename WriterType::Pointer writer = WriterType::New();
//    writer->SetInput( thresholder->GetOutput() );
//    writer->SetFileName( argv[3] );
//    writer->Update();


		typedef  itk::ImageFileWriter<InternalImageType> WriterType;
		typename WriterType::Pointer writer = WriterType::New();
		writer->SetInput( fastMarching->GetOutput() );
		writer->SetFileName( argv[3] );
		writer->Update();

  if( argc > 9 )
    {
//     {
//     std::string filename = std::string( argv[9] ) +
//       std::string( "LevelSet.nii.gz" );
//     typedef  itk::ImageFileWriter<InternalImageType> WriterType;
//     typename WriterType::Pointer writer = WriterType::New();
//     writer->SetInput( fastMarching->GetOutput() );
//     writer->SetFileName( filename.c_str() );
//     writer->Update();
//     }

    {
    std::string filename = std::string( argv[9] ) +
      std::string( "LabelMap.nii.gz" );
    typedef itk::ImageFileWriter< LabelImageType > LabelImageWriterType;
    typename LabelImageWriterType::Pointer mapWriter = LabelImageWriterType::New();
    mapWriter->SetInput( fastMarching->GetLabelImage() );
    mapWriter->SetFileName( filename.c_str() );
    mapWriter->Update();
    }

    if( fastMarching->GetConnectedComponentImage() )
						{
						std::string filename = std::string( argv[9] ) +
								std::string( "ConnectedComponents.nii.gz" );
						typedef itk::ImageFileWriter<typename FastMarchingFilterType::ConnectedComponentImageType> LabelImageWriterType;
						typename LabelImageWriterType::Pointer writer = LabelImageWriterType::New();
						writer->SetInput( fastMarching->GetConnectedComponentImage() );
						writer->SetFileName( filename.c_str() );
						writer->Update();
						}
    }

  return 0;
}

int main( int argc, char *argv[] )
{
  if( argc < 6 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension";
    std::cerr << " speedImage outputImage seedImage ";
    std::cerr << " stoppingValue [checkTopology] [useWellComposed]";
    std::cerr << " [simplePointConnectivity] [otherFilePrefix]" << std::endl;
    return 1;
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     FastMarchingImageFilter<2>( argc, argv );
     break;
   case 3:
     FastMarchingImageFilter<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


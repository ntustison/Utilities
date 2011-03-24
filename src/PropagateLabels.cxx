#include "itkBinaryContourImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkFastMarchingImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"

template <unsigned int ImageDimension>
int PropagateLabels( unsigned int argc, char *argv[] )
{
  typedef int LabelType;
  typedef float RealType;

  typedef itk::Image<LabelType, ImageDimension> LabelImageType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

		typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
		typename LabelReaderType::Pointer labelImageReader = LabelReaderType::New();
		labelImageReader->SetFileName( argv[2] );
		labelImageReader->Update();

  bool useEuclidean = false;
  if( argc > 4 && atoi( argv[4] ) )
    {
    useEuclidean = true;
    }

  typename RealImageType::Pointer mask = NULL;

  typedef itk::ImageFileReader<RealImageType> ReaderType;
  typename ReaderType::Pointer imageReader = ReaderType::New();
  if( argc > 5 )
    {
    imageReader->SetFileName( argv[5] );
    imageReader->Update();

    mask = imageReader->GetOutput();
    }

  typename RealImageType::Pointer minimumDistance = RealImageType::New();
  minimumDistance->SetOrigin( labelImageReader->GetOutput()->GetOrigin() );
  minimumDistance->SetSpacing( labelImageReader->GetOutput()->GetSpacing() );
  minimumDistance->SetDirection( labelImageReader->GetOutput()->GetDirection() );
  minimumDistance->SetRegions( labelImageReader->GetOutput()->GetRequestedRegion() );
  minimumDistance->Allocate();
  minimumDistance->FillBuffer( itk::NumericTraits<RealType>::max() );

  typename LabelImageType::Pointer minimumLabels = LabelImageType::New();
  minimumLabels->SetOrigin( labelImageReader->GetOutput()->GetOrigin() );
  minimumLabels->SetSpacing( labelImageReader->GetOutput()->GetSpacing() );
  minimumLabels->SetDirection( labelImageReader->GetOutput()->GetDirection() );
  minimumLabels->SetRegions( labelImageReader->GetOutput()->GetRequestedRegion() );
  minimumLabels->Allocate();
  minimumLabels->FillBuffer( 0 );

  std::vector<LabelType> labels;
  itk::ImageRegionIterator<LabelImageType> It( labelImageReader->GetOutput(),
    labelImageReader->GetOutput()->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if( It.Get() != 0 &&
      std::find( labels.begin(), labels.end(), It.Get() ) == labels.end() )
      {
      labels.push_back( It.Get() );
      }
    }
  std::sort( labels.begin(), labels.end() );

  std::vector<LabelType>::const_iterator it;
  for( it = labels.begin(); it != labels.end(); ++it )
    {
    typedef itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType> ThresholderType;
    typename ThresholderType::Pointer thresholder = ThresholderType::New();
    thresholder->SetInput( labelImageReader->GetOutput() );
    thresholder->SetLowerThreshold( *it );
    thresholder->SetUpperThreshold( *it );
    thresholder->SetInsideValue( 1 );
    thresholder->SetOutsideValue( 0 );
    thresholder->Update();

    typename RealImageType::Pointer distanceImage = RealImageType::New();

    if( useEuclidean )
      {
      typedef itk::SignedMaurerDistanceMapImageFilter
        <LabelImageType, RealImageType> DistancerType;
      typename DistancerType::Pointer distancer = DistancerType::New();
      distancer->SetInput( thresholder->GetOutput() );
      distancer->SetSquaredDistance( true );
      distancer->SetUseImageSpacing( true );
      distancer->SetInsideIsPositive( false );
      distancer->Update();

      distanceImage = distancer->GetOutput();
      }
    else
      {
      typedef itk::BinaryContourImageFilter<LabelImageType, RealImageType>
        ContourFilterType;
      typename ContourFilterType::Pointer contour = ContourFilterType::New();
      contour->SetInput( thresholder->GetOutput() );
      contour->FullyConnectedOff();
      contour->SetBackgroundValue( 0 );
      contour->SetForegroundValue( 1 );
      contour->Update();

      typedef itk::FastMarchingImageFilter<RealImageType, RealImageType>
        FastMarchingFilterType;
      typename FastMarchingFilterType::Pointer fastMarching
        = FastMarchingFilterType::New();

      if( mask )
        {
        fastMarching->SetInput( mask );
        }
      else
        {
        fastMarching->SetSpeedConstant( 1.0 );
        fastMarching->SetOverrideOutputInformation( true );
        fastMarching->SetOutputOrigin( labelImageReader->GetOutput()->GetOrigin() );
        fastMarching->SetOutputSpacing( labelImageReader->GetOutput()->GetSpacing() );
        fastMarching->SetOutputRegion( labelImageReader->GetOutput()->GetRequestedRegion() );
        fastMarching->SetOutputDirection( labelImageReader->GetOutput()->GetDirection() );
        }

      typedef typename FastMarchingFilterType::NodeContainer NodeContainer;
      typedef typename FastMarchingFilterType::NodeType NodeType;
      typename NodeContainer::Pointer trialPoints = NodeContainer::New();
      trialPoints->Initialize();
      typename NodeContainer::Pointer alivePoints = NodeContainer::New();
      alivePoints->Initialize();

      unsigned long trialCount = 0;
      unsigned long aliveCount = 0;

      itk::ImageRegionIteratorWithIndex<RealImageType> ItC(
        contour->GetOutput(), contour->GetOutput()->GetRequestedRegion() );
      itk::ImageRegionIteratorWithIndex<LabelImageType> ItL(
        labelImageReader->GetOutput(), labelImageReader->GetOutput()->GetRequestedRegion() );
      for( ItC.GoToBegin(), ItL.GoToBegin(); !ItC.IsAtEnd(); ++ItC, ++ItL )
        {
        if( ItC.Get() == 1 )
          {
          NodeType node;
          node.SetValue( 0.0 );
          node.SetIndex( ItC.GetIndex() );
          trialPoints->InsertElement( trialCount++, node );
          }
        else if( ItL.Get() == *it )
          {
          NodeType node;
          node.SetValue( 0.0 );
          node.SetIndex( ItC.GetIndex() );
          alivePoints->InsertElement( aliveCount++, node );
          }
        }
      fastMarching->SetTrialPoints( trialPoints );
      fastMarching->SetAlivePoints( alivePoints );
      fastMarching->SetStoppingValue( itk::NumericTraits<RealType>::max() );
//           fastMarching->SetTopologyCheck( FastMarchingFilterType::None );
      fastMarching->Update();

      distanceImage = fastMarching->GetOutput();
      }

    itk::ImageRegionIteratorWithIndex<RealImageType> ItD(
      distanceImage, distanceImage->GetRequestedRegion() );
    itk::ImageRegionIterator<RealImageType> ItM(
      minimumDistance, minimumDistance->GetRequestedRegion() );
    itk::ImageRegionIterator<LabelImageType> ItL(
      minimumLabels, minimumLabels->GetRequestedRegion() );

    ItD.GoToBegin();
    ItM.GoToBegin();
    ItL.GoToBegin();

    while( !ItD.IsAtEnd() )
      {
      if( mask && mask->GetPixel( ItD.GetIndex() ) != 0 )
        {
        if( ItD.Get() < ItM.Get() )
          {
          ItM.Set( ItD.Get() );
          ItL.Set( *it );
          }
        }
      ++ItD;
      ++ItM;
      ++ItL;
      }
    }

  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( minimumLabels );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " imageDimension labelImage outputImage [useEuclidean] "
      << "[speedImage]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     PropagateLabels<2>( argc, argv );
     break;
   case 3:
     PropagateLabels<3>( argc, argv );
     break;
   case 4:
     PropagateLabels<4>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


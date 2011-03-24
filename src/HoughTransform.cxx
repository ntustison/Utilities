#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkHoughTransform2DCirclesImageFilter.h"


template <unsigned int ImageDimension>
int HoughTransform( unsigned int argc, char *argv[] )
{
  typedef float PixelType;
  typedef float RealType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::HoughTransform2DCirclesImageFilter<
    RealType, RealType> HoughFilterType;
  HoughFilterType::Pointer hough = HoughFilterType::New();
  hough->SetInput( reader->GetOutput() );
  hough->SetMinimumRadius( atoi( argv[4] ) );
  hough->SetMaximumRadius( atoi( argv[5] ) );
  hough->SetSigmaGradient(  atof( argv[6] )  );
  hough->SetVariance( atof( argv[7] ) );
  hough->SetNumberOfCircles(  atoi( argv[8] )  );
  hough->SetSweepAngle(  atof( argv[9] )  );
  hough->SetThreshold( atof( argv[10] ) );
  hough->Update();

  typename HoughFilterType::CirclesListType circles =
    hough->GetCircles( hough->GetNumberOfCircles() );
  typename HoughFilterType::CirclesListType::iterator iter;

  unsigned count = 0;
  for ( iter = circles.begin(); iter != circles.end(); ++iter )
    {
    std::cout << count++ << ": ("
      << (*iter)->GetObjectToParentTransform()->GetOffset()[0] << ", "
      << (*iter)->GetObjectToParentTransform()->GetOffset()[1] << ") "
      << (*iter)->GetRadius()[0] << std::endl;
    }

  typename ImageType::Pointer output = ImageType::New();
  output->SetOrigin( reader->GetOutput()->GetOrigin() );
  output->SetSpacing( reader->GetOutput()->GetSpacing() );
  output->SetDirection( reader->GetOutput()->GetDirection() );
  output->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  output->Allocate();
  output->FillBuffer( 0 );

  itk::ImageRegionIteratorWithIndex<ImageType> It( reader->GetOutput(),
    reader->GetOutput()->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    unsigned count = 0;
    PixelType label = 0;
				for ( iter = circles.begin(); iter != circles.end(); ++iter )
						{
      count++;
      RealType centerX = (*iter)->GetObjectToParentTransform()->GetOffset()[0];
      RealType centerY = (*iter)->GetObjectToParentTransform()->GetOffset()[1];
      RealType radius = (*iter)->GetRadius()[0];
      RealType distance = vcl_sqrt( vnl_math_sqr( It.GetIndex()[0] - centerX ) +
        vnl_math_sqr( It.GetIndex()[1] - centerY ) );
      if( distance <= radius )
        {
        label = static_cast<PixelType>( count );
        break;
        }
						}
    output->SetPixel( It.GetIndex(), label );
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( output );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 11 )
    {
    std::cout << argv[0] << " imageDimension inputImage outputImage "
      << "minRadiusInVoxels maxRadiusInVoxels sigmaGradient numberOfCircles variance "
      << "sweepAngle threshold" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     HoughTransform<2>( argc, argv );
     break;
//    case 3:
//      HoughTransform<3>( argc, argv );
//      break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


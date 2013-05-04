#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkArray.h"
#include "itkConstantPadImageFilter.h"

template <unsigned int ImageDimension>
int PadImageForIsotropicBSplineMesh( unsigned int argc, char *argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );

  typename ImageType::Pointer inputImage = reader->GetOutput();
  inputImage->Update();
  inputImage->DisconnectPipeline();

  // the user wants to specify things in terms of spline distance.
  //  1. need to pad the images to get as close to possible to the
  //     requested domain size.
  float splineDistance = 200;
  if( argc > 4 )
    {
    splineDistance = atof( argv[4] );
    }
  unsigned int splineOrder = 3;
  if( argc > 5 )
    {
    splineOrder = static_cast<unsigned int>( atoi( argv[5] ) );
    }
  PixelType padValue = 0;
  if( argc > 6 )
    {
    padValue = static_cast<PixelType>( atof( argv[6] ) );
    }
  bool writeImage = true;
  if( argc > 7 )
    {
    writeImage = static_cast<bool>( atoi( argv[7] ) );
    }

  unsigned long lowerBound[ImageDimension];
  unsigned long upperBound[ImageDimension];

  itk::Array<unsigned int> numberOfControlPoints( ImageDimension );
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    float domain = static_cast<float>(
      inputImage->GetLargestPossibleRegion().GetSize()[d] - 1 ) * inputImage->GetSpacing()[d];
    unsigned int numberOfSpans = static_cast<unsigned int>( vcl_ceil( domain / splineDistance ) );
    unsigned long extraPadding = static_cast<unsigned long>(
      ( numberOfSpans * splineDistance - domain ) / inputImage->GetSpacing()[d] + 0.5 );
    lowerBound[d] = static_cast<unsigned long>( 0.5 * extraPadding );
    upperBound[d] = extraPadding - lowerBound[d];
    numberOfControlPoints[d] = numberOfSpans + splineOrder;
    }

  if( writeImage )
    {
    typedef itk::ConstantPadImageFilter<ImageType, ImageType> PadderType;
    typename PadderType::Pointer padder = PadderType::New();
    padder->SetInput( inputImage );
    padder->SetPadLowerBound( lowerBound );
    padder->SetPadUpperBound( upperBound );
    padder->SetConstant( padValue );
    padder->Update();

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( argv[3] );
    writer->SetInput( padder->GetOutput() );
    writer->Update();
    }

  for( unsigned int d = 0; d < ImageDimension-1; d++ )
    {
    std::cout << numberOfControlPoints[d] - splineOrder << "x";
    }
  std::cout << numberOfControlPoints[ImageDimension-1] - splineOrder << std::endl;

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " imageDimension inputImage outputImage [splineDistance=200] [splineOrder=3] [padValue=0] [writeImage=1]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     PadImageForIsotropicBSplineMesh<2>( argc, argv );
     break;
   case 3:
     PadImageForIsotropicBSplineMesh<3>( argc, argv );
     break;
   case 4:
     PadImageForIsotropicBSplineMesh<4>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


#include <stdio.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkGradientRecursiveGaussianImageFilter.h"

template <unsigned int ImageDimension>
int Gradient( int argc, char *argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::Vector<PixelType, ImageDimension> CovariantVectorType;
  typedef itk::Image<CovariantVectorType, ImageDimension> CovariantVectorImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  bool sigma = 0.0;
  if( argc > 5 )
    {
    sigma = atof( argv[5] );
    }

  bool useImageDirection = true;
  if( argc > 6 )
    {
    useImageDirection = static_cast<bool>( atoi( argv[6] ) );
    }

  typedef itk::GradientRecursiveGaussianImageFilter<ImageType,
    CovariantVectorImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->SetSigma( sigma );
  filter->SetUseImageDirection( useImageDirection );
  filter->Update();

  if( argc > 4 )
    {
    typename ReaderType::Pointer maskReader = ReaderType::New();
    maskReader->SetFileName( argv[4] );

    typename FilterType::OutputPixelType zeroVector;
    zeroVector.Fill( 0.0 );

    itk::ImageRegionIterator<CovariantVectorImageType>
      ItG( filter->GetOutput(),
      filter->GetOutput()->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<ImageType> ItM( maskReader->GetOutput(),
      maskReader->GetOutput()->GetLargestPossibleRegion() );

    for( ItG.GoToBegin(), ItM.GoToBegin(); !ItG.IsAtEnd(); ++ItM, ++ItG )
      {
      if( ItM.Get() == 0 )
        {
        ItG.Set( zeroVector );
        }
      }
    }

  typedef itk::ImageFileWriter<CovariantVectorImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " imageDimension inputImage outputImage [maskImage] [sigma=0.0] [useImageDirection=1] " << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     Gradient<2>( argc, argv );
     break;
   case 3:
     Gradient<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

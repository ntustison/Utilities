#include <stdio.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkGradientImageFilter.h"

template <unsigned int ImageDimension>
int Gradient( int argc, char *argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  bool useImageDirection = true;
  if( argc > 4 )
    {
    useImageDirection = static_cast<bool>( atoi( argv[4] ) );
    }

  typedef itk::GradientImageFilter<ImageType, PixelType, PixelType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->SetUseImageSpacing( true );
  filter->SetUseImageDirection( useImageDirection );
  filter->Update();

  if( argc > 5 )
    {
    typename ReaderType::Pointer maskReader = ReaderType::New();
    maskReader->SetFileName( argv[5] );

    typename FilterType::OutputPixelType zeroVector;
    zeroVector.Fill( 0.0 );

    itk::ImageRegionIterator<typename FilterType::OutputImageType>
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



  typedef itk::ImageFileWriter<typename FilterType::OutputImageType> WriterType;
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
    std::cout << argv[0] << " imageDimension inputImage outputImage [useImageDirection] [maskImage]" << std::endl;
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

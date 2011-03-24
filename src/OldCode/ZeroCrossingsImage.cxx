#include <stdio.h>

#include "itkDiscreteGaussianImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLaplacianImageFilter.H"
#include "itkZeroCrossingImageFilter.h"

template <unsigned int ImageDimension>
int main( int argc, char *argv[] )
{

  typedef float PixelType;
  
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::Image<int, ImageDimension> LabelImageType;

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  float variance = 0.0;
  if( argc > 5 )
    {
    variance = atof( argv[5] );
    } 

  typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> GaussianerType;
  typename GaussianerType::Pointer gaussianer = GaussianerType::New();
  gaussianer->SetInput( reader->GetOutput() );
  gaussianer->SetUseImageSpacing( true );
  gaussianer->SetVariance( variance );
  gaussianer->Update();

  typename ImageType::Pointer image = ImageType::New();
  
  if( atoi( argv[4] ) == 0 )
    {
    typedef itk::GradientMagnitudeImageFilter<ImageType, ImageType> GradientType;
    typename GradientType::Pointer gradient = GradientType::New();
    gradient->SetInput( gaussian->GetOutput() );
    gradient->Update(); 
    image = gradient->GetOutput();
    }   
  else
    {
    typedef itk::LaplacianImageFilter<ImageType, ImageType> LaplacianType;
    typename LaplacianType::Pointer laplacian = LaplacianType::New();
    gradient->SetInput( gaussian->GetOutput() );
    gradient->Update(); 
    image = gradient->GetOutput();
    }



  typedef itk::ZeroCrossingImageFilter<ImageType, LabelImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  filter->SetInput( reader->GetOutput() );
  filter->SetBackgroundValue( 0 );
  filter->SetForegroundValue( 1 );

  FilterType::ArrayType array;

  if ( argc > 3 )
    {
    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      array[d] = atof( argv[3+d] );
      }
    filter->SetVariance( array );
    }
  filter->Update();
  
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( argv[2] );                                          
  writer->Update();
    
  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] 
              << " imageDimension inputImage outputImage type [variance]" << std::endl;
    std::cerr << "   Type:  0. gradient magnitude" << std::endl;
    std::cerr << "          1. Laplacian" << std::endl;          
    exit( 1 );
    }

  switch( atoi( argv[1] ) ) 
   {
   case 2:
     ZeroCrossingsImage<2>( argc, argv );
     break;
   case 3:
     ZeroCrossingsImage<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


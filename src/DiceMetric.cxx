#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkDiceMetricImageFilter.h"

template <unsigned int ImageDimension>
int DiceMetric( int argc, char * argv[] )
{
  typedef float PixelType;
  
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType>  ReaderType;
  typename ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[2] );
  typename ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[3] );

  PixelType label = itk::NumericTraits<PixelType>::One;
  if ( argc > 4 )
    {
    label = static_cast<PixelType>( atof( argv[4] ) );
    } 

  typedef itk::DiceMetricImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput1( reader1->GetOutput() );
  filter->SetInput2( reader2->GetOutput() );
  filter->SetPixelLabel( label );
  filter->Update();
  
  std::cout << label << ": " << filter->GetDiceMetric() << std::endl;

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " ImageDimension image1 "
      << "image2 [label]" << std::endl;
    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) ) 
   {
   case 2:
     DiceMetric<2>( argc, argv );
     break;
   case 3:
     DiceMetric<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


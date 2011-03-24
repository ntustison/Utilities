#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkDiceMetricImageFilter.h"

#include "global.h"

int main( int argc, char * argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " image1 image2 [label]" << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType>  ReaderType;
  ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[1] );
  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[2] );

  PixelType label = itk::NumericTraits<PixelType>::One;
  if ( argc > 3 )
    {
    label = static_cast<PixelType>( atof( argv[3] ) );
    } 

  typedef itk::DiceMetricImageFilter<ImageType, ImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput1( reader1->GetOutput() );
  filter->SetInput2( reader2->GetOutput() );
  filter->SetPixelLabel( label );
  filter->Update();
  
  std::cout << "Dice metric = " << filter->GetDiceMetric() << std::endl;

  

  return EXIT_SUCCESS;
}



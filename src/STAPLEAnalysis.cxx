#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMultiLabelSTAPLEImageFilter.h"
// #include "itkSTAPLEImageFilter.h"
#include "itkTimeProbe.h"

#include <string>
#include <fstream.h>

template <unsigned int ImageDimension>
int STAPLE( unsigned int argc, char *argv[] )
{

  itk::TimeProbe timer;
  timer.Start();

  typedef float RealType;

  typedef itk::Image<unsigned int, ImageDimension> LabelImageType;

  typedef itk::MultiLabelSTAPLEImageFilter<LabelImageType, LabelImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  for ( unsigned int i = 3; i < argc; i++ )
    {
    typedef itk::ImageFileReader<LabelImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[i] );
    reader->Update();

    filter->SetInput( i-3, reader->GetOutput() );
    }
  filter->Update();

  for ( unsigned int i = 3; i < argc; i++ )
    {
//    std::cout << "Confusion matrix assessment: " << argv[i] << std::endl;
//    std::cout << filter->GetConfusionMatrix( i-3 ) << std::endl << std::endl;

    typename FilterType::ConfusionMatrixType conf =
      filter->GetConfusionMatrix( i-3 );
    for( unsigned int j = 0; j < conf.rows(); j++ )
      {
      for( unsigned int k = 0; k < conf.cols(); k++ )
        {
        std::cout << conf( j, k );
        if( j == conf.rows()-1 && k == conf.cols()-1 )
          {
          std::cout << std::endl;
          }
        else
          {
          std::cout << ",";
          }
        }
      }
    }

  typename FilterType::PriorProbabilitiesType priors =
    filter->GetPriorProbabilities();

//  std::cout << "Prior probabilities" << std::endl;
//  for ( unsigned int i = 0; i < priors.GetSize(); ++i )
//    {
//    std::cout << "\t" << i << ": " << priors[i] << "\n";
//    }

  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( filter->GetOutput() );
  writer->Update();

  timer.Stop();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension outputImage "
     << "segmentationImage1 ... segmentationImageN" << std::endl;
    exit( 0 );
    }

   switch( atoi( argv[1] ) )
    {
    case 2:
      STAPLE<2>( argc, argv );
      break;
    case 3:
      STAPLE<3>( argc, argv );
      break;
    default:
       std::cerr << "Unsupported dimension" << std::endl;
       exit( EXIT_FAILURE );
   }
}



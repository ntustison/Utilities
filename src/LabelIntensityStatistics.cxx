#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkLabelStatisticsImageFilter.h"

#include <vector>
#include <algorithm>
#include <iomanip>

template <unsigned int ImageDimension>
int LabelIntensityStatistics( int argc, char *argv[] )
{
  typedef int LabelType;
  typedef float RealType;

  typedef itk::Image<LabelType, ImageDimension> LabelImageType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<RealImageType> ReaderType;
  typename ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName( argv[2] );
  imageReader->Update();

		typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
		typename LabelReaderType::Pointer labelImageReader = LabelReaderType::New();
		labelImageReader->SetFileName( argv[3] );
		labelImageReader->Update();

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

  typedef itk::LabelStatisticsImageFilter<RealImageType, LabelImageType> HistogramGeneratorType;
  typename HistogramGeneratorType::Pointer stats = HistogramGeneratorType::New();
  stats->SetInput( imageReader->GetOutput() );
  stats->SetLabelInput( labelImageReader->GetOutput() );
  stats->Update();

//   std::cout << "                                       "
//             << "************ Individual Labels *************" << std::endl;
  std::cout << std::setw( 10 ) << "Label"
            << std::setw( 17 ) << "Mean"
            << std::setw( 17 ) << "Sigma"
            << std::setw( 17 ) << "Sum"
            << std::setw( 17 ) << "Min"
            << std::setw( 17 ) << "Max" << std::endl;

  std::vector<LabelType>::iterator it;
  for( it = labels.begin(); it != labels.end(); ++it )
    {
    std::cout << std::setw( 10 ) << *it;
    std::cout << std::setw( 17 ) << stats->GetMean( *it );
    std::cout << std::setw( 17 ) << stats->GetSigma( *it );
    std::cout << std::setw( 17 ) << stats->GetSum( *it );
    std::cout << std::setw( 17 ) << stats->GetMinimum( *it );
    std::cout << std::setw( 17 ) << stats->GetMaximum( *it );
    std::cout << std::endl;
    }

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension inputImage labelImage" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     LabelIntensityStatistics<2>( argc, argv );
     break;
   case 3:
     LabelIntensityStatistics<3>( argc, argv );
     break;
   case 4:
     LabelIntensityStatistics<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

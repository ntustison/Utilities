#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkImagePCAShapeModelEstimator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkMultiplyImageFilter.h"
#include "itkNumericSeriesFileNames.h"

#include <string>

template <unsigned int ImageDimension>
int CreatePCAImageDecomposition( int argc, char* argv[] )
{

  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  std::string format = std::string( argv[2] );

  // Get the filenames from the directory
  itk::NumericSeriesFileNames::Pointer names
    = itk::NumericSeriesFileNames::New();
  names->SetSeriesFormat( format.c_str() );
  names->SetStartIndex( atoi( argv[3] ) );
  names->SetEndIndex( atoi( argv[4] ) );
  names->SetIncrementIndex( 1 );

  typedef itk::ImagePCAShapeModelEstimator<ImageType, ImageType> ImagePCAType;
  typename ImagePCAType::Pointer pca = ImagePCAType::New();
  pca->SetNumberOfTrainingImages( names->GetFileNames().size() );
//  pca->SetNumberOfPrincipalComponentsRequired( ( ac > 6 ) ? atoi( av[6] ) : 5 );

  pca->SetNumberOfPrincipalComponentsRequired( names->GetFileNames().size() - 1 );

  for( unsigned int n = 0; n < names->GetFileNames().size(); n++ )
    {
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( ( names->GetFileNames()[n] ).c_str() );
    reader->Update();

    pca->SetInput( n, reader->GetOutput() );
    }


  try
    {
    pca->Update();
    }
  catch ( itk::ExceptionObject &ex )
    {
    std::cout << ex;
    return EXIT_FAILURE;
    }


  std::string outputFormat = std::string( argv[5] );

  // Get the filenames from the directory
  itk::NumericSeriesFileNames::Pointer outputNames
    = itk::NumericSeriesFileNames::New();
  outputNames->SetSeriesFormat( outputFormat.c_str() );
  outputNames->SetStartIndex( 0 );
  outputNames->SetEndIndex( pca->GetNumberOfOutputs() );
  outputNames->SetIncrementIndex( 1 );

  // write out mean shape
  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( pca->GetOutput() );
  writer->SetFileName( ( outputNames->GetFileNames()[0] ).c_str() );
  writer->Update();

  vnl_vector<double> eigenValues = pca->GetEigenValues();

  unsigned int numberOfOutputs = pca->GetNumberOfOutputs();
  if( argc > 6 && atof( argv[6] ) < 1.0 )
    {
    double runningTotal = 0.0;
    double eigenValueTotal = eigenValues.sum();
    for( unsigned int n = 0; n < eigenValues.size(); n++ )
      {
      runningTotal += eigenValues[n];
      if( runningTotal / eigenValueTotal > atof( argv[6] ) )
        {
        numberOfOutputs = n+2;
        break;
        }
      }
    }
  else if( argc > 6 )
    {
    numberOfOutputs = atoi( argv[6] );
    }

  double runningTotal = 0.0;
  double eigenValueTotal = eigenValues.sum();
  for( unsigned int n = 1; n < numberOfOutputs; n++ )
    {
    runningTotal += eigenValues[n-1];
    std::cout << "Eigen " << n << ": "
              << 100 * runningTotal / eigenValueTotal
              << "% (lambda = " << eigenValues[n-1] << ")" << std::endl;

    typedef itk::MultiplyImageFilter
      <ImageType,ImageType,ImageType> MultiplierType;
    typename MultiplierType::Pointer multiplier = MultiplierType::New();
    multiplier->SetInput1( pca->GetOutput( n ) );
    multiplier->SetConstant2( vcl_sqrt( pca->GetEigenValues()[n-1] ) );
    multiplier->Update();

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( multiplier->GetOutput() );
    writer->SetFileName( ( outputNames->GetFileNames()[n] ).c_str() );
    writer->Update();
    }

  return EXIT_SUCCESS;

}

int main( int argc, char *argv[] )
{
  if ( argc < 6 )
    {
    std::cout << argv[0]
      << " imageDimension inputImageSeriesFormat startIndex endIndex "
      << " outputImageSeriesFormat [numberOfPrincipalComponents or percentage<1]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     CreatePCAImageDecomposition<2>( argc, argv );
     break;
   case 3:
     CreatePCAImageDecomposition<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

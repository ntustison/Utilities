#include "itkImageDuplicator.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkImagePCAShapeModelEstimator.h"
#include "itkImagePCADecompositionCalculator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkAddImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkNumericSeriesFileNames.h"

#include <string>

template <unsigned int ImageDimension>
int CreatePCAImageDecompositionModel( int argc, char* argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImagePCAShapeModelEstimator<ImageType, ImageType> ImagePCAType;
  typename ImagePCAType::Pointer pca = ImagePCAType::New();
  pca->SetNumberOfTrainingImages( argc - 4 );
  pca->SetNumberOfPrincipalComponentsRequired( argc - 5 );

  for( unsigned int n = 4; n < argc; n++ )
    {
    std::cout << "Reading " << argv[n] << std::endl;
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[n] );
    reader->Update();

    pca->SetInput( n-4, reader->GetOutput() );
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

  std::string outputFormat = std::string( argv[2] );

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
  double eigenValueTotal = eigenValues.sum();

  unsigned int numberOfOutputs = pca->GetNumberOfOutputs();
  float percentage = 0.0;
  if( atof( argv[3] ) < 1.0 )
    {
    percentage = atof( argv[3] );

    double runningTotal = 0.0;
    for( unsigned int n = 0; n < eigenValues.size(); n++ )
      {
      runningTotal += eigenValues[n];
      if( runningTotal / eigenValueTotal > percentage )
        {
        numberOfOutputs = n + 1;
        break;
        }
      }
    }
  else
    {
    numberOfOutputs = atoi( argv[3] );

    double runningTotal = 0.0;
    for( unsigned int n = 0; n < numberOfOutputs; n++ )
      {
      runningTotal += eigenValues[n];
      }
    percentage = runningTotal / eigenValueTotal;
    }

  std::cout << "\n===========================================" << std::endl;
  std::cout << "Producing " << numberOfOutputs << " basis vectors (" << percentage * 100 << "%)." << std::endl;
  std::cout << "===========================================\n" << std::endl;

  double runningTotal = 0.0;
  for( unsigned int n = 0; n < numberOfOutputs; n++ )
    {
    runningTotal += eigenValues[n];
    std::cout << "Eigen " << n << ": "
              << 100 * runningTotal / eigenValueTotal
              << "% (lambda = " << eigenValues[n] << ")" << std::endl;


    typedef itk::MultiplyImageFilter
      <ImageType,ImageType,ImageType> MultiplierType;
    typename MultiplierType::Pointer multiplier = MultiplierType::New();
    multiplier->SetInput( pca->GetOutput( n ) );
    multiplier->SetConstant( vcl_sqrt( eigenValues[n] / eigenValues[0] ) );
    multiplier->Update();

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( multiplier->GetOutput() );
    writer->SetFileName( ( outputNames->GetFileNames()[n] ).c_str() );
    writer->Update();
    }

  return EXIT_SUCCESS;
}

template <unsigned int ImageDimension>
int PCADecomposeImage( int argc, char* argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImagePCADecompositionCalculator<ImageType, ImageType> PCAType;
  typename PCAType::Pointer pca = PCAType::New();

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  pca->SetImage( reader->GetOutput() );

  typename ReaderType::Pointer meanImageReader = ReaderType::New();
  meanImageReader->SetFileName( argv[3] );
  meanImageReader->Update();

  pca->SetMeanImage( meanImageReader->GetOutput() );

  itk::NumericSeriesFileNames::Pointer basisNames
    = itk::NumericSeriesFileNames::New();
  basisNames->SetSeriesFormat( argv[4] );
  basisNames->SetStartIndex( 1 );
  basisNames->SetEndIndex( atoi( argv[5] ) );
  basisNames->SetIncrementIndex( 1 );

  typename PCAType::BasisImagePointerVector basisImages;
  for( unsigned int n = 0; n < ( basisNames->GetFileNames() ).size(); n++ )
    {
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( basisNames->GetFileNames()[n] );
    reader->Update();

    basisImages.push_back( reader->GetOutput() );
    }

  pca->SetBasisImages( basisImages );
  pca->Compute();
  typename PCAType::BasisVectorType projection = pca->GetProjection();

  typedef typename itk::ImageDuplicator<ImageType> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage( meanImageReader->GetOutput() );
  duplicator->Update();

  typename ImageType::Pointer reconstructedImage = duplicator->GetOutput();

  if( argc > 6 )
    {
    for( unsigned int n = 0; n < projection.size(); n++ )
      {
      typedef itk::MultiplyImageFilter<ImageType, ImageType, ImageType> MultiplierType;
      typename MultiplierType::Pointer multiplier = MultiplierType::New();
      multiplier->SetInput( basisImages[n] );
      multiplier->SetConstant( projection[n] );
      multiplier->Update();

      typedef typename itk::AddImageFilter<ImageType, ImageType, ImageType> AdderType;
      typename AdderType::Pointer adder = AdderType::New();
      adder->SetInput1( reconstructedImage );
      adder->SetInput2( multiplier->GetOutput() );
      adder->Update();

      typename DuplicatorType::Pointer duplicator2 = DuplicatorType::New();
      duplicator2->SetInputImage( adder->GetOutput() );
      duplicator2->Update();

      reconstructedImage = duplicator2->GetOutput();
      }

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( reconstructedImage );
    writer->SetFileName( argv[6] );
    writer->Update();
    }
  else
    {
    std::cout << projection << std::endl;
    }

  return EXIT_SUCCESS;
}


int main( int argc, char *argv[] )
{
  if ( argc < 6 )
    {
    std::cout << "Usage 1: " << argv[0]
      << " imageDimension outputImageSeriesFormat"
      << " numberOfPrincipalComponents_or_percentage<1  inputImages" << std::endl;
    std::cout << "Usage 2: " << argv[0]
      << " imageDimension inputImage meanImage basisImageSeriesFormat numberOfBases [outputImage]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     if( argc == 6 || argc == 7 )
       {
       PCADecomposeImage<2>( argc, argv );
       }
     else
       {
       CreatePCAImageDecompositionModel<2>( argc, argv );
       }
     break;
   case 3:
     if( argc == 6 || argc == 7 )
       {
       PCADecomposeImage<3>( argc, argv );
       }
     else
       {
       CreatePCAImageDecompositionModel<3>( argc, argv );
       }
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkExpImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkSquareImageFilter.h"
#include "itkSubtractImageFilter.h"

#include <string>
#include <vector>

template<class TValue>
TValue Convert( std::string optionString )
{
  TValue value;
  std::istringstream iss( optionString );
  iss >> value;
  return value;
}

template<class TValue>
std::vector<TValue> ConvertVector( std::string optionString )
{
  std::vector<TValue> values;
  std::string::size_type crosspos = optionString.find( 'x', 0 );

  if ( crosspos == std::string::npos )
    {
    values.push_back( Convert<TValue>( optionString ) );
    }
  else
    {
    std::string element = optionString.substr( 0, crosspos ) ;
    TValue value;
    std::istringstream iss( element );
    iss >> value;
    values.push_back( value );
    while ( crosspos != std::string::npos )
      {
      std::string::size_type crossposfrom = crosspos;
      crosspos = optionString.find( 'x', crossposfrom + 1 );
      if ( crosspos == std::string::npos )
        {
        element = optionString.substr( crossposfrom + 1, optionString.length() );
        }
      else
        {
        element = optionString.substr( crossposfrom + 1, crosspos ) ;
        }
      std::istringstream iss2( element );
      iss2 >> value;
      values.push_back( value );
      }
    }
  return values;
}

template <unsigned int ImageDimension>
int CreateGMMProbabilityImages( unsigned int argc, char *argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );

  typename ImageType::Pointer input = reader->GetOutput();
  input->Update();
  input->DisconnectPipeline();

  std::vector<float> means = ConvertVector<float>( std::string( argv[4] ) );
  std::vector<float> sigmas = ConvertVector<float>( std::string( argv[5] ) );

  if( means.size() != sigmas.size() )
    {
    std::cerr << "Error:  means vector size does not match the sigmas vector size." << std::endl;
    return EXIT_FAILURE;
    }

  for( unsigned int n = 0; n < means.size(); n++ )
    {
    float prefactor = 1.0 / ( vcl_sqrt( 2.0 * vnl_math::pi ) * sigmas[n] );

    typedef itk::SubtractImageFilter<ImageType, ImageType, ImageType> SubtracterType;
    typename SubtracterType::Pointer subtracter = SubtracterType::New();
    subtracter->SetInput1( input );
    subtracter->SetConstant2( means[n] );

    typedef itk::SquareImageFilter<ImageType, ImageType> SquarerType;
    typename SquarerType::Pointer squarer = SquarerType::New();
    squarer->SetInput( subtracter->GetOutput() );

    typedef itk::MultiplyImageFilter<ImageType, ImageType, ImageType> MultiplierType;
    typename MultiplierType::Pointer multiplier = MultiplierType::New();
    multiplier->SetInput1( squarer->GetOutput() );
    multiplier->SetConstant2( -1.0 / ( 2.0 * vnl_math_sqr( sigmas[n] ) ) );

    typedef itk::ExpImageFilter<ImageType, ImageType> ExpType;
    typename ExpType::Pointer exp = ExpType::New();
    exp->SetInput( multiplier->GetOutput() );

    typename MultiplierType::Pointer multiplier2 = MultiplierType::New();
    multiplier2->SetInput1( exp->GetOutput() );
    multiplier2->SetConstant2( prefactor );

    std::stringstream ss;
    ss << argv[3] << n+1 << ".nii.gz";

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( ( ss.str() ).c_str() );
    writer->SetInput( multiplier2->GetOutput() );
    writer->Update();
    }

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " imageDimension inputImage outputImagePrefix meanVector stdVector" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     CreateGMMProbabilityImages<2>( argc, argv );
     break;
   case 3:
     CreateGMMProbabilityImages<3>( argc, argv );
     break;
   case 4:
     CreateGMMProbabilityImages<4>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


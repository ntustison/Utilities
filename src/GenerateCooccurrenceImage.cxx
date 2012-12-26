#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkRescaleIntensityImageFilter.h"
#include "itkTextureFeaturesImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

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
int GenerateCooccurrenceImage( int argc, char *argv[] )
{
  typedef float PixelType;
  typedef float RealType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName( argv[2] );
  imageReader->Update();

  // we have to rescale the input image since the cooccurrence filter
  // adds 1 to the upper bound of the joint histogram and small ranges
  // will be greatly affected by this.

  typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescalerType;
  typename RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetInput( imageReader->GetOutput() );
  rescaler->SetOutputMinimum( 0 );
  rescaler->SetOutputMaximum( 10000 );
  rescaler->Update();

  typedef itk::Statistics::TextureFeaturesImageFilter<RealImageType> TextureFilterType;
  typename TextureFilterType::Pointer textureFilter = TextureFilterType::New();
  textureFilter->SetInput( rescaler->GetOutput() );

  if( argc > 5 )
    {
    try
      {
      typedef typename TextureFilterType::MaskImageType MaskImageType;

      typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
      typename MaskReaderType::Pointer labelImageReader = MaskReaderType::New();
      labelImageReader->SetFileName( argv[5] );
      labelImageReader->Update();

      textureFilter->SetMaskImage( labelImageReader->GetOutput() );
  //     textureFilter->SetInsidePixelValue( label );
      }
    catch( ... ) {}
    }


  typename TextureFilterType::RadiusType radius;
  radius.Fill( 10 );
  if( argc > 6 )
    {
    std::vector<unsigned int> rad = ConvertVector<unsigned int>( std::string( argv[6] ) );
    if( rad.size() != ImageDimension )
      {
      radius.Fill( rad[0] );
      }
    else
      {
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        radius[d] = rad[d];
        }
      }
    radius.Fill( atoi( argv[6] ) );
    }
  textureFilter->SetNeighborhoodRadius( radius );

  unsigned int numberOfBins = 64;
  if ( argc > 7 )
    {
    numberOfBins = static_cast<PixelType>( atoi( argv[7] ) );
    }
  textureFilter->SetNumberOfBinsPerAxis( numberOfBins );

  itk::ImageRegionIteratorWithIndex<ImageType> ItI( rescaler->GetOutput(),
    rescaler->GetOutput()->GetLargestPossibleRegion() );

  PixelType maxValue = itk::NumericTraits<PixelType>::NonpositiveMin();
  PixelType minValue = itk::NumericTraits<PixelType>::max();

  typedef typename TextureFilterType::MaskImageType MaskImageType;
  const MaskImageType * mask = textureFilter->GetMaskImage();

  for( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
    {
    if ( !mask || ( mask->GetPixel( ItI.GetIndex() ) == textureFilter->GetInsidePixelValue() ) )
      {
      if ( ItI.Get() < minValue )
        {
        minValue = ItI.Get();
        }
      else if ( ItI.Get() > maxValue )
        {
        maxValue = ItI.Get();
        }
      }
    }
  textureFilter->SetPixelValueMinMax( minValue, maxValue );

  textureFilter->Update();

  typedef typename TextureFilterType::OutputImageType VectorImageType;

  typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType, ImageType> IndexSelectionType;
  typename IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
  indexSelectionFilter->SetInput( textureFilter->GetOutput() );

  int operation = atoi( argv[3] );

  switch( operation )
    {
    case 0:
      {
      indexSelectionFilter->SetIndex( 0 );
      break;
      }
    case 1:
      {
      indexSelectionFilter->SetIndex( 1 );
      break;
      }
    case 2:
      {
      indexSelectionFilter->SetIndex( 2 );
      break;
      }
    case 3:
      {
      indexSelectionFilter->SetIndex( 3 );
      break;
      }
    case 4:
      {
      indexSelectionFilter->SetIndex( 4 );
      break;
      }
    case 5:
      {
      indexSelectionFilter->SetIndex( 5 );
      break;
      }
    case 6:
      {
      indexSelectionFilter->SetIndex( 6 );
      break;
      }
    case 7:
      {
      indexSelectionFilter->SetIndex( 7 );
      break;
      }
    default:
      {
      std::cerr << "Unrecognized option: " << operation << std::endl;
      return EXIT_FAILURE;
      }
    }

  indexSelectionFilter->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[4] );
  writer->SetInput( indexSelectionFilter->GetOutput() );
  writer->Update();

  return 0;
}


int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension inputImage operation "
     << "outputImage [maskImage] [neighborhoodRadius=10] [numberOfBinsPerAxis=64]" << std::endl;
    std::cerr << "operation:  " << std::endl;
    std::cerr << "   0. energy " << std::endl;
    std::cerr << "   1. entropy " << std::endl;
    std::cerr << "   2. correlation " << std::endl;
    std::cerr << "   3. inverse difference moment " << std::endl;
    std::cerr << "   4. inertia " << std::endl;
    std::cerr << "   5. cluster shade " << std::endl;
    std::cerr << "   6. cluster prominence " << std::endl;
    std::cerr << "   7. haralick's correlation " << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     GenerateCooccurrenceImage<2>( argc, argv );
     break;
   case 3:
     GenerateCooccurrenceImage<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


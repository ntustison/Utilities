#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkEuler2DTransform.h"
#include "itkEuler3DTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

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
							std::istringstream iss( element );
							iss >> value;
							values.push_back( value );
							}
					}
			return values;
			}

template <unsigned int ImageDimension>
int RigidTransform3D( int argc, char *argv[] )
{

  typedef float PixelType;

  typedef double RealType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  if( ImageDimension == 3 )
    {
    typedef itk::Euler3DTransform<RealType> TransformType;
    typename TransformType::Pointer transform = TransformType::New();

    std::vector<RealType> rotation =
      ConvertVector<RealType>( std::string( argv[4] ) );

    if( rotation.size() != ImageDimension )
      {
      std::cerr << "Incorrect rotation angle specification." << std::endl;
      return EXIT_FAILURE;
      }
    transform->SetRotation( rotation[0] * vnl_math::pi / 180.0,
                            rotation[1] * vnl_math::pi / 180.0,
                            rotation[2] * vnl_math::pi / 180.0 );

    if( argc > 5 )
      {
      typename TransformType::OutputVectorType translation;
      std::vector<RealType> t = ConvertVector<RealType>( std::string( argv[5] ) );
      if( t.size() != ImageDimension )
        {
        std::cerr << "Incorrect translation specification." << std::endl;
        return EXIT_FAILURE;
        }
      translation[0] = t[0];
      translation[1] = t[1];
      translation[2] = t[2];
      transform->SetTranslation( translation );
      }
    else
      {
      typename TransformType::OutputVectorType translation;
      translation[0] = 0.0;
      translation[1] = 0.0;
      translation[2] = 0.0;
      transform->SetTranslation( translation );
      }
    if( argc > 6 )
      {
      typename TransformType::InputPointType center;
      std::vector<RealType> t = ConvertVector<RealType>( std::string( argv[6] ) );
      if( t.size() != ImageDimension )
        {
        std::cerr << "Incorrect center specification." << std::endl;
        return EXIT_FAILURE;
        }
      center[0] = t[0];
      center[1] = t[1];
      center[2] = t[2];
      transform->SetCenter( center );
      }
    else
      {
      typename TransformType::InputPointType center;

      typename ImageType::SizeType size = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
      typename ImageType::DirectionType direction = reader->GetOutput()->GetDirection();
      typename ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing();
      typename ImageType::PointType origin = reader->GetOutput()->GetOrigin();

      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        center[d] = origin[d] + 0.5 * spacing[d] * ( size[d] - 1.0 );
        }
      transform->SetCenter( direction * center );
      }


    typedef itk::LinearInterpolateImageFunction<ImageType, RealType> LinearInterpolatorType;
    typename LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();
    interpolator->SetInputImage( reader->GetOutput() );

    typedef itk::ResampleImageFilter<ImageType, ImageType, RealType> ResamplerType;
    typename ResamplerType::Pointer resampler = ResamplerType::New();
    resampler->SetTransform( transform );
    resampler->SetInterpolator( interpolator );
    resampler->SetInput( reader->GetOutput() );
    resampler->SetOutputSpacing( reader->GetOutput()->GetSpacing() );
    resampler->SetOutputOrigin( reader->GetOutput()->GetOrigin() );
    resampler->SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
    resampler->Update();

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( argv[3] );
    writer->SetInput( resampler->GetOutput() );
    writer->Update();
    }


 return 0;
}


template <unsigned int ImageDimension>
int RigidTransform2D( int argc, char *argv[] )
{

  typedef float PixelType;

  typedef double RealType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  if( ImageDimension == 2 )
    {
    typedef itk::Euler2DTransform<RealType> TransformType;
    typename TransformType::Pointer transform = TransformType::New();

    std::vector<RealType> rotation =
      ConvertVector<RealType>( std::string( argv[4] ) );

    if( rotation.size() != 1 )
      {
      std::cerr << "Incorrect rotation angle specification." << std::endl;
      return EXIT_FAILURE;
      }
    transform->SetRotation( rotation[0] * vnl_math::pi / 180.0 );

    if( argc > 5 )
      {
      typename TransformType::OutputVectorType translation;
      std::vector<RealType> t = ConvertVector<RealType>( std::string( argv[5] ) );
      if( t.size() != ImageDimension )
        {
        std::cerr << "Incorrect translation specification." << std::endl;
        return EXIT_FAILURE;
        }
      translation[0] = t[0];
      translation[1] = t[1];
      transform->SetTranslation( translation );
      }
    else
      {
      typename TransformType::OutputVectorType translation;
      translation[0] = 0.0;
      translation[1] = 0.0;
      transform->SetTranslation( translation );
      }
    if( argc > 6 )
      {
      typename TransformType::InputPointType center;
      std::vector<RealType> t = ConvertVector<RealType>( std::string( argv[6] ) );
      if( t.size() != ImageDimension )
        {
        std::cerr << "Incorrect center specification." << std::endl;
        return EXIT_FAILURE;
        }
      center[0] = t[0];
      center[1] = t[1];
      transform->SetCenter( center );
      }
    else
      {
      typename TransformType::InputPointType center;

      typename ImageType::SizeType size = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
      typename ImageType::DirectionType direction = reader->GetOutput()->GetDirection();
      typename ImageType::SpacingType spacing = reader->GetOutput()->GetSpacing();
      typename ImageType::PointType origin = reader->GetOutput()->GetOrigin();

      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        center[d] = origin[d] + 0.5 * spacing[d] * ( size[d] - 1.0 );
        }
      transform->SetCenter( direction * center );
      }


    typedef itk::LinearInterpolateImageFunction<ImageType, RealType> LinearInterpolatorType;
    typename LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();
    interpolator->SetInputImage( reader->GetOutput() );

    typedef itk::ResampleImageFilter<ImageType, ImageType, RealType> ResamplerType;
    typename ResamplerType::Pointer resampler = ResamplerType::New();
    resampler->SetTransform( transform );
    resampler->SetInterpolator( interpolator );
    resampler->SetInput( reader->GetOutput() );
    resampler->SetOutputSpacing( reader->GetOutput()->GetSpacing() );
    resampler->SetOutputOrigin( reader->GetOutput()->GetOrigin() );
    resampler->SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
    resampler->Update();

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( argv[3] );
    writer->SetInput( resampler->GetOutput() );
    writer->Update();
    }


 return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension inputImage outputImage "
      << "rotation[0]xrotation[1]x... [translation[0]xtranslation[1]x...] "
      << "[center[0]xcenter[1]x... ]" << std::endl;
    std::cerr << "If center isn't specified, it is assumed that the center is "
      << "defined by the center of the input image.  Angles are in degrees." << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     RigidTransform2D<2>( argc, argv );
     break;
   case 3:
     RigidTransform3D<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

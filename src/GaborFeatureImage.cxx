#include "itkConvolutionImageFilter.h"
#include "itkGaussianInterpolateImageFunction.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkTimeProbe.h"

#include <string>

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



int GaborFeatureImage3D( int argc, char *argv[] )
{
  unsigned int ImageDimension = 3;

  itk::TimeProbe timer;
  timer.Start();

  typedef float RealType;

  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Image<unsigned int, ImageDimension> LabelImageType;

  typedef itk::ImageFileReader<RealImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  /**
   * Construct the gabor image kernel
   *   The idea is that we construct a gabor image kernel and
   *   downsample to a smaller size for image convolution.
   */
  typedef itk::GaborImageSource<ImageType> GaborSourceType;
  GaborSourceType::Pointer gabor = GaborSourceType::New();
  ImageType::SpacingType spacing;
  ImageType::RegionType::SizeType size;
  ImageType::PointType origin;
  ImageType::DirectionType direction;

  /**
   * The following parameter values were based on some empirical testing.
   * May want to change.
   */

  origin.Fill( 0.0 );
  spacing.Fill( 1.0 );
  size.Fill( 255 );
  direction.SetIdentity();

  GaborSourceType::ArrayType mean;
  GaborSourceType::ArrayType sigma;

  for( unsigned int d = 0; d < 3; d++ )
    {
    mean[d] = origin[d] + 0.5 * spacing[d] *
      static_cast<RealType>( size[d] - 1 );
    }
  sigma[0] = 50.0;
  sigma[1] = 75.0;
  sigma[2] = 75.0;

  gabor->SetSpacing( spacing );
  gabor->SetOrigin( origin );
  gabor->SetSize( size );
  gabor->SetDirection( direction );
  gabor->SetFrequency( 0.001 );
  gabor->SetCalculateImaginaryPart( true );
  gabor->SetSigma( sigma );
  gabor->SetMean( mean );
  gabor->Update();

  /**
   * Construct the Gaussian interpolater for the gabor filter resampling
   */
  typedef itk::GaussianInterpolateImageFunction<ImageType, RealType>
    GaussianInterpolatorType;
  typename GaussianInterpolatorType::Pointer g_interpolator
    = GaussianInterpolatorType::New();
  g_interpolator->SetInputImage( reader->GetOutput() );

  double sigma[ImageDimension];
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    sigma[d] = 0.8;
    }
  double alpha = 1.0;
  g_interpolator->SetParameters( sigma, alpha );

  /**
   * Read in the size of the kernel specified by the user.
   */
  ImageType::RegionType::SizeType kernelSize;
  std::vector<unsigned int> ksize = ConvertVector<unsigned int>(
    std::string( argv[4] ) );
  kernelSize[0] = ksize[0];
  kernelSize[1] = ksize[1];
  kernelSize[2] = ksize[2];

  /**
   * Allocate memory for the maximal response image.
   */
  RealImageType::Pointer maxResponse = RealImageType::New();
  maxResponse->SetOrigin( reader->GetOutput()->GetOrigin() );
  maxResponse->SetSpacing( reader->GetOutput()->GetSpacing() );
  maxResponse->SetRegions( reader->GetOutput()->GetBufferedRegion() );
  maxResponse->SetDirection( reader->GetOutput()->GetDirection() );
  maxResponse->Allocate();
  maxResponse->FillBuffer( 0 );

  /**
   * Rotate the gabor kernel around the z axis for the user number of specified
   * steps.
   */
  for( RealType theta = 0.0; theta < 180.0; theta += atof( argv[3] ) )
    {
    typedef itk::Euler3DTransform<RealType> TransformType;
    TransformType::Pointer transform = TransformType::New();
    TransformType::OutputVectorType translation;
    translation.Fill( 0.0 );
    TransformType::InputPointType center;

    for( unsigned int d = 0; d < 3; d++ )
      {
      center[0] = gabor->GetOutput()->GetOrigin()[d] +
        gabor->GetOutput()->GetSpacing()[d] *
        ( gabor->GetOutput()->GetBufferedRegion().GetSize()[d] - 1 );
      }
    transform->SetRotation( 0.0, 0.0, theta );
    transform->SetTranslation( translation );
    transform->SetCenter( center );

    typedef itk::ResampleImageFilter<ImageType, ImageType, RealType> ResamplerType;
    ResamplerType::Pointer resampler = ResamplerType::New();
    resampler->SetTransform( transform );
    resampler->SetInterpolator( g_interpolator );
    resampler->SetInput( gabor->GetOutput() );
    // The output spacing and origin are irrelevant in the convolution
    // calculation.
    resampler->SetOutputSpacing( gabor->GetOutput()->GetSpacing() );
    resampler->SetOutputOrigin( reader->GetOutput()->GetOrigin() );
    resampler->SetSize( kernelSize );
    resampler->Update();

    /**
     * Convolve the input image with the resampled gabor image kernel.
     */
    typedef itk::ConvolutionImageFilter<ImageType> ConvolutionFilterType;
    typename ConvolutionFilterType::Pointer convoluter
      = ConvolutionFilterType::New();
    convoluter->SetInput( reader->GetOutput() );
    convoluter->SetImageKernelInput( resampler->GetOutput() );
    convoluter->NormalizeOn();
    convoluter->Update();

    /**
     *
     */
    ImageRegionIterator<ImageType> ItM( maximalResponseImage,
      maximalResponseImage->GetBufferedRegion() );
    ImageRegionIterator<ImageType> ItG( convoluter->GetOutput(),
      convoluter->GetOutput()->GetBufferedRegion() );
    for( ItM.GoToBegin(), ItG.GoToBegin(); !ItM.IsAtEnd(); ++ItG, ++ItM )
      {
      ImageType::PixelType responseValue = vnl_math_abs( ItG.Get() );
      if( responseValue > ItM.Get() )
        {
        ItM.Set( responseValue );
        }
      }
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( maximalResponseImage );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 7 )
    {
    std::cerr << "Usage: " << argv[0] << " inputImage"
      << "outputImage thetaStepSize kernelSize" << std::endl;
    exit( 0 );
    }

  return GaborFeatureImage3D( argc, argv );
}


#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkScalarToFractalImageFilter.h"
#include "itkTimeProbe.h"

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


class CommandProgressUpdate2D : public itk::Command
{
public:
  typedef  CommandProgressUpdate2D                      Self;
  typedef  itk::Command                                 Superclass;
  typedef  itk::SmartPointer<CommandProgressUpdate2D>  Pointer;
  itkNewMacro( CommandProgressUpdate2D );
protected:

  CommandProgressUpdate2D() : m_CurrentProgress( 0 ) {};

  typedef itk::Image<float, 2> ImageType;

  typedef itk::ScalarToFractalImageFilter<ImageType, ImageType>
    FilterType;

  unsigned int m_CurrentProgress;

public:

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    itk::ProcessObject *po = dynamic_cast<itk::ProcessObject *>( caller );
    if (! po) return;
//    std::cout << po->GetProgress() << std::endl;
    if( typeid( event ) == typeid ( itk::ProgressEvent )  )
      {
      if( this->m_CurrentProgress < 99 )
        {
        this->m_CurrentProgress++;
        if( this->m_CurrentProgress % 10 == 0 )
          {
          std::cout << this->m_CurrentProgress << std::flush;
          }
        else
          {
          std::cout << "*" << std::flush;
          }
        }
      }
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
    itk::ProcessObject *po = dynamic_cast<itk::ProcessObject *>(
      const_cast<itk::Object *>( object ) );
    if (! po) return;

    if( typeid( event ) == typeid ( itk::ProgressEvent )  )
      {
      if( this->m_CurrentProgress < 99 )
        {
        this->m_CurrentProgress++;
        if( this->m_CurrentProgress % 10 == 0 )
          {
          std::cout << this->m_CurrentProgress << std::flush;
          }
        else
          {
          std::cout << "*" << std::flush;
          }
        }
      }
    }
};

class CommandProgressUpdate3D : public itk::Command
{
public:
  typedef  CommandProgressUpdate3D                      Self;
  typedef  itk::Command                                 Superclass;
  typedef  itk::SmartPointer<CommandProgressUpdate3D>   Pointer;
  itkNewMacro( CommandProgressUpdate3D );
protected:

  CommandProgressUpdate3D() : m_CurrentProgress( 0 ) {};

  typedef itk::Image<float, 3> ImageType;

  typedef itk::ScalarToFractalImageFilter<ImageType, ImageType>
    FilterType;

  unsigned int m_CurrentProgress;

public:

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    itk::ProcessObject *po = dynamic_cast<itk::ProcessObject *>( caller );
    if (! po) return;
//    std::cout << po->GetProgress() << std::endl;
    if( typeid( event ) == typeid ( itk::ProgressEvent )  )
      {
      if( this->m_CurrentProgress < 100 )
        {
        this->m_CurrentProgress++;
        if( this->m_CurrentProgress % 10 == 0 )
          {
          std::cout << this->m_CurrentProgress << std::flush;
          }
        else
          {
          std::cout << "*" << std::flush;
          }
        }
      }
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
    itk::ProcessObject *po = dynamic_cast<itk::ProcessObject *>(
      const_cast<itk::Object *>( object ) );
    if (! po) return;

    if( typeid( event ) == typeid ( itk::ProgressEvent )  )
      {
      if( this->m_CurrentProgress < 99 )
        {
        this->m_CurrentProgress++;
        if( this->m_CurrentProgress % 10 == 0 )
          {
          std::cout << this->m_CurrentProgress << std::flush;
          }
        else
          {
          std::cout << "*" << std::flush;
          }
        }
      }
    }
};

template <unsigned int ImageDimension>
int itkScalarToFractalImageFilterTest( int argc, char *argv[] )
{
  typedef float                                 PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer imageReader = ReaderType::New();
  imageReader->SetFileName( argv[2] );
  imageReader->Update();

  typedef itk::ScalarToFractalImageFilter<ImageType, ImageType>
    FractalFilterType;
  typename FractalFilterType::Pointer fractal = FractalFilterType::New();
  fractal->SetInput( imageReader->GetOutput() );

  if( argc > 4 )
    {
    typename FractalFilterType::RadiusType radius;

    std::vector<unsigned int> vector
      = ConvertVector<unsigned int>( std::string( argv[4] ) );
    if( vector.size() > 0 && vector.size() != ImageDimension )
      {
      radius.Fill( vector[0] );
      }
    else
      {
      for ( unsigned int d = 0; d < ImageDimension; d++ )
        {
        radius[d] = vector[d];
        }
      }
    fractal->SetNeighborhoodRadius( radius );
    }
  if( argc > 5 )
    {
    PixelType maskLabel = 1.0;
    if( argc > 6 )
      {
      maskLabel = static_cast<PixelType>( atof( argv[6] ) );
      }

    typename ReaderType::Pointer labelImageReader = ReaderType::New();
    labelImageReader->SetFileName( argv[5] );
    labelImageReader->Update();

    typedef typename FractalFilterType::MaskImageType MaskImageType;

    typedef itk::BinaryThresholdImageFilter<ImageType, MaskImageType>
      ThresholderType;
    typename ThresholderType::Pointer thresholder = ThresholderType::New();
    thresholder->SetInput( labelImageReader->GetOutput() );
    thresholder->SetInsideValue( 1 );
    thresholder->SetOutsideValue( 0 );
    thresholder->SetLowerThreshold( maskLabel );
    thresholder->SetUpperThreshold( maskLabel );
    thresholder->Update();

    fractal->SetMaskImage( thresholder->GetOutput() );
    }

  if( ImageDimension == 2 )
    {
    typename CommandProgressUpdate2D::Pointer observer
      = CommandProgressUpdate2D::New();
    fractal->AddObserver( itk::ProgressEvent(), observer );
    }
  else
    {
    typename CommandProgressUpdate3D::Pointer observer
      = CommandProgressUpdate3D::New();
    fractal->AddObserver( itk::ProgressEvent(), observer );
    }

  try
    {
    itk::TimeProbe timer;

    timer.Start();
    std::cout << "/" << std::flush;
    fractal->Update();
    std::cout << "/" << std::flush;
    timer.Stop();

    std::cout << "   (elapsed time: " << timer.GetMean()
      << ")" << std::endl;
    }
  catch( ... )
    {
    std::cerr << "Exception caught." << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( fractal->GetOutput() );
  writer->SetFileName( argv[3] );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension "
      << "inputImage outputImage [radius] [labelImage] [label]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     itkScalarToFractalImageFilterTest<2>( argc, argv );
     break;
   case 3:
     itkScalarToFractalImageFilterTest<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


#include "antsCommandLineParser.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkDiReCTImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "ItkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTimeProbe.h"

#include <string>
#include <algorithm>
#include <vector>

//template<class TFilter>
//class CommandIterationUpdate : public itk::Command
//{
//public:
//  typedef CommandIterationUpdate   Self;
//  typedef itk::Command             Superclass;
//  typedef itk::SmartPointer<Self>  Pointer;
//  itkNewMacro( Self );
//protected:
//  CommandIterationUpdate() {};
//public:
//
//  void Execute(itk::Object *caller, const itk::EventObject & event)
//    {
//    Execute( (const itk::Object *) caller, event);
//    }
//
//  void Execute(const itk::Object * object, const itk::EventObject & event)
//    {
//    const TFilter * filter =
//      dynamic_cast< const TFilter * >( object );
//    if( typeid( event ) != typeid( itk::IterationEvent ) )
//      { return; }
//    if( filter->GetElapsedIterations() == 1 )
//      {
//      std::cout << "Current level = " << filter->GetCurrentLevel() + 1
//        << std::endl;
//      }
//    std::cout << "  Iteration " << filter->GetElapsedIterations()
//      << " (of "
//      << filter->GetMaximumNumberOfIterations()[filter->GetCurrentLevel()]
//      << ").  ";
//    std::cout << " Current convergence value = "
//      << filter->GetCurrentConvergenceMeasurement()
//      << " (threshold = " << filter->GetConvergenceThreshold()
//      << ")" << std::endl;
//    }
//};

template <unsigned int ImageDimension>
int DiReCT( itk::ants::CommandLineParser *parser )
{
  typedef double RealType;

  typedef itk::Image<unsigned int, ImageDimension> LabelImageType;
  typename LabelImageType::Pointer segmentationImage = NULL;

  typedef itk::Image<RealType, ImageDimension> ImageType;
  typename ImageType::Pointer greyMatterProbabilityImage = NULL;
  typename ImageType::Pointer whiteMatterProbabilityImage = NULL;

  typedef itk::DiReCTImageFilter<LabelImageType, ImageType> DiReCTFilterType;
  typename DiReCTFilterType::Pointer direct = DiReCTFilterType::New();

  //
  // segmentation image
  //
  typename itk::ants::CommandLineParser::OptionType::Pointer
    segmentationImageOption = parser->GetOption( "segmentation-image" );
  if( segmentationImageOption )
    {
    typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
    typename LabelReaderType::Pointer labelReader = LabelReaderType::New();

    std::string inputFile = segmentationImageOption->GetValue();
    labelReader->SetFileName( inputFile.c_str() );

    segmentationImage = labelReader->GetOutput();
    segmentationImage->Update();
    segmentationImage->DisconnectPipeline();
    }
  else
    {
    std::cerr << "Segmentation image not specified." << std::endl;
    return EXIT_FAILURE;
    }
  direct->SetSegmentationImage( segmentationImage );

  //
  // grey matter probability image
  //
  typename itk::ants::CommandLineParser::OptionType::Pointer
    greyMatterOption = parser->GetOption( "grey-matter-probability-image" );
  if( greyMatterOption )
    {
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer gmReader = ReaderType::New();

    std::string gmFile = greyMatterOption->GetValue();
    gmReader->SetFileName( gmFile.c_str() );

    greyMatterProbabilityImage = gmReader->GetOutput();
    greyMatterProbabilityImage->Update();
    greyMatterProbabilityImage->DisconnectPipeline();
    }
  else
    {
    std::cout << "Grey matter probability image not specified. "
      << "Creating one from the segmentation image." << std::endl;

    typedef itk::BinaryThresholdImageFilter<LabelImageType, ImageType>
      ThresholderType;
    typename ThresholderType::Pointer thresholder = ThresholderType::New();
    thresholder->SetInput( segmentationImage );
    thresholder->SetLowerThreshold( 2 );
    thresholder->SetUpperThreshold( 2 );
    thresholder->SetInsideValue( 1 );
    thresholder->SetOutsideValue( 0 );

    typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> SmootherType;
    typename SmootherType::Pointer smoother = SmootherType::New();
    smoother->SetVariance( 1.0 );
    smoother->SetUseImageSpacingOn();
    smoother->SetMaximumError( 0.01 );
    smoother->SetInput( thresholder->GetOutput() );
    smoother->Update();

    greyMatterProbabilityImage = smoother->GetOutput();
    greyMatterProbabilityImage->DisconnectPipeline();
    }
  direct->SetGreyMatterProbabilityImage( greyMatterProbabilityImage );

  //
  // white matter probability image
  //
  typename itk::ants::CommandLineParser::OptionType::Pointer
    whiteMatterOption = parser->GetOption( "white-matter-probability-image" );
  if( whiteMatterOption )
    {
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer wmReader = ReaderType::New();

    std::string wmFile = whiteMatterOption->GetValue();
    wmReader->SetFileName( wmFile.c_str() );

    whiteMatterProbabilityImage = wmReader->GetOutput();
    whiteMatterProbabilityImage->Update();
    whiteMatterProbabilityImage->DisconnectPipeline();
    }
  else
    {
    std::cout << "Grey matter probability image not specified. "
      << "Creating one from the segmentation image." << std::endl;

    typedef itk::BinaryThresholdImageFilter<LabelImageType, ImageType>
      ThresholderType;
    typename ThresholderType::Pointer thresholder = ThresholderType::New();
    thresholder->SetInput( segmentationImage );
    thresholder->SetLowerThreshold( 3 );
    thresholder->SetUpperThreshold( 3 );
    thresholder->SetInsideValue( 1 );
    thresholder->SetOutsideValue( 0 );

    typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> SmootherType;
    typename SmootherType::Pointer smoother = SmootherType::New();
    smoother->SetVariance( 1.0 );
    smoother->SetUseImageSpacingOn();
    smoother->SetMaximumError( 0.01 );
    smoother->SetInput( thresholder->GetOutput() );
    smoother->Update();

    whiteMatterProbabilityImage = smoother->GetOutput();
    whiteMatterProbabilityImage->DisconnectPipeline();
    }
  direct->SetWhiteMatterProbabilityImage( whiteMatterProbabilityImage );

  //
  // thickness prior estimate
  //
  typename itk::ants::CommandLineParser::OptionType::Pointer
    thicknessPriorOption = parser->GetOption( "thickness-prior-estimate" );
  if( thicknessPriorOption )
    {
    direct->SetThicknessPriorEstimate( parser->Convert<RealType>(
      thicknessPriorOption->GetValue() ) );
    }

  //
  // gradient step
  //
  typename itk::ants::CommandLineParser::OptionType::Pointer
    gradientStepOption = parser->GetOption( "gradient-step" );
  if( gradientStepOption )
    {
    direct->SetGradientStep( parser->Convert<RealType>(
      gradientStepOption->GetValue() ) );
    }

  //
  // smoothing sigma
  //
  typename itk::ants::CommandLineParser::OptionType::Pointer
    smoothingSigmaOption = parser->GetOption( "smoothing-sigma" );
  if( smoothingSigmaOption )
    {
    direct->SetSmoothingSigma( parser->Convert<RealType>(
      smoothingSigmaOption->GetValue() ) );
    }


  try
    {
    //direct->DebugOn();
    direct->Update();
    }
  catch( itk::ExceptionObject &e )
    {
    std::cerr << "Exception caught: " << e << std::endl;
    return EXIT_FAILURE;
    }
  direct->Print( std::cout, 3 );

  /**
   * output
   */
  typename itk::ants::CommandLineParser::OptionType::Pointer outputOption =
    parser->GetOption( "output" );
  if( outputOption )
    {
    typedef  itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( direct->GetOutput() );
    writer->SetFileName( ( outputOption->GetValue() ).c_str() );
    writer->Update();
    }

  return EXIT_SUCCESS;
}

void InitializeCommandLineOptions( itk::ants::CommandLineParser *parser )
{
  typedef itk::ants::CommandLineParser::OptionType OptionType;

  {
  std::string description =
    std::string( "This option forces the image to be treated as a specified-" ) +
    std::string( "dimensional image.  If not specified, N4 tries to " ) +
    std::string( "infer the dimensionality from the input image." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "image-dimensionality" );
  option->SetShortName( 'd' );
  option->SetUsageOption( 0, "2/3" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "A segmentation image must be supplied labeling the " ) +
    std::string( "cerebrospinal fluid (csf), grey matter (gm), and " ) +
    std::string( "white matter (wm) with values of '1', '2', and '3'," ) +
    std::string( "respectively." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "segmentation-image" );
  option->SetShortName( 's' );
  option->SetUsageOption( 0, "imageFilename" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "In addition to the segmentation image, a grey matter " ) +
    std::string( "probability image can be used. If no such image is " ) +
    std::string( "supplied, one is created using the segmentation image " ) +
    std::string( "and the specified sigma." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "grey-matter-probability-image" );
  option->SetShortName( 'g' );
  option->SetUsageOption( 0, "imageFilename" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "In addition to the segmentation image, a white matter " ) +
    std::string( "probability image can be used. If no such image is " ) +
    std::string( "supplied, one is created using the segmentation image " ) +
    std::string( "and the specified sigma." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "white-matter-probability-image" );
  option->SetShortName( 'w' );
  option->SetUsageOption( 0, "imageFilename" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Provides a constraint on the final thickness measurement. " ) +
    std::string( "References in the literature give a normal thickness of " ) +
    std::string( "typically 3 mm with normal range from ~2 mm in the " ) +
    std::string( "calcarine cortex to ~4 mm in the precentral gyrus. " ) +
    std::string( "Default = 6 mm." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "thickness-prior-estimate" );
  option->SetShortName( 't' );
  option->SetUsageOption( 0, "thicknessPriorEstimate" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "Gradient step size for the registration algorithm." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "gradient-step" );
  option->SetShortName( 'r' );
  option->SetUsageOption( 0, "stepSize" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "smoothing-sigma" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "smoothing-sigma" );
  option->SetShortName( 'm' );
  option->SetUsageOption( 0, "sigma" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "The output consists of a thickness map defined in the " ) +
    std::string( "cortical gray matter. " );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "output" );
  option->SetShortName( 'o' );
  option->SetUsageOption( 0, "imageFilename" );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Print the help menu (short version)." );

  OptionType::Pointer option = OptionType::New();
  option->SetShortName( 'h' );
  option->SetDescription( description );
  option->AddValue( std::string( "0" ) );
  parser->AddOption( option );
  }

  {
  std::string description = std::string( "Print the help menu." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "help" );
  option->SetDescription( description );
  option->AddValue( std::string( "0" ) );
  parser->AddOption( option );
  }

}

int main( int argc, char *argv[] )
{
  itk::ants::CommandLineParser::Pointer parser =
    itk::ants::CommandLineParser::New();
  parser->SetCommand( argv[0] );

  std::string commandDescription =
    std::string( "DiReCT is a registration based estimate of cortical " ) +
    std::string( "thickness.  It was published in S. R. Das, B. B. " ) +
    std::string( "Avants, M. Grossman, and J. C. Gee, Registration based " ) +
    std::string( "cortical thickness measurement, Neuroimage 2009, " ) +
    std::string( "45:867--879." );

  parser->SetCommandDescription( commandDescription );
  InitializeCommandLineOptions( parser );

  parser->Parse( argc, argv );

  if( argc < 2 || parser->Convert<bool>(
    parser->GetOption( "help" )->GetValue() ) )
    {
    parser->PrintMenu( std::cout, 5, false );
    exit( EXIT_FAILURE );
    }
  else if( parser->Convert<bool>(
    parser->GetOption( 'h' )->GetValue() ) )
    {
    parser->PrintMenu( std::cout, 5, true );
    exit( EXIT_FAILURE );
    }

  // Get dimensionality
  unsigned int dimension = 3;

  itk::ants::CommandLineParser::OptionType::Pointer dimOption =
    parser->GetOption( "image-dimensionality" );
  if( dimOption && dimOption->GetNumberOfValues() > 0 )
    {
      dimension = parser->Convert<unsigned int>( dimOption->GetValue() );
    }
  else
    {
    // Read in the first intensity image to get the image dimension.
    std::string filename;

    itk::ants::CommandLineParser::OptionType::Pointer imageOption =
      parser->GetOption( "input-image" );
    if( imageOption && imageOption->GetNumberOfValues() > 0 )
      {
      if( imageOption->GetNumberOfParameters( 0 ) > 0 )
        {
        filename = imageOption->GetParameter( 0, 0 );
        }
      else
        {
        filename = imageOption->GetValue( 0 );
        }
      }
    else
      {
      std::cerr << "No input images were specified.  Specify an input image"
        << " with the -i option" << std::endl;
      return( EXIT_FAILURE );
      }
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
        filename.c_str(), itk::ImageIOFactory::ReadMode );
    dimension = imageIO->GetNumberOfDimensions();
    }

  std::cout << std::endl << "Running DiReCT for "
    << dimension << "-dimensional images." << std::endl << std::endl;

  switch( dimension )
   {
   case 2:
     DiReCT<2>( parser );
     break;
   case 3:
     DiReCT<3>( parser );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


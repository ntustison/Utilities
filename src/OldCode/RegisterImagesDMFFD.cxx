/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: RegisterImagesDMFFD.cxx,v $
  Language:  C++
  Date:      $Date: 2009/04/22 19:18:27 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itkArray.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkCommandLineOption.h"
#include "itkCommandLineParser.h"
#include "itkDMFFDRegistrationFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkTimeProbe.h"
#include "itkVectorImageFileReader.h"
#include "itkVectorImageFileWriter.h"
#include "itkWarpImageFilter.h"

#include <string>
#include <fstream.h>
#include <iostream>

// Image similarity metrics
#include "itkAvantsMutualInformationRegistrationFunction.h"
#include "itkMeanSquareRegistrationFunction.h"
#include "itkAvantsPDEDeformableRegistrationFunction.h"
#include "itkProbabilisticRegistrationFunction.h"
#include "itkRobustDemonsRegistrationFunction.h"
#include "itkRobustOpticalFlow.h"
#include "itkSectionMutualInformationRegistrationFunction.h"
#include "itkSyNDemonsRegistrationFunction.h"

template <unsigned int Dimension>
int RegisterImages( itk::CommandLineParser *parser )
  {
  bool verbose = parser->template Convert<bool>(
    parser->GetOption( "verbose" )->GetValue() );

  typedef float RealType;

  typedef itk::Image<RealType, Dimension> ImageType;
  typedef itk::Image<RealType, Dimension> RealImageType;

  // Set up the registration filter
  std::cout << "Set up the registration filter (ImageDimension = " 
    << Dimension << ")." << std::endl;

  typedef itk::DMFFDRegistrationFilter<RealImageType,
    RealImageType, RealImageType> RegistrationFilterType;
  typename RegistrationFilterType::Pointer registrationFilter
    = RegistrationFilterType::New();
  registrationFilter->SetVerbose( verbose );
  typedef typename RegistrationFilterType::DeformationFieldType 
    DeformationFieldType;  

  /**
   * Read in fixed and moving images
   */
  if( verbose )
    {
    std::cout << "Reading images and setting up metric." << std::endl;
    }

  typename itk::CommandLineParser::OptionType::Pointer metricOption =
    parser->GetOption( "metric" );

  if( metricOption && metricOption->GetNumberOfParameters() < 2 )
    {
    std::cerr << "Incorrect metric option specification." << std::endl; 
    std::cerr << "   " << metricOption->GetDescription() << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::ImageFileReader<RealImageType> ReaderType;
  typedef itk::ImageFileReader<RealImageType> ReaderType;

  typename ReaderType::Pointer fixedReader = ReaderType::New();
  fixedReader->SetFileName( ( metricOption->GetParameter( 0 ) ).c_str() );
  fixedReader->Update();

  typename ReaderType::Pointer movingReader = ReaderType::New();
  movingReader->SetFileName( ( metricOption->GetParameter( 1 ) ).c_str() );
  movingReader->Update();
  
  bool histogramEqualization = parser->template Convert<bool>( 
    parser->GetOption( "histogram-equalization" )->GetValue() );
  if( histogramEqualization )
    {
    if( verbose )
      { 
      std::cout << "Histogram match the images." << std::endl;
      }

    typedef itk::HistogramMatchingImageFilter<RealImageType, RealImageType> 
      HEFilterType;
    typename HEFilterType::Pointer intensityEqualizeFilter = HEFilterType::New();

    intensityEqualizeFilter->SetReferenceImage( fixedReader->GetOutput() );
    intensityEqualizeFilter->SetInput( movingReader->GetOutput() );
    intensityEqualizeFilter->SetNumberOfHistogramLevels( 255 );
    intensityEqualizeFilter->SetNumberOfMatchPoints( 12 );
    intensityEqualizeFilter->ThresholdAtMeanIntensityOn();
    intensityEqualizeFilter->Update();

    registrationFilter->SetInput( 0, fixedReader->GetOutput() );
    registrationFilter->SetInput( 1, intensityEqualizeFilter->GetOutput() );
    }
  else
    {
    registrationFilter->SetInput( 0, fixedReader->GetOutput() );
    registrationFilter->SetInput( 1, movingReader->GetOutput() );
    }  

  std::string whichMetric = metricOption->GetValue();  

  unsigned int radius = 0;
  if( metricOption->GetNumberOfParameters() >= 3 )
    {
    radius = parser->template Convert<unsigned int>( 
      metricOption->GetParameter( 2 ) ); 
    }

  /**
   * Set up the image metric
   */
  typedef itk::AvantsPDEDeformableRegistrationFunction<RealImageType,
    RealImageType, DeformationFieldType> MetricType;
  typename MetricType::Pointer metric[2];
  typename MetricType::RadiusType metricRadius;
  metricRadius.Fill( radius );

  unsigned int numberOfMISamples = 7500;
  unsigned int numberOfHistogramBins = 64;
  if ( fixedReader->GetOutput()->GetLargestPossibleRegion().GetSize()[0] 
    < 80 || Dimension == 2 )
    {
    numberOfHistogramBins = 32;
    }
  else if ( fixedReader->GetOutput()->GetLargestPossibleRegion().GetSize()[0] 
    > 256 )
    {
    numberOfHistogramBins = 128;
    }

  if( strcmp( "MSQ", whichMetric.c_str() ) )
    {
    typedef itk::SyNDemonsRegistrationFunction
      <RealImageType, RealImageType, DeformationFieldType> Metric0Type;
    typename Metric0Type::Pointer metric0 = Metric0Type::New();
    metric0->SetIntensityDifferenceThreshold( 1e-12 );
//    metric0->SetRobustnessParameter( -1.e12 );
    metric0->SetNormalizeGradient( false );
    metric0->SetGradientStep( 1e6 );
    metric[0] = metric0;
    }
  else if( strcmp( "PR", whichMetric.c_str() ) 
    || strcmp( "CC", whichMetric.c_str() ) )  
    {
    typedef itk::ProbabilisticRegistrationFunction
      <RealImageType, RealImageType, DeformationFieldType> Metric1Type;
    typename Metric1Type::Pointer metric1 = Metric1Type::New();
    metric1->SetFullyRobust( true );
    metric1->SetNormalizeGradient( false );
    metric1->SetGradientStep( 1e6 );
    metric[0] = metric1;
    }
  else if( strcmp( "MI", whichMetric.c_str() ) )  
    {
    typedef itk::AvantsMutualInformationRegistrationFunction
      <RealImageType, RealImageType, DeformationFieldType> Metric2Type;
    typename Metric2Type::Pointer metric2 = Metric2Type::New();
    metric2->SetNumberOfSpatialSamples( numberOfMISamples );
    metric2->SetNumberOfHistogramBins( numberOfHistogramBins );
    metric2->SetNormalizeGradient( true );
    metric2->SetGradientStep( 1e6 );
    metric[0] = metric2;
    }

  if( verbose )
    {
    std::cout << "Employing " << whichMetric 
      << " metric (radius = " << radius << ")." << std::endl; 
    }

  metric[0]->SetRadius( metricRadius );
  metric[0]->SetFixedPointSet( NULL );
  metric[0]->SetMovingPointSet( NULL );

  registrationFilter->SetPDEDeformableMetric( metric[0], 0 );  

  /**
   * Set transformation variables
   */
  if( verbose )
    {
    std::cout << "Setting transformation variables." << std::endl;
    }

  typename itk::CommandLineParser::OptionType::Pointer transformationOption =
    parser->GetOption( "transformation" );
  if( !transformationOption ||
    transformationOption->GetNumberOfParameters() < 2 )
    {
    std::cerr << transformationOption->GetDescription() << std::endl;
    return EXIT_FAILURE;
    }
  registrationFilter->SetSplineOrder( parser->template Convert<unsigned int>(
    transformationOption->GetParameter( 0 ) ) );

  std::vector<unsigned int> meshResolution = parser->template
    ConvertVector<unsigned int>( transformationOption->GetParameter( 1 ) );
  typename RegistrationFilterType::ArrayType initialMeshResolution;
  for ( unsigned int i = 0; i < Dimension; i++ )
    {
    initialMeshResolution[i] = meshResolution[i];
    }
  registrationFilter->SetInitialMeshResolution( initialMeshResolution );

  // Read in initial control point lattice.
  if( transformationOption->GetNumberOfParameters() > 2 )
    {
    typedef itk::VectorImageFileReader<
      typename RegistrationFilterType::RealImageType,
      typename RegistrationFilterType::ControlPointLatticeType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( ( transformationOption->GetParameter( 2 ) ).c_str() );
    try
      {
      reader->Update();
      registrationFilter->SetInitialDeformationFieldControlPoints(
        reader->GetOutput() );
      }
    catch( ... )
      {
      }
    }

  if( transformationOption->GetNumberOfParameters() > 3 )
    {
    registrationFilter->SetEnforceDiffeomorphism( true );
    registrationFilter->SetMinimumJacobian( parser->template 
      Convert<RealType>( transformationOption->GetParameter( 3 ) ) );
    }
  if( transformationOption->GetNumberOfParameters() > 4 )
    {
    registrationFilter->SetMaximumJacobian( parser->template 
      Convert<RealType>( transformationOption->GetParameter( 4 ) ) );
    }
  if( transformationOption->GetNumberOfParameters() > 5 )
    {
    registrationFilter->SetInteriorPenaltyParameter( parser->template 
      Convert<RealType>( transformationOption->GetParameter( 5 ) ) );
    }
    

  /**
   * Set optimization variables
   */
  if( verbose )
    {
    std::cout << "Setting optimization variables." << std::endl;
    }

  typename itk::CommandLineParser::OptionType::Pointer optimizationOption =
    parser->GetOption( "optimization" );
  if( !optimizationOption || optimizationOption->GetNumberOfParameters() < 1 )
    {
    std::cerr << "Incorrect optimization option specification." << std::endl; 
    std::cerr << "   " << optimizationOption->GetDescription() << std::endl;
    return EXIT_FAILURE;
    }
  std::vector<unsigned int> numIterations = parser->template
    ConvertVector<unsigned int>( optimizationOption->GetParameter( 0 ) );
  typename RegistrationFilterType::ResizableUIntArrayType numberOfIterations;
  numberOfIterations.SetSize( numIterations.size() );
  for ( unsigned int i = 0; i < numIterations.size(); i++ )
    {
    numberOfIterations[i] = numIterations[i];
    }
  registrationFilter->SetMaximumNumberOfIterations( numberOfIterations );
  if( optimizationOption->GetNumberOfParameters() > 1 )
    {
    std::vector<typename RegistrationFilterType::RealType> gradFactors = 
      parser->template ConvertVector<typename RegistrationFilterType::RealType>( 
      optimizationOption->GetParameter( 1 ) );
    typename RegistrationFilterType::ResizableRealArrayType 
      gradientScalingFactors;
    gradientScalingFactors.SetSize( numIterations.size() );
    if( gradFactors.size() != numIterations.size() )
      {
      gradientScalingFactors.Fill( gradFactors[0] );
      }
    else
      {
      for ( unsigned int i = 0; i < gradFactors.size(); i++ )
        {
        gradientScalingFactors[i] = gradFactors[i];
        }
      }  
    registrationFilter->SetGradientScalingFactor( gradientScalingFactors );
    }
  if( optimizationOption->GetNumberOfParameters() > 2 )
    {
    registrationFilter->SetLineSearchMaximumIterations( parser->template
      Convert<unsigned int>( optimizationOption->GetParameter( 2 ) ) );
    }
  if( optimizationOption->GetNumberOfParameters() > 3 )
    {
    registrationFilter->SetLineSearchMaximumStepSize( parser->template
      Convert<float>( optimizationOption->GetParameter( 3 ) ) );
    }

  typename RegistrationFilterType::ArrayType fixedImageShrinkFactors;
  typename RegistrationFilterType::ArrayType movingImageShrinkFactors;
  for ( unsigned int i = 0; i < Dimension; i++ )
    {
    RealType fixedSize = static_cast<RealType>(
      fixedReader->GetOutput()->GetLargestPossibleRegion().GetSize()[i] );
    RealType movingSize = static_cast<RealType>(
      movingReader->GetOutput()->GetLargestPossibleRegion().GetSize()[i] );
    RealType fixedFactor = vnl_math_min( vcl_pow( 
      2.0, static_cast<int>( numIterations.size()-1 ) ), fixedSize/32.0 );
    RealType movingFactor = vnl_math_min( vcl_pow( 
      2.0, static_cast<int>( numIterations.size()-1 ) ), movingSize/32.0 );

    fixedImageShrinkFactors[i] = vnl_math_max( 1u, 
      static_cast<unsigned int>( fixedFactor ) );
    movingImageShrinkFactors[i] = vnl_math_max( 1u, 
      static_cast<unsigned int>( movingFactor ) );
    }
  registrationFilter->SetFixedImageShrinkFactors( fixedImageShrinkFactors );
  registrationFilter->SetMovingImageShrinkFactors( movingImageShrinkFactors );

  itk::TimeProbe timer;
  timer.Start();
  std::cout << "START" << std::endl;
  registrationFilter->Update();
  timer.Stop();
  std::cout << "Registration filter run time = " 
    << timer.GetMeanTime() << std::endl;

  if( verbose )
    {
    std::cout << "Writing the output." << std::endl;
    }

  typedef itk::BSplineControlPointImageFilter
    <DeformationFieldType, DeformationFieldType> BSplineControlPointsFilterType;
  typename BSplineControlPointsFilterType::Pointer bspliner 
    = BSplineControlPointsFilterType::New();
  typename BSplineControlPointsFilterType::ArrayType close;

  close.Fill( false );
  bspliner->SetSplineOrder( registrationFilter->GetSplineOrder() );
  bspliner->SetCloseDimension( close );
  bspliner->SetInput( registrationFilter->GetTotalDeformationFieldControlPoints() );
  bspliner->SetOrigin( fixedReader->GetOutput()->GetOrigin() );
  bspliner->SetSize( fixedReader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  bspliner->SetSpacing( fixedReader->GetOutput()->GetSpacing() );
  bspliner->Update();

  typedef itk::VectorImageFileWriter<DeformationFieldType, 
    RealImageType> DeformationFieldWriterType;

  {
  std::string fileName = parser->GetOption( "output" )->GetValue()
    + std::string( "Warp.nii.gz" );

  typename DeformationFieldWriterType::Pointer fieldwriter 
    = DeformationFieldWriterType::New();
  fieldwriter->SetFileName( fileName.c_str() );
  fieldwriter->SetInput( bspliner->GetOutput() );
  fieldwriter->Update();
  }


  typedef itk::WarpImageFilter<RealImageType, 
    RealImageType, DeformationFieldType> WarperType;
  typename WarperType::Pointer warper = WarperType::New();

  warper->SetInput( movingReader->GetOutput() );
  warper->SetDeformationField( bspliner->GetOutput() );
  warper->SetOutputSpacing( fixedReader->GetOutput()->GetSpacing() );
  warper->SetOutputOrigin( fixedReader->GetOutput()->GetOrigin() );
  warper->SetEdgePaddingValue( 0.0 );
  warper->Update();

  {
  std::string fileName = parser->GetOption( "output" )->GetValue()
    + std::string( "Warped.nii.gz" );

  typedef itk::ImageFileWriter<RealImageType> ImageWriterType;
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName( fileName.c_str() );
  writer->SetInput( warper->GetOutput() );
  writer->Update();
  }

  {
  itk::OStringStream buf;
  buf << registrationFilter->GetSplineOrder();

  std::string fileName = parser->GetOption( "output" )->GetValue()
    + std::string( "CPLattice" ) + buf.str()
    + std::string( ".nii.gz" );

  typename DeformationFieldWriterType::Pointer cpwriter = DeformationFieldWriterType::New();
  cpwriter->SetFileName( fileName.c_str() );
  cpwriter->SetInput( registrationFilter->GetTotalDeformationFieldControlPoints() );
  cpwriter->Update();
  }

  return EXIT_SUCCESS;  
}

void InitializeCommandLineOptions( itk::CommandLineParser *parser )
{
  typedef itk::CommandLineParser::OptionType OptionType;

  {
  /**
   * --transformation
   */
  std::string description =
    std::string( "[splineOrder,meshResolution,<initialControlPointLattice>," ) +
    std::string( "<minJacobian>,<maxJacobian>,<interiorPenaltyParameter>]" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "transformation" );
  option->SetShortName( 't' );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  /**
   * --metric
   */
  std::string description =
    std::string( "MSQ,CC,MI[fixedImage,movingImage,radius]" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "metric" );
  option->SetShortName( 'm' );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  /**
   * --optimization
   */
  std::string description =
    std::string( "[maximumNumberOfIterationsAtEachLevel," ) +
    std::string( "<gradientScalingFactor(s)>,<lineSearchIterations>," ) +
    std::string( "<lineSearchMaximumStepSize>]" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "optimization" );
  option->SetShortName( 'z' );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  /**
   * --histogram-equalization
   */
  std::string description = 
    std::string( "Histogram equalization of moving image." );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "histogram-equalization" );
  option->SetShortName( 'h' );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "output" );
  option->SetShortName( 'o' );
  option->SetDescription( "File prefix" );
  parser->AddOption( option );
  }

  {
  /**
   * --verbose
   */
  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "verbose" );
  option->SetShortName( 'v' );
  option->SetDescription( "" );
  parser->AddOption( option );
  }

  {
  /**
   * --help
   */
  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "help" );
  option->SetShortName( 'h' );
  option->SetDescription( "Print menu." );
  option->AddValue( "0" );
  parser->AddOption( option );
  }
}

int main( int argc, char *argv[] )
{
  if ( argc < 2 )
    {
    std::cout << "Usage: " << argv[0] << " ImageDimension args" << std::endl;
    exit( 1 );
    }

  itk::CommandLineParser::Pointer parser = itk::CommandLineParser::New();
  parser->SetCommand( argv[0] );
  parser->SetCommandDescription( "DMFFD Image Registration" );
  InitializeCommandLineOptions( parser );

  parser->Parse( argc, argv );

  if( argc < 3 || parser->Convert<bool>( 
    parser->GetOption( "help" )->GetValue() ) )
    {
    parser->PrintMenu( std::cout, 5 );
    exit( EXIT_FAILURE );
    }
    
  switch( atoi( argv[1] ) )
   {
   case 2:
     {
     RegisterImages<2>( parser );
     break;
     }
   case 3:
     {
     RegisterImages<3>( parser );
     break;
     }
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

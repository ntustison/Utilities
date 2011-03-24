/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: RegisterPointSetsDMFFD.cxx,v $
  Language:  C++
  Date:      $Date: 2010/01/09 04:02:50 $
  Version:   $Revision: 1.6 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itkArray.h"
#include "itkCommandLineOption.h"
#include "itkCommandLineParser.h"
#include "itkDMFFDLabeledPointSetRegistrationFilter.h"
#include "itkImageFileReader.h"
#include "itkLabeledPointSetFileReader.h"
#include "itkLabeledPointSetFileWriter.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkTimeProbe.h"
#include "itkTransformFileWriter.h"
#include "itkVectorImageFileReader.h"
#include "itkVectorImageFileWriter.h"

#include <string>
#include <fstream.h>
#include <iostream>

template <unsigned int Dimension>
int itkDMFFDLabeledPointSetRegistrationFilterTest( itk::CommandLineParser *parser )
  {
  bool verbose = parser->template Convert<bool>(
    parser->GetOption( "verbose" )->GetValue() );

  typedef long LabelType;
  typedef float RealType;

  typedef itk::PointSet<LabelType, Dimension> PointSetType;
  typedef itk::DMFFDLabeledPointSetRegistrationFilter<PointSetType>
    RegistrationFilterType;
  typename RegistrationFilterType::Pointer registrationFilter
    = RegistrationFilterType::New();
  registrationFilter->SetVerbose( verbose );

  std::cout << verbose << std::endl;

  /**
   * Read in point-sets
   */
  if( verbose )
    {
    std::cout << "Reading point sets." << std::endl;
    }

  typename itk::CommandLineParser::OptionType::Pointer pointSetOption =
    parser->GetOption( "point-sets" );

  if( pointSetOption && pointSetOption->GetNumberOfParameters() < 2 )
    {
    std::cerr << "Incorrect point set option specification." << std::endl;
    std::cerr << "   " << pointSetOption->GetDescription() << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::LabeledPointSetFileReader<PointSetType> ReaderType;

  typename ReaderType::Pointer fixedPointsReader = ReaderType::New();
  fixedPointsReader->SetFileName(
    ( pointSetOption->GetParameter( 0 ) ).c_str() );
  fixedPointsReader->Update();

  typename ReaderType::Pointer movingPointsReader = ReaderType::New();
  movingPointsReader->SetFileName(
    ( pointSetOption->GetParameter( 1 ) ).c_str() );
  movingPointsReader->Update();

  if( pointSetOption->GetNumberOfParameters() < 4 )
    {
    registrationFilter->SetUseInputAsSamples( true );
    }
  if( pointSetOption->GetNumberOfParameters() > 2 )
    {
    registrationFilter->SetNumberOfFixedSamples( parser->template
      Convert<unsigned long>( pointSetOption->GetParameter( 2 ) ) );
    }
  if( pointSetOption->GetNumberOfParameters() > 3 )
    {
    registrationFilter->SetNumberOfMovingSamples( parser->template
      Convert<unsigned long>( pointSetOption->GetParameter( 3 ) ) );
    }

  /**
   * If the labels option has not been invoked, assign the fixed and moving
   * points appropriately.  Otherwise, revise the fixed and moving
   * point sets according to the label option parameters.
   */

  typename itk::CommandLineParser::OptionType::Pointer labelOption =
    parser->GetOption( "labels" );

  if( !labelOption || labelOption->GetNumberOfParameters() == 0 )
    {
    registrationFilter->SetInput( 0, fixedPointsReader->GetOutput() );
    registrationFilter->SetInput( 1, movingPointsReader->GetOutput() );

    if( verbose )
      {
      std::cout << "Number of fixed points: "
        << fixedPointsReader->GetOutput()->GetNumberOfPoints() << std::endl;
      std::cout << "    Number of fixed labels: "
        << fixedPointsReader->GetNumberOfLabels() << std::endl;
      std::cout << "    Distinct fixed labels: ";
      for( unsigned int n = 0;
        n < fixedPointsReader->GetNumberOfLabels(); n++ )
        {
        std::cout << fixedPointsReader->GetLabelSet()->operator[]( n ) << " ";
        }
      std::cout << std::endl;

      std::cout << "Number of moving points: "
        << movingPointsReader->GetOutput()->GetNumberOfPoints() << std::endl;
      std::cout << "    Number of moving labels: "
        << movingPointsReader->GetNumberOfLabels() << std::endl;
      std::cout << "    Distinct moving labels: ";
      for( unsigned int n = 0;
        n < movingPointsReader->GetNumberOfLabels(); n++ )
        {
        std::cout << movingPointsReader->GetLabelSet()->operator[]( n ) << " ";
        }
      std::cout << std::endl;
      }
    }
  else
    {
    typename PointSetType::Pointer movingPoints = PointSetType::New();
    typename PointSetType::Pointer fixedPoints = PointSetType::New();

    std::vector<LabelType> labels = parser->template ConvertVector<LabelType>(
      labelOption->GetParameter( 0 ) );
    std::vector<float> labelPercentages;
    if( labelOption->GetNumberOfParameters() > 1 )
      {
      labelPercentages = parser->template ConvertVector<float>(
        labelOption->GetParameter( 1 ) );
      }
    if( labelOption->GetNumberOfParameters() > 2 )
      {
      std::vector<float> labelGradientWeights = parser->template
        ConvertVector<float>( labelOption->GetParameter( 2 ) );

      if( labelGradientWeights.size() != labels.size() )
        {
        std::cerr << "Size of weights does not match size of labels."
          << std::endl;
        return EXIT_FAILURE;
        }
      typename RegistrationFilterType::LabelWeightsMapType labelWeights;
      for( unsigned int i = 0; i < labels.size(); i++ )
        {
        labelWeights.insert( std::pair<long, float>(
          static_cast<long>( labels[i] ), labelGradientWeights[i] ) );
        }
      registrationFilter->SetLabelWeights( labelWeights );
      }

    movingPoints->Initialize();
    typename PointSetType::PointsContainerIterator ItM
      = movingPointsReader->GetOutput()->GetPoints()->Begin();
    typename PointSetType::PointDataContainerIterator ItMD
      = movingPointsReader->GetOutput()->GetPointData()->Begin();

    unsigned long count = 0;
    while( ItM != movingPointsReader->GetOutput()->GetPoints()->End() )
      {
      std::vector<LabelType>::iterator loc
        = std::find( labels.begin(), labels.end(), ItMD.Value() );
      if( loc != labels.end() )
        {
        int index = static_cast<int>( std::distance( labels.begin(), loc ) );
        if( labelPercentages.size() > 0 && labelPercentages[index] < 1.0 )
          {
          typedef itk::Statistics::MersenneTwisterRandomVariateGenerator
            GeneratorType;
          typename GeneratorType::Pointer generator = GeneratorType::New();
          generator->SetSeed();

          if( generator->GetVariateWithClosedRange()
            <= labelPercentages[index] )
            {
            movingPoints->SetPoint( count, ItM.Value() );
            movingPoints->SetPointData( count, ItMD.Value() );
            count++;
            }
          }
        else
          {
          movingPoints->SetPoint( count, ItM.Value() );
          movingPoints->SetPointData( count, ItMD.Value() );
          count++;
          }
        }
      ++ItM;
      ++ItMD;
      }


    fixedPoints->Initialize();
    typename PointSetType::PointsContainerIterator ItF
      = fixedPointsReader->GetOutput()->GetPoints()->Begin();
    typename PointSetType::PointDataContainerIterator ItFD
      = fixedPointsReader->GetOutput()->GetPointData()->Begin();

    count = 0;
    while( ItF != fixedPointsReader->GetOutput()->GetPoints()->End() )
      {
      std::vector<LabelType>::iterator loc
        = std::find( labels.begin(), labels.end(), ItFD.Value() );
      if( loc != labels.end() )
        {
        int index = static_cast<int>( std::distance( labels.begin(), loc ) );

        if( labelPercentages.size() > 0 && labelPercentages[index] < 1.0 )
          {
          typedef itk::Statistics::MersenneTwisterRandomVariateGenerator
            GeneratorType;
          typename GeneratorType::Pointer generator = GeneratorType::New();
          generator->SetSeed();

          if( generator->GetVariateWithClosedRange()
            <= labelPercentages[index] )
            {
            fixedPoints->SetPoint( count, ItF.Value() );
            fixedPoints->SetPointData( count, ItFD.Value() );
            count++;
            }
          }
        else
          {
          fixedPoints->SetPoint( count, ItF.Value() );
          fixedPoints->SetPointData( count, ItFD.Value() );
          count++;
          }
        }
      ++ItF;
      ++ItFD;
      }

    std::cout << "Number of fixed points: "
      << fixedPoints->GetNumberOfPoints() << std::endl;
    std::cout << "    Number of fixed labels: "
      << labels.size() << std::endl;
    std::cout << "    Distinct fixed labels: ";
    for ( unsigned int n = 0; n < labels.size(); n++ )
      {
      std::cout << labels[n] << " ";
      }
    std::cout << std::endl;
    std::cout << "Number of moving points: "
      << movingPoints->GetNumberOfPoints() << std::endl;
    std::cout << "    Number of moving labels: "
      << labels.size() << std::endl;
    std::cout << "    Distinct moving labels: ";
    for ( unsigned int n = 0; n < labels.size(); n++ )
      {
      std::cout << labels[n] << " ";
      }
    std::cout << std::endl;

    registrationFilter->SetInput( 0, fixedPoints );
    registrationFilter->SetInput( 1, movingPoints );
    }

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

  // Read in domain image for size, origin, and spacing information.
  if( transformationOption->GetNumberOfParameters() > 3 )
    {
    typedef itk::Image<char, Dimension> ImageType;
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( ( transformationOption->GetParameter( 3 ) ).c_str() );
    try
      {
      reader->Update();

      registrationFilter->SetSize(
        reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
      registrationFilter->SetOrigin( reader->GetOutput()->GetOrigin() );
      registrationFilter->SetSpacing( reader->GetOutput()->GetSpacing() );
      }
    catch( ... )
      {
      }
    }

  std::vector<unsigned int> meshResolution = parser->template
    ConvertVector<unsigned int>( transformationOption->GetParameter( 1 ) );
  typename RegistrationFilterType::ArrayType initialMeshResolution;
  if( meshResolution.size() == 1 &&
    transformationOption->GetNumberOfParameters() > 3 )
    {
    for( unsigned int i = 0; i < Dimension; i++ )
      {
      initialMeshResolution[i] = static_cast<unsigned int>( 0.5 +
        ( registrationFilter->GetSize()[i] - 1 ) *
        registrationFilter->GetSpacing()[i] / meshResolution[0] );
      }
    }
  else
    {
    for ( unsigned int i = 0; i < Dimension; i++ )
      {
      initialMeshResolution[i] = meshResolution[i];
      }
    }
  registrationFilter->SetInitialMeshResolution( initialMeshResolution );

  // Should initial similarity transform be calculated?
  if( transformationOption->GetNumberOfParameters() > 2 )
    {
    registrationFilter->SetCalculateInitialSimilarityTransform(
      parser->template Convert<bool>(
      transformationOption->GetParameter( 2 ) ) );
    }


  // Read in initial control point lattice.
  if( transformationOption->GetNumberOfParameters() > 4 )
    {
    typedef itk::VectorImageFileReader<
      typename RegistrationFilterType::RealImageType,
      typename RegistrationFilterType::ControlPointLatticeType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( ( transformationOption->GetParameter( 4 ) ).c_str() );
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

  // Read in directionality vector.
  if( transformationOption->GetNumberOfParameters() > 5 )
    {
    std::vector<unsigned int> dir = parser->template
      ConvertVector<unsigned int>( transformationOption->GetParameter( 5 ) );
    typename RegistrationFilterType::ArrayType directionality;
    for( unsigned int d = 0; d < Dimension; d++ )
      {
      directionality[d] = dir[d];
      }
    registrationFilter->SetDirectionality( directionality );
    }

  /**
   * Set similarity variables
   */
  if( verbose )
    {
    std::cout << "Setting similarity variables." << std::endl;
    }

  typename itk::CommandLineParser::OptionType::Pointer similarityOption =
    parser->GetOption( "similarity" );
  if( !similarityOption || similarityOption->GetNumberOfParameters() < 3 )
    {
    std::cerr << "Incorrect similarity option specification." << std::endl;
    std::cerr << "   " << similarityOption->GetDescription() << std::endl;
    return EXIT_FAILURE;
    }
  registrationFilter->SetAlpha( parser->template Convert<float>(
    similarityOption->GetParameter( 0 ) ) );
  registrationFilter->SetAnnealingRate( parser->template Convert<float>(
    similarityOption->GetParameter( 1 ) ) );

  std::vector<float> sigmas = parser->template ConvertVector<float>(
    similarityOption->GetParameter( 2 ) );
  typename RegistrationFilterType::ResizableRealArrayType pointSetSigmas;
  pointSetSigmas.SetSize( sigmas.size() );
  for( unsigned int d = 0; d < sigmas.size(); d++ )
    {
    pointSetSigmas[d] = sigmas[d];
    }
  registrationFilter->SetPointSetSigmas( pointSetSigmas );
  unsigned int parameterNumber = 3;
  if( similarityOption->GetNumberOfParameters() > parameterNumber )
    {
    registrationFilter->SetUseRegularizationTerm( parser->template
      Convert<bool>( similarityOption->GetParameter( parameterNumber++ ) ) );
    }
  if( similarityOption->GetNumberOfParameters() > parameterNumber )
    {
    registrationFilter->SetEvaluationKNeighborhood( parser->template
      Convert<unsigned long>(
      similarityOption->GetParameter( parameterNumber++ ) ) );
    }
  if( similarityOption->GetNumberOfParameters() > parameterNumber )
    {
    registrationFilter->SetUseAnisotropicCovariances( parser->template
      Convert<bool>( similarityOption->GetParameter( parameterNumber++ ) ) );
    }
  if( similarityOption->GetNumberOfParameters() > parameterNumber )
    {
    registrationFilter->SetCovarianceKNeighborhood( parser->template
      Convert<unsigned long>(
      similarityOption->GetParameter( parameterNumber++ ) ) );
    }
  if( similarityOption->GetNumberOfParameters() > parameterNumber )
    {
    registrationFilter->SetKernelSigma( parser->template
      Convert<float>( similarityOption->GetParameter( parameterNumber++ ) ) );
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


  /**
   * Run registration
   */
  if( verbose )
    {
    std::cout << "Starting registration" << std::endl;
    }
  itk::TimeProbe timer;
  timer.Start();
  registrationFilter->Update();
  timer.Stop();
  if( verbose )
    {
    std::cout << "Registration filter run time = "
      << timer.GetMeanTime() << std::endl;
    }

  /**
   * Write the outputs
   */
  if( verbose )
    {
    std::cout << "Writing the output." << std::endl;
    }

  // Write warped points to vtk file
  std::string fileName = parser->GetOption( "output" )->GetValue()
    + std::string( "WarpedPoints.vtk" );

  typedef itk::LabeledPointSetFileWriter<PointSetType> PointSetWriterType;
  typename PointSetWriterType::Pointer pswriter = PointSetWriterType::New();
  pswriter->SetFileName( fileName.c_str() );
  pswriter->SetInput( registrationFilter->GetOutput() );
  pswriter->Update();

  // Write similarity transform to file
  {
  std::string fileName = parser->GetOption( "output" )->GetValue()
    + std::string( "Similarity.txt" );

  typedef itk::MatrixOffsetTransformBase<double,
    Dimension, Dimension> TransformType;
  typename TransformType::Pointer transform = TransformType::New();
  typename TransformType::OutputVectorType offset;
  for( unsigned int d = 0; d < Dimension; d++ )
    {
    offset[d] = registrationFilter->GetFixedPointSetCentroid()[d]
      - registrationFilter->GetMovingPointSetCentroid()[d]
      * registrationFilter->GetScaleFactor();
    }
  transform->SetOffset( offset );
  typename TransformType::MatrixType matrix;
  matrix.SetIdentity();
  matrix *= registrationFilter->GetScaleFactor();
  transform->SetMatrix( matrix );

  typedef itk::TransformFileWriter WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( fileName.c_str() );
  writer->SetInput( transform );
  writer->Update();
  }

  // Write deformation field to component nifti images.

  {
  std::string fileName = parser->GetOption( "output" )->GetValue()
    + std::string( "Warp.nii.gz" );
  typedef typename RegistrationFilterType::ControlPointFilterType
    ControlPointFilterType;
  typename ControlPointFilterType::Pointer bspliner
    = ControlPointFilterType::New();
  bspliner->SetSplineOrder( registrationFilter->GetSplineOrder() );
  bspliner->SetInput(
    registrationFilter->GetTotalDeformationFieldControlPoints() );
  bspliner->SetOrigin( registrationFilter->GetOrigin() );
  bspliner->SetSize( registrationFilter->GetSize() );
  bspliner->SetSpacing( registrationFilter->GetSpacing() );

  typename ControlPointFilterType::OutputImageType::DirectionType direction;
  direction.SetIdentity();

  bspliner->SetDirection( direction );
  bspliner->Update();

  typedef itk::VectorImageFileWriter<
    typename RegistrationFilterType::DeformationFieldType,
    typename RegistrationFilterType::RealImageType> DeformationFieldWriterType;
  typename DeformationFieldWriterType::Pointer writer
    = DeformationFieldWriterType::New();
  writer->SetInput( bspliner->GetOutput() );
  writer->SetFileName( fileName.c_str() );
  writer->Update();
  }

  // Write control points to component nifti images.

  {
  std::string fileName = parser->GetOption( "output" )->GetValue()
    + std::string( "WarpControlPointLattice.nii.gz" );

  typedef itk::VectorImageFileWriter<
    typename RegistrationFilterType::ControlPointLatticeType,
    typename RegistrationFilterType::RealImageType> DeformationFieldWriterType;
  typename DeformationFieldWriterType::Pointer writer
    = DeformationFieldWriterType::New();
  writer->SetInput(
    registrationFilter->GetTotalDeformationFieldControlPoints() );
  writer->SetFileName( fileName.c_str() );
  writer->Update();
  }

  return EXIT_SUCCESS;

}

void InitializeCommandLineOptions( itk::CommandLineParser *parser )
{
  typedef itk::CommandLineParser::OptionType OptionType;

  {
  std::string description =
    std::string( "[fixedPointSet,movingPointSet," ) +
    std::string( "<numberOfFixedSamples>,<numberOfMovingSamples>]" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "point-sets" );
  option->SetShortName( 'p' );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "[splineOrder,meshResolution or isotropicSpacing,<calculateInitialSimilarityTransform>," ) +
    std::string( "<domainImage>,<initialControlPointLattice>,<directionality>]" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "transformation" );
  option->SetShortName( 't' );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description = std::string(
    "[<whichLabels>,<labelPercentages>,<labelGradientWeights>]" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "labels" );
  option->SetShortName( 'l' );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
  std::string description =
    std::string( "[alpha,annealingRate,pointSetSigma," ) +
    std::string( "<useRegularizationTerm>,<evaluationKNeighborhood>," ) +
    std::string( "<useAnisotropicCovariances>,<covarianceKNeighborhood>," ) +
    std::string( "<kernelSigma>]" );

  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "similarity" );
  option->SetShortName( 's' );
  option->SetDescription( description );
  parser->AddOption( option );
  }

  {
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
  OptionType::Pointer option = OptionType::New();
  option->SetLongName( "output" );
  option->SetShortName( 'o' );
  option->SetDescription( "filePrefix" );
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
    std::cout << "Usage: " << argv[0]
      << " PointSetDimension args" << std::endl;
    exit( 1 );
    }

  itk::CommandLineParser::Pointer parser = itk::CommandLineParser::New();
  parser->SetCommand( argv[0] );
  parser->SetCommandDescription( "DMFFD Labeled Point-Set Registration" );
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
     itkDMFFDLabeledPointSetRegistrationFilterTest<2>( parser );
     break;
     }
   case 3:
     {
     itkDMFFDLabeledPointSetRegistrationFilterTest<3>( parser );
     break;
     }
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

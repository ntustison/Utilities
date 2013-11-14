// General includes
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include "itksys/SystemTools.hxx"

// ITK includes
#include "itkNumericTraits.h"
#include "itkTimeProbe.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkPolyLineParametricPath.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkArrivalFunctionToPathFilter.h"
#include "itkSpeedFunctionToPathFilter.h"
#include "itkPathIterator.h"
#include "itkGradientDescentOptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkIterateNeighborhoodOptimizer.h"
#include "itkTubeSpatialObject.h"
#include "itkTubeSpatialObjectPoint.h"
#include "itkSpatialObjectPoint.h"
#include "itkSpatialObjectWriter.h"


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

/////////////////////////////////////////////////////////////
// Template for SpeedToPath with RegularStepGradientDescentOptimizer
template <int VDimension>
int SpeedToPath_RegularStepGradientDescent_ND(int argc, char* argv[])
{
		const unsigned int Dimension = VDimension;
		typedef float PixelType;
		typedef unsigned char OutputPixelType;
		typedef itk::Image< PixelType, Dimension > ImageType;
		typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
		typedef itk::ImageFileReader< ImageType > ReaderType;
		typedef itk::ImageFileWriter< OutputImageType > WriterType;
		typedef itk::PolyLineParametricPath< Dimension > PathType;
		typedef itk::SpeedFunctionToPathFilter< ImageType, PathType > PathFilterType;
		typedef typename PathFilterType::CostFunctionType::CoordRepType CoordRepType;
		typedef itk::PathIterator< OutputImageType, PathType > PathIteratorType;

		// Get arguments
		unsigned int argi = 1;
		char* OutputFilename = argv[argi++];
		char* SpeedFilename = argv[argi++];
		float TerminationValue = atof( argv[argi++] );
		unsigned int NumberOfIterations = atoi( argv[argi++] );
		float StepLengthFactor = atof( argv[argi++] );
		float StepLengthRelax = atof( argv[argi++] );

		// Read speed function
		typename ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( SpeedFilename );
		reader->Update();
		typename ImageType::Pointer speed = reader->GetOutput();
		speed->DisconnectPipeline();

		// Compute the minimum spacing
		typename ImageType::SpacingType spacing = speed->GetSpacing();
		double minspacing = spacing[0];
		for (unsigned int dim=0; dim<Dimension; dim++)
						if (spacing[dim] < minspacing) minspacing = spacing[dim];

		// Create Interpolator
		typedef itk::LinearInterpolateImageFunction<ImageType, CoordRepType>
				InterpolatorType;
		typename InterpolatorType::Pointer interp = InterpolatorType::New();

		// Create Cost Function
		typename PathFilterType::CostFunctionType::Pointer cost =
						PathFilterType::CostFunctionType::New();
		cost->SetInterpolator( interp );

		// Create RegularStepGradientDescentOptimizer
		typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
		typename OptimizerType::Pointer optimizer = OptimizerType::New();
		optimizer->SetNumberOfIterations( NumberOfIterations );
		optimizer->SetMaximumStepLength( 1.0*StepLengthFactor*minspacing );
		optimizer->SetMinimumStepLength( 0.5*StepLengthFactor*minspacing );
		optimizer->SetRelaxationFactor( StepLengthRelax );

		// Create path filter
		typename PathFilterType::Pointer pathFilter = PathFilterType::New();
		pathFilter->SetInput( speed );
		pathFilter->SetCostFunction( cost );
		pathFilter->SetOptimizer( optimizer );
		pathFilter->SetTerminationValue( TerminationValue );

  typename PathFilterType::PathInfo info;
  for( unsigned int n = argi; n < argc; n++ )
    {
    std::vector<unsigned int> index = ConvertVector<unsigned int>( std::string( argv[n] ) );

    typename ImageType::IndexType imageIndex;
    for( unsigned int d = 0; d < VDimension; d++ )
      {
      imageIndex[d] = index[d];
      }
    typename ImageType::PointType imagePoint;
    speed->TransformIndexToPhysicalPoint( imageIndex, imagePoint );

    typename PathFilterType::PointType pathPoint;
    pathPoint.CastFrom( imagePoint );

    if( n == argi )
      {
      info.SetStartPoint( pathPoint );
      }
    else if( n == argc - 1 )
      {
      info.SetEndPoint( pathPoint );
      }
    else
      {
      info.AddWayPoint( pathPoint );
      }
    }
  pathFilter->AddPathInfo( info );

		// Compute the path
		std::cout << "Computing path..." << std::endl;
		itk::TimeProbe time;
		time.Start( );
		pathFilter->Update( );
		time.Stop( );
		std::cout << std::setprecision(3) << "Path computed in: " << time.GetMean() << " seconds" << std::endl;

		// Allocate output image
		typename OutputImageType::Pointer output = OutputImageType::New();
		output->SetRegions( speed->GetLargestPossibleRegion() );
		output->SetSpacing( speed->GetSpacing() );
		output->SetOrigin( speed->GetOrigin() );
		output->Allocate( );
		output->FillBuffer( itk::NumericTraits<OutputPixelType>::Zero );

		// Rasterize path
		for (unsigned int i=0; i<pathFilter->GetNumberOfOutputs(); i++)
		  {
				// Get the path
				typename PathType::Pointer path = pathFilter->GetOutput( i );

				// Check path is valid
				if ( path->GetVertexList()->Size() == 0 )
  				{
						std::cout << "WARNING: Path " << (i+1) << " contains no points!" << std::endl;
						continue;
		  		}

				// Iterate path and convert to image
				std::cout << "Rasterizing path..." << std::endl;
				PathIteratorType it( output, path );
				for (it.GoToBegin(); !it.IsAtEnd(); ++it)
				  {
						it.Set( itk::NumericTraits<OutputPixelType>::One );
				  }
		  }

		// Write output
		std::cout << "Output: " << OutputFilename << std::endl;
		typename WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( OutputFilename );
		writer->SetInput( output );
		writer->Update();

		//Return
		return EXIT_SUCCESS;
}

/////////////////////////////////////////////////////////////
// Template for SpeedToPath with IterateNeighborhoodOptimizer
template <int VDimension>
int SpeedToPath_IterateNeighborhood_ND(int argc, char* argv[])
{
 	const unsigned int Dimension = VDimension;
	 typedef float PixelType;
		typedef unsigned char OutputPixelType;
		typedef itk::Image< PixelType, Dimension > ImageType;
		typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
		typedef itk::ImageFileReader< ImageType > ReaderType;
		typedef itk::ImageFileWriter< OutputImageType > WriterType;
		typedef itk::PolyLineParametricPath< Dimension > PathType;
		typedef itk::SpeedFunctionToPathFilter< ImageType, PathType > PathFilterType;
		typedef typename PathFilterType::CostFunctionType::CoordRepType CoordRepType;
		typedef itk::PathIterator< OutputImageType, PathType > PathIteratorType;

		// Get arguments
		unsigned int argi = 1;
		char* OutputFilename = argv[argi++];
		char* SpeedFilename = argv[argi++];
		float TerminationValue = atof( argv[argi++] );
		float StepLengthFactor = atof( argv[argi++] );
		// NOTE: Points will be read from the command line later

		// Read speed function
		std::cout << "Speed: " << SpeedFilename << std::endl;
		typename ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( SpeedFilename );
		reader->Update();
		typename ImageType::Pointer speed = reader->GetOutput();
		speed->DisconnectPipeline();

		// Create Interpolator
		typedef itk::LinearInterpolateImageFunction<ImageType, CoordRepType>
				InterpolatorType;
		typename InterpolatorType::Pointer interp = InterpolatorType::New();

		// Create Cost Function
		typename PathFilterType::CostFunctionType::Pointer cost =
						PathFilterType::CostFunctionType::New();
		cost->SetInterpolator( interp );

		// Create IterateNeighborhoodOptimizer
		typedef itk::IterateNeighborhoodOptimizer OptimizerType;
		typename OptimizerType::Pointer optimizer = OptimizerType::New();
		optimizer->MinimizeOn( );
		optimizer->FullyConnectedOn( );
		typename OptimizerType::NeighborhoodSizeType size( Dimension );
		for (unsigned int i=0; i<Dimension; i++)
				size[i] = speed->GetSpacing()[i] * StepLengthFactor;
		optimizer->SetNeighborhoodSize( size );

		// Create path filter
		typename PathFilterType::Pointer pathFilter = PathFilterType::New();
		pathFilter->SetInput( speed );
		pathFilter->SetCostFunction( cost );
		pathFilter->SetOptimizer( optimizer );
		pathFilter->SetTerminationValue( TerminationValue );

  typename PathFilterType::PathInfo info;
  for( unsigned int n = argi; n < argc; n++ )
    {
    std::vector<unsigned int> index = ConvertVector<unsigned int>( std::string( argv[n] ) );

    typename ImageType::IndexType imageIndex;
    for( unsigned int d = 0; d < VDimension; d++ )
      {
      imageIndex[d] = index[d];
      }
    typename ImageType::PointType imagePoint;
    speed->TransformIndexToPhysicalPoint( imageIndex, imagePoint );

    typename PathFilterType::PointType pathPoint;
    pathPoint.CastFrom( imagePoint );

    if( n == argi )
      {
      info.SetStartPoint( pathPoint );
      }
    else if( n == argc - 1 )
      {
      info.SetEndPoint( pathPoint );
      }
    else
      {
      info.AddWayPoint( pathPoint );
      }
    }
  pathFilter->AddPathInfo( info );

		// Compute the path
		std::cout << "Computing path..." << std::endl;
		itk::TimeProbe time;
		time.Start( );
		pathFilter->Update( );
		time.Stop( );
		std::cout << std::setprecision(3) << "Path computed in: " << time.GetMean() << " seconds" << std::endl;

		// Allocate output image
		typename OutputImageType::Pointer output = OutputImageType::New();
		output->SetRegions( speed->GetLargestPossibleRegion() );
		output->SetSpacing( speed->GetSpacing() );
		output->SetOrigin( speed->GetOrigin() );
		output->Allocate( );
		output->FillBuffer( itk::NumericTraits<OutputPixelType>::Zero );

// Rasterize path
		for (unsigned int i=0; i<pathFilter->GetNumberOfOutputs(); i++)
  		{
				// Get the path
				typename PathType::Pointer path = pathFilter->GetOutput( i );

				// Check path is valid
				if ( path->GetVertexList()->Size() == 0 )
				  {
						std::cout << "WARNING: Path " << (i+1) << " contains no points!" << std::endl;
						continue;
		  		}

				// Iterate path and convert to image
				std::cout << "Rasterizing path..." << std::endl;
				PathIteratorType it( output, path );
				for (it.GoToBegin(); !it.IsAtEnd(); ++it)
				  {
						it.Set( itk::NumericTraits<OutputPixelType>::One );
				  }
		}

		// Write output
		std::cout << "Output: " << OutputFilename << std::endl;
		typename WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( OutputFilename );
		writer->SetInput( output );
		writer->Update();

		//Return
		return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if( argc < 3 )
    {
    std::cout << "Usage: "<< argv[0]
      << " imageDimension optimizerType optimizerTypeSpecificArguments" << std::endl;
    std::cout << "  optimizerType = 1 (IterateNeighborhood)" << std::endl;
    std::cout << "  optimizerType = 2 (RegularGradientDescent)" << std::endl;
    exit( 1 );
    }

  int optimizerType = atoi( argv[2] );

		// Print header info
		if ( optimizerType == 1 && argc < 7 )
		  {
				std::cerr << "Usage: " << std::endl;
				std::cerr << argv[0];
				std::cerr << " OutputFilename";
				std::cerr << " SpeedFilename";
				std::cerr << " TerminationValue(Good default = 2.0)";    // Good default = 2.0
				std::cerr << " StepLengthFactor(Good default = 1.0)";    // Good default = 1.0
				std::cerr << " index0 ... indexN";
				std::cerr << std::endl;
				return EXIT_FAILURE;
  		}
  else if ( optimizerType == 2 && argc < 9 )
    {
				std::cerr << "Usage: " << std::endl;
				std::cerr << argv[0];
				std::cerr << " OutputFilename";
				std::cerr << " SpeedFilename";
				std::cerr << " TerminationValue(Good default=2.0)";    // Good default = 2.0
				std::cerr << " NumberOfIterations(Good default=1000)";  // Good default = 1000
				std::cerr << " StepLengthFactor(Good default=1.0)";    // Good default = 1.0
				std::cerr << " StepLengthRelax(Good default=0.999)";     // Good default = 0.999
				std::cerr << " index0 ... indexN";
				std::cerr << std::endl;
				return EXIT_FAILURE;
    }


  switch( atoi( argv[1] ) )
   {
   case 2:
     if( optimizerType == 1 )
       {
       return SpeedToPath_RegularStepGradientDescent_ND<2>( argc, argv );
       }
     else
       {
       return SpeedToPath_IterateNeighborhood_ND<2>( argc, argv );
       }
     break;
   case 3:
     if( optimizerType == 1 )
       {
       return SpeedToPath_RegularStepGradientDescent_ND<3>( argc, argv );
       }
     else
       {
       return SpeedToPath_IterateNeighborhood_ND<3>( argc, argv );
       }
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}







// #include "itkImageFileReader.h"
// #include "itkImageFileWriter.h"
// #include "itkPathIterator.h"
// #include "itkMinimalPathImageFunction.h"
//
// #include <string>
// #include <vector>
//
// #include <fstream>
//
// template<class TValue>
// TValue Convert( std::string optionString )
// {
// 		TValue value;
// 		std::istringstream iss( optionString );
// 		iss >> value;
// 		return value;
// }
//
// template<class TValue>
// std::vector<TValue> ConvertVector( std::string optionString )
// {
// 		std::vector<TValue> values;
// 		std::string::size_type crosspos = optionString.find( 'x', 0 );
//
// 		if ( crosspos == std::string::npos )
// 				{
// 				values.push_back( Convert<TValue>( optionString ) );
// 				}
// 		else
// 				{
// 				std::string element = optionString.substr( 0, crosspos ) ;
// 				TValue value;
// 				std::istringstream iss( element );
// 				iss >> value;
// 				values.push_back( value );
// 				while ( crosspos != std::string::npos )
// 						{
// 						std::string::size_type crossposfrom = crosspos;
// 						crosspos = optionString.find( 'x', crossposfrom + 1 );
// 						if ( crosspos == std::string::npos )
// 								{
// 								element = optionString.substr( crossposfrom + 1, optionString.length() );
// 								}
// 						else
// 								{
// 								element = optionString.substr( crossposfrom + 1, crosspos ) ;
// 								}
// 						std::istringstream iss( element );
// 						iss >> value;
// 						values.push_back( value );
// 						}
// 				}
// 		return values;
//   }
//
// template<unsigned int ImageDimension>
// int MinimalPath( unsigned int argc, char *argv[] )
// {
//   typedef float PixelType;
//   typedef unsigned char PathPixelType;
//
//   typedef itk::Image<PixelType, ImageDimension> ImageType;
//   typedef itk::Image<PathPixelType, ImageDimension> PathImageType;
//   typedef typename ImageType::IndexType IndexType;
//
//   typedef itk::ImageFileReader<ImageType> ReaderType;
//   typename ReaderType::Pointer reader = ReaderType::New();
//   reader->SetFileName( argv[2] );
//   reader->Update();
//
//   typename PathImageType::Pointer pathImage = PathImageType::New();
//   pathImage->SetDirection( reader->GetOutput()->GetDirection() );
//   pathImage->SetOrigin( reader->GetOutput()->GetOrigin() );
//   pathImage->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
//   pathImage->SetSpacing( reader->GetOutput()->GetSpacing() );
//   pathImage->Allocate();
//   pathImage->FillBuffer( itk::NumericTraits<PathPixelType>::Zero );
//
//   std::vector<unsigned int> anchor
//     = ConvertVector<unsigned int>( std::string( argv[4] ) );
//   std::vector<unsigned int> free
//     = ConvertVector<unsigned int>( std::string( argv[5] ) );
//
//   IndexType anchorIndex;
//   IndexType freeIndex;
//   for( unsigned int d = 0; d < ImageDimension; d++ )
//     {
//     anchorIndex[d] = anchor[d];
//     freeIndex[d] = free[d];
//     }
//
//   typedef itk::MinimalPathImageFunction<ImageType> FunctionType;
//   typename FunctionType::Pointer function = FunctionType::New();
//   function->SetUseFaceConnectedness( argc > 6 ?
//     static_cast<bool>( atoi( argv[6] ) ) : true );
//   function->SetUseImageSpacing( true );
//   function->SetInputImage( reader->GetOutput() );
//   function->SetAnchorSeed( anchorIndex );
//
//   std::string txtFileName = std::string( argv[3] ) + std::string( ".txt" );
//   std::ofstream str( txtFileName.c_str() );
//   str << "0 0 0 0" << std::endl;
//
//   typedef itk::PathIterator<ImageType,
//     typename FunctionType::OutputType> IteratorType;
//
//   typename FunctionType::OutputType::Pointer path
//     = function->EvaluateAtIndex( freeIndex );
//   IteratorType It( reader->GetOutput(), path );
//   It.GoToBegin();
//   while ( !It.IsAtEnd() )
//     {
//     typename PathImageType::PointType point;
//     pathImage->TransformIndexToPhysicalPoint( It.GetIndex(), point );
//     str << point[0] << " " << point[1];
//     if( ImageDimension == 3 )
//       {
//       str << " " << point[2];
//       }
//     else
//       {
//       str << " 0";
//       }
//     str << " 1" << std::endl;
//
//     pathImage->SetPixel( It.GetIndex(),
//       itk::NumericTraits<PathPixelType>::One );
//     ++It;
//     }
//   str << "0 0 0 0" << std::endl;
//   str.close();
//
//
//   std::string imageFileName = std::string( argv[3] ) + std::string( ".nii.gz" );
//   typedef itk::ImageFileWriter<PathImageType> WriterType;
//   typename WriterType::Pointer writer = WriterType::New();
//   writer->SetFileName( imageFileName.c_str() );
//   writer->SetInput( pathImage );
//   writer->Update();
//
//   return EXIT_SUCCESS;
// }
//
// int main( int argc, char *argv[] )
// {
//   if( argc < 6 )
//     {
//     std::cout << "Usage: "<< argv[0]
//       << " imageDimension inputImage outputPrefix anchorIndex"
//       << " freeIndex [useFaceConnectedness]" << std::endl;
//
//     exit( 1 );
//     }
//
//   switch( atoi( argv[1] ) )
//    {
//    case 2:
//      MinimalPath<2>( argc, argv );
//      break;
//    case 3:
//      MinimalPath<3>( argc, argv );
//      break;
//    default:
//       std::cerr << "Unsupported dimension" << std::endl;
//       exit( EXIT_FAILURE );
//    }
// }
//

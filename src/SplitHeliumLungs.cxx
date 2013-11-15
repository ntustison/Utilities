// ITK includes
#include "itkNumericTraits.h"
#include "itkExtractImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkPolyLineParametricPath.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkArrivalFunctionToPathFilter.h"
#include "itkSpeedFunctionToPathFilter.h"
#include "itkPathIterator.h"
#include "itkGradientDescentOptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkIterateNeighborhoodOptimizer.h"

// General includes
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include "itksys/SystemTools.hxx"

int main( int argc, char *argv[] )
{
  if( argc < 3 )
    {
    std::cout << "Usage: "<< argv[0] << " inputSpeedImage outputImage [stepLengthFactor=0.1] [terminationValue=2.0]" << std::endl;
    exit( 1 );
    }

  const unsigned int ImageDimension = 3;

  float StepLengthFactor = 0.1;
  if( argc > 3 )
    {
    StepLengthFactor = atof( argv[3] );
    }
  float TerminationValue = 2.0;
  if( argc > 4 )
    {
    TerminationValue = atof( argv[4] );
    }

  typedef itk::Image<float, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  reader->Update();

  ImageType::Pointer outputImage = ImageType::New();
  outputImage->CopyInformation( reader->GetOutput() );
  outputImage->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  outputImage->Allocate();
  outputImage->FillBuffer( 0 );

  typedef itk::Image<unsigned int, 2> OutputSliceType;
  typedef itk::Image<float, 2> SliceType;

		typedef itk::PolyLineParametricPath<2> PathType;
		typedef itk::SpeedFunctionToPathFilter<SliceType, PathType> PathFilterType;
		typedef PathFilterType::CostFunctionType::CoordRepType CoordRepType;
		typedef itk::PathIterator<OutputSliceType, PathType> PathIteratorType;

  ImageType::RegionType region;
  ImageType::RegionType::SizeType imageSize
    = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
  ImageType::IndexType imageIndex
    = reader->GetOutput()->GetLargestPossibleRegion().GetIndex();

  ImageType::IndexType index;
  ImageType::RegionType::SizeType size;
  size[0] = imageSize[0];
  size[1] = imageSize[1];
  size[2] = 0;
  region.SetSize( size );
  index[0] = imageIndex[0];
  index[1] = imageIndex[1];

  long lastIndex = imageIndex[2] + imageSize[2]-1;

  for( long n = imageIndex[2]; n < lastIndex; n++ )
    {
    std::cout << "Processing slice " << n-imageIndex[2]+1 << " out of " << imageSize[2] - 1 << std::endl;

    index[2] = n;
    region.SetIndex( index );

    typedef itk::ExtractImageFilter<ImageType, SliceType> ExtracterType;
    ExtracterType::Pointer extracter = ExtracterType::New();
    extracter->SetInput( reader->GetOutput() );
    extracter->SetDirectionCollapseToIdentity();
    extracter->SetExtractionRegion( region );

    SliceType::Pointer speed = extracter->GetOutput();
    speed->Update();
    speed->DisconnectPipeline();

    // Compute the minimum spacing
    SliceType::SpacingType spacing = speed->GetSpacing();
    double minspacing = vnl_math_min( spacing[0], spacing[1] );

    // Create Interpolator
    typedef itk::LinearInterpolateImageFunction<SliceType, CoordRepType>
      InterpolatorType;
    InterpolatorType::Pointer interp = InterpolatorType::New();

    // Create Cost Function
    PathFilterType::CostFunctionType::Pointer cost =
        PathFilterType::CostFunctionType::New();
    cost->SetInterpolator( interp );

    // Create IterateNeighborhoodOptimizer
    typedef itk::IterateNeighborhoodOptimizer OptimizerType;
    OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->MinimizeOn( );
    optimizer->FullyConnectedOn( );
    OptimizerType::NeighborhoodSizeType size( 2 );
    for (unsigned int i=0; i<2; i++)
      size[i] = speed->GetSpacing()[i] * StepLengthFactor;
    optimizer->SetNeighborhoodSize( size );

    // Create path filter
    PathFilterType::Pointer pathFilter = PathFilterType::New();
    pathFilter->SetInput( speed );
    pathFilter->SetCostFunction( cost );
    pathFilter->SetOptimizer( optimizer );
    pathFilter->SetTerminationValue( TerminationValue );

    // add start and end points
    SliceType::IndexType anchor;
    SliceType::IndexType free;

    PathFilterType::PathInfo info;

    anchor[0] = static_cast<long>(
      vcl_floor( imageIndex[0] + 0.5*imageSize[0] + 0.5 ) );
    anchor[1] = imageIndex[1];

    SliceType::PointType anchorPoint;
    speed->TransformIndexToPhysicalPoint( anchor, anchorPoint );
    PathFilterType::PointType anchorPathPoint;
    anchorPathPoint.CastFrom( anchorPoint );
    info.SetStartPoint( anchorPathPoint );

    free[0]   = static_cast<long>(
      vcl_floor( imageIndex[0] + 0.5*imageSize[0] + 0.5 ) );
    free[1]   = imageIndex[1] + imageSize[1]-1;

    SliceType::PointType freePoint;
    speed->TransformIndexToPhysicalPoint( free, freePoint );
    PathFilterType::PointType freePathPoint;
    freePathPoint.CastFrom( freePoint );
    info.SetEndPoint( freePathPoint );

    pathFilter->AddPathInfo( info );
    pathFilter->Update( );

    OutputSliceType::Pointer output = OutputSliceType::New();
    output->SetRegions( speed->GetLargestPossibleRegion() );
    output->SetSpacing( speed->GetSpacing() );
    output->SetOrigin( speed->GetOrigin() );
    output->Allocate( );
    output->FillBuffer( 0 );

    for (unsigned int i=0; i<pathFilter->GetNumberOfOutputs(); i++)
      {
      // Get the path
      PathType::Pointer path = pathFilter->GetOutput( i );

      // Check path is valid
      if ( path->GetVertexList()->Size() == 0 )
        {
        std::cout << "WARNING: Path " << (i+1) << " contains no points!" << std::endl;
        continue;
        }

      // Iterate path and convert to image
      PathIteratorType it( output, path );
      for (it.GoToBegin(); !it.IsAtEnd(); ++it)
        {
        it.Set( 1 );
        }
      }

    itk::ImageRegionIteratorWithIndex<OutputSliceType> ItS( output, output->GetLargestPossibleRegion() );
    for( ItS.GoToBegin(); !ItS.IsAtEnd(); ++ItS )
      {
      if( ItS.Get() > 0 )
        {
        OutputSliceType::IndexType sliceIndex = ItS.GetIndex();
        ImageType::IndexType outputIndex;
        outputIndex[0] = sliceIndex[0];
        outputIndex[1] = sliceIndex[1];
        outputIndex[2] = n;
        outputImage->SetPixel( outputIndex, 1 );
        }
      }
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( outputImage );
  writer->SetFileName( argv[2] );
  writer->Update();

  return EXIT_SUCCESS;
}


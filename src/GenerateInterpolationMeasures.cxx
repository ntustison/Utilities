#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkLabelPerimeterEstimationCalculator.h"

#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>


template <unsigned int ImageDimension>
int GenerateInterpolationMeasures( int argc, char *argv[] )
{
  typedef int PixelType;
  typedef float RealType;

  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Image<PixelType, ImageDimension> LabelImageType;
  typedef itk::ImageFileReader<LabelImageType> ReaderType;

  typename ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[2] );
  reader1->Update();

  typename ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[3] );
  reader2->Update();

  typedef itk::LabelGeometryImageFilter<LabelImageType, RealImageType> FilterType;

  typename FilterType::Pointer filter1 = FilterType::New();
  filter1->SetInput( reader1->GetOutput() );
  filter1->CalculatePixelIndicesOff();
  filter1->CalculateOrientedBoundingBoxOff();
  filter1->CalculateOrientedLabelRegionsOff();
  filter1->Update();

  typename FilterType::Pointer filter2 = FilterType::New();
  filter2->SetInput( reader2->GetOutput() );
  filter2->CalculatePixelIndicesOff();
  filter2->CalculateOrientedBoundingBoxOff();
  filter2->CalculateOrientedLabelRegionsOff();
  filter2->Update();

  typedef itk::LabelPerimeterEstimationCalculator<LabelImageType> AreaFilterType;

  typename AreaFilterType::Pointer areaFilter1 = AreaFilterType::New();
  areaFilter1->SetImage( reader1->GetOutput() );
  areaFilter1->Compute();

  typename AreaFilterType::Pointer areaFilter2 = AreaFilterType::New();
  areaFilter2->SetImage( reader2->GetOutput() );
  areaFilter2->Compute();


  std::cout << std::left << std::setw( 7 )  << "Label"
            << std::left << std::setw( 20 ) << "PercentVolumeDiff"
            << std::left << std::setw( 25 ) << "PercentSurfAreaDiff"
            << std::left << std::setw( 25 ) << "Distance";
  std::cout << std::endl;


  typename FilterType::LabelsType allLabels1 = filter1->GetLabels();
  std::sort( allLabels1.begin(), allLabels1.end() );

  typename FilterType::LabelsType allLabels2 = filter2->GetLabels();

  typename FilterType::LabelsType::iterator allLabelsIt;
  for( allLabelsIt = allLabels1.begin(); allLabelsIt != allLabels1.end(); allLabelsIt++ )
    {
    if( *allLabelsIt == 0 || std::find( allLabels2.begin(), allLabels2.end(), *allLabelsIt ) == allLabels2.end() )
      {
      continue;
      }
    RealType volume = filter1->GetVolume( *allLabelsIt );
    RealType volumeRotated = filter2->GetVolume( *allLabelsIt );

    RealType percentVolumeDifference = ( volume - volumeRotated ) / volume;

    RealType surfaceArea = areaFilter1->GetPerimeter( *allLabelsIt );
    RealType surfaceAreaRotated = areaFilter2->GetPerimeter( *allLabelsIt );

    RealType percentSurfaceAreaDifference = ( surfaceArea - surfaceAreaRotated ) / surfaceArea;

    typename LabelImageType::PointType centroidPoint = filter1->GetCentroid( *allLabelsIt );
    typename LabelImageType::PointType centroidPointRotated = filter2->GetCentroid( *allLabelsIt );

    RealType distance = centroidPoint.EuclideanDistanceTo( centroidPointRotated );

    std::cout << std::left << std::setw( 7 )  << *allLabelsIt
              << std::left << std::setw( 20 ) << percentVolumeDifference
              << std::left << std::setw( 25 ) << percentSurfaceAreaDifference
              << std::left << std::setw( 25 ) << distance;
    std::cout << std::endl;
    }
  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension labelImage labelImageRotated " << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     GenerateInterpolationMeasures<2>( argc, argv );
     break;
   case 3:
     GenerateInterpolationMeasures<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


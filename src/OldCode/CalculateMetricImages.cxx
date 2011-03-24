#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNeighborhoodIterator.h"
#include "itkVector.h"
#include "itkVectorImageFileReader.h"
#include "itkVectorImageFileWriter.h"

// Image similarity metrics
#include "itkAvantsMutualInformationRegistrationFunction.h"
#include "itkMeanSquareRegistrationFunction.h"
#include "itkPDEDeformableRegistrationFunction.h"
#include "itkProbabilisticRegistrationFunction.h"
#include "itkRobustDemonsRegistrationFunction.h"
#include "itkRobustOpticalFlow.h"
#include "itkSectionMutualInformationRegistrationFunction.h"

#include "vnl/vnl_math.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << "Usage: " << argv[0]
              << " fixedImage movingImage whichMetric outputMetricImage" << std::endl;
   std::cout << "  Cost/similarity Options: " << std::endl;
   std::cout << "     0: optical flow" << std::endl;
   std::cout << "     1: n.a." << std::endl;
   std::cout << "     2: Mutual Information (not good)" << std::endl;
   std::cout << "     3: MutualInformation (good - 2nd most recommended) " << std::endl;
   std::cout << "     4: Cross Correlation of radius 5 (recommended) "     << std::endl;
   std::cout << "     5: Cross Correlation of radius 2 (similar to 4) "    << std::endl;
    exit( 1 );
    }
  typedef float                                                   RealType;

  typedef itk::Image<RealType, ImageDimension>                    ImageType;
  typedef itk::Image<RealType, ImageDimension>                    RealImageType;
  typedef itk::ImageFileReader<ImageType>                         ImageReaderType;
  typedef itk::Vector<RealType, ImageDimension>                   VectorType;
  typedef itk::Image<VectorType, ImageDimension>                  DeformationFieldType;

  ImageReaderType::Pointer fixedImage = ImageReaderType::New();
  fixedImage->SetFileName( argv[1] );
  fixedImage->Update();

  ImageReaderType::Pointer movingImage = ImageReaderType::New();
  movingImage->SetFileName( argv[2] );
  movingImage->Update();

  typedef itk::PDEDeformableRegistrationFunction<RealImageType, RealImageType, DeformationFieldType> MetricType;
  MetricType::Pointer metric;
  MetricType::RadiusType radius;
  bool maximizeMetric;

  unsigned int numberOfMISamples = 7500;
  unsigned int numberOfHistogramBins = 64;
  if ( fixedImage->GetOutput()->GetLargestPossibleRegion().GetSize()[0] < 80 || ImageDimension == 2 )
    {
    numberOfHistogramBins = 32;
    }
  else if ( fixedImage->GetOutput()->GetLargestPossibleRegion().GetSize()[0] > 256 )
    {
    numberOfHistogramBins = 128;
    }

  switch( atoi( argv[3] ) )
    {
    case 0:
      {
      typedef itk::MeanSquareRegistrationFunction
        <RealImageType, RealImageType, DeformationFieldType> Metric0Type;
       Metric0Type::Pointer metric0 = Metric0Type::New();
      metric0->SetIntensityDifferenceThreshold( 0 );
      metric0->SetGradientStep( 1e6 );
      metric0->SetRobust( true );
      metric0->SetSymmetric( false );
      metric0->SetNormalizeGradient( false );
      maximizeMetric = false;
      radius.Fill( 0 );
      metric = metric0;
      break;
      }
    case 1:
      {
      typedef itk::ProbabilisticRegistrationFunction
        <RealImageType, RealImageType, DeformationFieldType> Metric1Type;
       Metric1Type::Pointer metric1 = Metric1Type::New();
      metric1->SetFullyRobust( true );
      metric1->SetNormalizeGradient( false );
      metric1->SetGradientStep( 1e6 );
      radius.Fill( 2 );
      maximizeMetric = true;
      metric = metric1;
      break;
      }
    case 2:
      {
      typedef itk::AvantsMutualInformationRegistrationFunction
        <RealImageType, RealImageType, DeformationFieldType> Metric2Type;
       Metric2Type::Pointer metric2 = Metric2Type::New();
      metric2->SetNumberOfSpatialSamples( numberOfMISamples );
      metric2->SetNumberOfHistogramBins( numberOfHistogramBins );
      metric2->SetNormalizeGradient( true );
      metric2->SetGradientStep( 1e6 );
      radius.Fill( 1 );
      maximizeMetric = true;
      metric = metric2;
      break;
      }
    case 3:
      {
      typedef itk::AvantsMutualInformationRegistrationFunction
        <RealImageType, RealImageType, DeformationFieldType> Metric3Type;
       Metric3Type::Pointer metric3 = Metric3Type::New();
      metric3->SetNumberOfSpatialSamples( numberOfMISamples );
      metric3->SetNumberOfHistogramBins( numberOfHistogramBins );
      metric3->SetNormalizeGradient( false );
      metric3->SetGradientStep( 1e5 );
      radius.Fill( 1 );
      maximizeMetric = true;
      metric = metric3;
      break;
      }
    case 4:
      {
      typedef itk::ProbabilisticRegistrationFunction
        <RealImageType, RealImageType, DeformationFieldType> Metric4Type;
       Metric4Type::Pointer metric4 = Metric4Type::New();
      metric4->SetNormalizeGradient( false );
      metric4->SetGradientStep( 1e6 );
      radius.Fill( 5 );
      maximizeMetric = true;
      metric = metric4;
      break;
      }
    case 5:
      {
      typedef itk::ProbabilisticRegistrationFunction
        <RealImageType, RealImageType, DeformationFieldType> Metric5Type;
       Metric5Type::Pointer metric5 = Metric5Type::New();
      metric5->SetNormalizeGradient( false );
      metric5->SetGradientStep( 1e6 );
      radius.Fill( 2 );
      maximizeMetric = true;
      metric = metric5;
      break;
      }
    case 6: default:
      {
      typedef itk::RobustOpticalFlow
        <RealImageType, RealImageType, DeformationFieldType> Metric6Type;
       Metric6Type::Pointer metric6 = Metric6Type::New();
      metric6->SetGradientStep( 1e6 );
      radius.Fill( 2 );
      maximizeMetric = false;
      metric = metric6;
      break;
      }
    case 7:
      {
      typedef itk::AvantsMutualInformationRegistrationFunction
        <RealImageType, RealImageType, DeformationFieldType> Metric7Type;
       Metric7Type::Pointer metric7 = Metric7Type::New();
      metric7->SetNumberOfSpatialSamples( numberOfMISamples );
      metric7->SetNumberOfHistogramBins( numberOfHistogramBins );
      metric7->SetOpticalFlow( true );
      metric7->SetGradientStep( 1e5 );
      maximizeMetric = true;
      radius.Fill( 1 );
      metric = metric7;
      break;
      }
    case 8:
      {
      typedef itk::SectionMutualInformationRegistrationFunction
        <RealImageType, RealImageType, DeformationFieldType> Metric8Type;
       Metric8Type::Pointer metric8 = Metric8Type::New();
      metric8->SetNumberOfSpatialSamples( 7000 );
      metric8->SetNumberOfHistogramBins( 26 );
      metric8->SetOpticalFlow( false );
      metric8->SetNormalizeGradient( true );
      metric8->SetZeroInZ( true );
      metric8->SetGradientStep( 1e2 );
      radius.Fill( 1 );
      maximizeMetric = true;
      metric = metric8;
      break;
      }
    }
  metric->SetRadius( radius );

  VectorType V;

  metric->SetRadius( radius );
  metric->SetMovingImage( movingImage->GetOutput() );
  metric->SetFixedImage( fixedImage->GetOutput() );
  metric->SetDeformationField( NULL );
  metric->InitializeIteration();

  ImageType::Pointer metricField = ImageType::New();
  metricField->SetOrigin( movingImage->GetOutput()->GetOrigin() );
  metricField->SetSpacing( movingImage->GetOutput()->GetSpacing() );
  metricField->SetRegions( movingImage->GetOutput()->GetLargestPossibleRegion() );
  metricField->Allocate();
  metricField->FillBuffer( 0 );

  DeformationFieldType::Pointer gradientField = DeformationFieldType::New();
  gradientField->SetOrigin( movingImage->GetOutput()->GetOrigin() );
  gradientField->SetSpacing( movingImage->GetOutput()->GetSpacing() );
  gradientField->SetRegions( movingImage->GetOutput()->GetLargestPossibleRegion() );
  gradientField->Allocate();
  V.Fill( 0 );
  gradientField->FillBuffer( V );

  itk::NeighborhoodIterator<DeformationFieldType>
    It( radius, gradientField, gradientField->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<ImageType> ItM
    ( metricField, metricField->GetLargestPossibleRegion() );

  RealType energy = 0.0;

  for ( ItM.GoToBegin(); !ItM.IsAtEnd(); ++ItM )
    {
    It.SetLocation( ItM.GetIndex() );
    metric->SetEnergy( 0.0 );
    VectorType grad = metric->ComputeUpdate( It, NULL );
    RealType value = metric->GetEnergy();
    if ( !vnl_math_isnan( value ) )
      {
      energy += value;
      ItM.Set( value );
      It.SetCenterPixel( grad );
      }
    }
  std::cout << "Image metric = " << energy << std::endl;

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[4] );
  writer->SetInput( metricField );
  writer->Update();

  typedef itk::VectorImageFileWriter<DeformationFieldType, ImageType> VectorImageWriterType;
  VectorImageWriterType::Pointer vectorImageWriter = VectorImageWriterType::New();
  vectorImageWriter->SetFileName( argv[4] );
  vectorImageWriter->SetInput( gradientField );
  vectorImageWriter->Update();
}

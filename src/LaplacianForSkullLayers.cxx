#include "itkBinaryThresholdImageFilter.h"
#include "itkBresenhamLine.h"
#include "itkGradientImageFilter.h"
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLabelContourImageFilter.h"
#include "itkNeighborhoodIterator.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkVariableSizeMatrix.h"
#include "itkVector.h"

#include <string>
#include <vector>


template <unsigned int ImageDimension>
int LaplacianLayers( int argc, char *argv[] )
{
  const unsigned int numberOfLayers = 3;

  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef unsigned int LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
  typename LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName( argv[3] );
  labelReader->Update();

  typedef itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( labelReader->GetOutput() );
  filter->SetLowerThreshold( 1 );
  filter->SetUpperThreshold( 3 );
  filter->SetInsideValue( 1 );
  filter->SetOutsideValue( 0 );
  filter->Update();

  typedef itk::LabelContourImageFilter<LabelImageType, LabelImageType> ContourFilterType;
  typename ContourFilterType::Pointer contours = ContourFilterType::New();
  contours->SetInput( labelReader->GetOutput() );
  contours->SetFullyConnected( true );
  contours->SetBackgroundValue( 0 );
  contours->Update();

  typename ContourFilterType::Pointer contours2 = ContourFilterType::New();
  contours2->SetInput( filter->GetOutput() );
  contours2->SetFullyConnected( true );
  contours2->SetBackgroundValue( 0 );
  contours2->Update();

  typename ImageType::Pointer laplacian = ImageType::New();
  laplacian->CopyInformation( labelReader->GetOutput() );
  laplacian->SetRegions( reader->GetOutput()->GetRequestedRegion() );
  laplacian->Allocate();
  laplacian->FillBuffer( 0 );

  itk::ImageRegionIteratorWithIndex<LabelImageType> ItF(
    contours2->GetOutput(), contours2->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<LabelImageType> ItC(
    contours->GetOutput(), contours->GetOutput()->GetLargestPossibleRegion() );
  for( ItF.GoToBegin(), ItC.GoToBegin(); !ItF.IsAtEnd(); ++ItF, ++ItC )
    {
    if( ItC.Get() != 1 )
      {
      ItF.Set( 0 );
      }
    if( ItF.Get() == 1 )
      {
      laplacian->SetPixel( ItF.GetIndex(), 10000 );
      }
    if( ItC.Get() == 3 )
      {
      laplacian->SetPixel( ItF.GetIndex(), 0 );
      }
    }



  typedef itk::NeighborhoodIterator<ImageType> NeighborhoodIteratorType;
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );

  unsigned int count = 0;

  float totalFieldEnergy = itk::NumericTraits<float>::Zero;
  float totalFieldEnergyNew = itk::NumericTraits<float>::max();

  bool isConverged = false;

  while( count++ < 300 && !isConverged )
    {
    totalFieldEnergyNew = 0.0;

    NeighborhoodIteratorType ItN( radius, laplacian, laplacian->GetRequestedRegion() );

    for( ItN.GoToBegin(); !ItN.IsAtEnd(); ++ItN )
      {
      if( labelReader->GetOutput()->GetPixel( ItN.GetIndex() ) == 1 )
        {
        PixelType localFieldEnergy = 0;
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          localFieldEnergy += vnl_math_sqr( 0.5 * ( ItN.GetNext( d ) - ItN.GetPrevious( d ) ) );
          }
        totalFieldEnergyNew += vcl_sqrt( localFieldEnergy );
        }

      if( contours2->GetOutput()->GetPixel( ItN.GetIndex() ) == 1 )
        {
        ItN.SetCenterPixel( 10000 );
        }
      else if( contours->GetOutput()->GetPixel( ItN.GetIndex() ) == 3 )
        {
        ItN.SetCenterPixel( 0 );
        }
      else if( labelReader->GetOutput()->GetPixel( ItN.GetIndex() ) == 1 )
        {
        PixelType localLaplacian = 0;
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          localLaplacian += ( ItN.GetNext( d ) + ItN.GetPrevious( d ) );
          }
        localLaplacian /= 6.0;
        ItN.SetCenterPixel( localLaplacian );
        }
      }
    float convergenceValue = vnl_math_abs( ( totalFieldEnergy - totalFieldEnergyNew ) / totalFieldEnergy );

    isConverged = ( convergenceValue < 1e-10 && count > 5 );
    totalFieldEnergy = totalFieldEnergyNew;
    std::cout << "Iteration " << count << ": energy = " << totalFieldEnergy << " (convergence value = " << convergenceValue << ")" << std::endl;
    }


  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[4] );
  writer->SetInput( laplacian );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << argv[0] << " imageDimension inputImage labeledMask outputImage" << std::endl;
    exit( 0 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     LaplacianLayers<2>( argc, argv );
     break;
   case 3:
     LaplacianLayers<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


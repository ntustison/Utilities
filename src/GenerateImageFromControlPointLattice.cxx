#include "itkBSplineControlPointImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageRegionIterator.h"
#include "itkTxtTransformIO.h"
#include "itkVector.h"

template <unsigned int ImageDimension>
int GenerateImageFromControlPointLattice( int argc, char *argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::Vector<PixelType, 1> ScalarType;
  typedef itk::Image<ScalarType, ImageDimension> ScalarImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typename ScalarImageType::Pointer image = ScalarImageType::New();
  image->SetOrigin( reader->GetOutput()->GetOrigin() );
  image->SetSpacing( reader->GetOutput()->GetSpacing() );
  image->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  image->SetDirection( reader->GetOutput()->GetDirection() );
  image->Allocate();
  ScalarType S;
  S.Fill( 0 );
  image->FillBuffer( S );

  itk::ImageRegionIterator<ImageType> ItR( reader->GetOutput(),
    reader->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ScalarImageType> ItS( image,
    image->GetLargestPossibleRegion() );

  for ( ItR.GoToBegin(), ItS.GoToBegin(); !ItR.IsAtEnd(); ++ItR, ++ItS )
    {
    ScalarType S;
    S.Fill( ItR.Get() );
    ItS.Set( S );
    }

  typedef itk::BSplineControlPointImageFilter
    <ScalarImageType, ScalarImageType> BSplineControlPointsFilterType;
  typename BSplineControlPointsFilterType::Pointer bspliner
    = BSplineControlPointsFilterType::New();

  typename BSplineControlPointsFilterType::ArrayType close;
  close.Fill( false );

  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  typename ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName( argv[4] );
  imageReader->Update();

  bspliner->SetSplineOrder( atoi( argv[5] ) );
  bspliner->SetCloseDimension( close );
  bspliner->SetInput( image );
  bspliner->SetOrigin( imageReader->GetOutput()->GetOrigin() );
  bspliner->SetSize( imageReader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  bspliner->SetSpacing( imageReader->GetOutput()->GetSpacing() );
  bspliner->SetDirection( imageReader->GetOutput()->GetDirection() );
  bspliner->Update();

  typename ImageType::Pointer output = ImageType::New();
  output->CopyInformation( imageReader->GetOutput() );
  output->SetRegions( imageReader->GetOutput()->GetLargestPossibleRegion() );
  output->Allocate();

  itk::ImageRegionIterator<ImageType> ItO( output,
    output->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ScalarImageType> ItB( bspliner->GetOutput(),
    bspliner->GetOutput()->GetLargestPossibleRegion() );

  for ( ItO.GoToBegin(), ItB.GoToBegin(); !ItO.IsAtEnd(); ++ItB, ++ItO )
    {
    ItO.Set( ItB.Get()[0] );
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( output );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 6 )
    {
    std::cout << "Usage: " << argv[0]
      << " imageDimension controlPointLattice outputImage referenceImage"
      << " bsplineOrder" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     GenerateImageFromControlPointLattice<2>( argc, argv );
     break;
   case 3:
     GenerateImageFromControlPointLattice<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


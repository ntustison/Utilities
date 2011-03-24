#include "itkBSplineControlPointImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageRegionIterator.h"
#include "itkVector.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc != 5 )
    {
    std::cout << "Usage: " << argv[0] << " controlPointLattice outputImage domainImage"
              << " bsplineOrder" << std::endl;
    exit( 1 );
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;   

  typedef itk::Vector<PixelType, 1> ScalarType;
  typedef itk::Image<ScalarType, ImageDimension> ScalarImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  ScalarImageType::Pointer image = ScalarImageType::New();
  image->SetOrigin( reader->GetOutput()->GetOrigin() );
  image->SetSpacing( reader->GetOutput()->GetSpacing() );
  image->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
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
  BSplineControlPointsFilterType::Pointer bspliner 
    = BSplineControlPointsFilterType::New();

  BSplineControlPointsFilterType::ArrayType close;
  close.Fill( false );

  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName( argv[3] );
  imageReader->Update();

  bspliner->SetSplineOrder( atoi( argv[4] ) );
  bspliner->SetCloseDimension( close );
  bspliner->SetInput( image );
  bspliner->SetOrigin( imageReader->GetOutput()->GetOrigin() );
  bspliner->SetSize( imageReader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  bspliner->SetSpacing( imageReader->GetOutput()->GetSpacing() );
  bspliner->Update();

  ImageType::Pointer output = ImageType::New();
  output->SetOrigin( imageReader->GetOutput()->GetOrigin() );
  output->SetSpacing( imageReader->GetOutput()->GetSpacing() );
  output->SetRegions( imageReader->GetOutput()->GetLargestPossibleRegion().GetSize() );
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
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( output );
  writer->Update();
  
  return 0;
}

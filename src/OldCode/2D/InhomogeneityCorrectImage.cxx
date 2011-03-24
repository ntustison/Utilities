#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkPointSet.h"
#include "itkStatisticsImageFilter.h"
#include "itkVector.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " inputImage outputImage [maskImage] " << std::endl;
    exit( 1 );
    }
  typedef float RealType;  

  typedef itk::Image<RealType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();
 
  ImageType::Pointer maskImage = ImageType::New();
  if ( argc > 3 )
    {
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer maskReader = ReaderType::New();
    maskReader->SetFileName( argv[3] );
    maskReader->Update();
    maskImage = maskReader->GetOutput();
    }
  else
    {
    maskImage->SetOrigin( reader->GetOutput()->GetOrigin() );
    maskImage->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
    maskImage->SetSpacing( reader->GetOutput()->GetSpacing() );
    maskImage->Allocate();
    maskImage->FillBuffer( 1.0 );
    }   

  typedef itk::Vector<RealType, 1> ScalarType;
  typedef itk::Image<ScalarType, ImageDimension> ScalarFieldType;
  typedef itk::PointSet<ScalarType, ImageDimension> PointSetType;
  PointSetType::Pointer fieldPoints = PointSetType::New();    

  itk::ImageRegionIteratorWithIndex<ImageType> It( reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<ImageType> ItM( maskImage, maskImage->GetLargestPossibleRegion() );

  unsigned long Npoints = 0;
  for ( It.GoToBegin(), ItM.GoToBegin(); !It.IsAtEnd(); ++It, ++ItM )
    {
    if ( ItM.Get() > 0 )
      {
      ScalarType scalar;
      scalar[0] = It.Get();
 
      PointSetType::PointType point;
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        point[i] = It.GetIndex()[i];
        } 
      fieldPoints->SetPoint( Npoints, point );
      fieldPoints->SetPointData( Npoints, scalar );
      Npoints++;
      }
    }   

  typedef itk::BSplineScatteredDataPointSetToImageFilter <PointSetType, ScalarFieldType> BSplineFilterType;
  BSplineFilterType::Pointer bspliner = BSplineFilterType::New();

  BSplineFilterType::ArrayType ncps;
  ncps.Fill( 64 );

//  bspliner->DebugOn();
  bspliner->SetInput( fieldPoints );
  bspliner->SetOrigin( reader->GetOutput()->GetOrigin() );
  bspliner->SetSpacing( reader->GetOutput()->GetSpacing() );
  bspliner->SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  bspliner->SetGenerateOutputImage( true );
  bspliner->SetNumberOfLevels( 1 );
  bspliner->SetSplineOrder( 3 );
  bspliner->SetNumberOfControlPoints( ncps );
  bspliner->SetInput( fieldPoints );
  bspliner->Update();

  ImageType::Pointer output = ImageType::New();
  output->SetOrigin( bspliner->GetOutput()->GetOrigin() );
  output->SetSpacing( bspliner->GetOutput()->GetSpacing() );
  output->SetRegions( bspliner->GetOutput()->GetLargestPossibleRegion() );
  output->Allocate();

  itk::ImageRegionIterator<ScalarFieldType> ItB( bspliner->GetOutput(), bspliner->GetOutput()->GetLargestPossibleRegion() );  
  itk::ImageRegionIterator<ImageType> ItI( reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );  
  itk::ImageRegionIterator<ImageType> ItO( output, output->GetLargestPossibleRegion() );  

  ItB.GoToBegin();
  ItO.GoToBegin();
  while ( !ItB.IsAtEnd() )
    {
    ItO.Set( ItB.Get()[0] ); 
    ++ItB;
    ++ItO;
    } 

  typedef itk::ImageFileWriter<ImageType> WriterType;
/*
  WriterType::Pointer bsplinewriter = WriterType::New();
  bsplinewriter->SetFileName( "bsplineField.nii" );
  bsplinewriter->SetInput( output );
  bsplinewriter->Update();
*/
  typedef itk::StatisticsImageFilter<ImageType> StatsType;
  StatsType::Pointer stats = StatsType::New();
  stats->SetInput( output );
  stats->Update();

  ItI.GoToBegin();
  ItO.GoToBegin();
  while ( !ItI.IsAtEnd() )
    {
    if ( ItO.Get() != 0 )
      {
      ItO.Set( stats->GetMean() * ItI.Get() / ItO.Get() ); 
      } 
    ++ItI;
    ++ItO;
    } 


  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( output );
  writer->Update();



  return 0;
}

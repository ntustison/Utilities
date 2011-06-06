#include <stdio.h>

#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkLaplacianRecursiveGaussianImageFilter.h"
#include "itkZeroCrossingImageFilter.h"

#include <string>

template <int ImageDimension>
int ZeroCrossingImage( int argc, char *argv[] )
{
  typedef float PixelType;
 
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<int, ImageDimension> LabelImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  int type = atoi( argv[4] );
  if( type < 1 || type > 2 )
    {
    std::cerr << "Invalid type." << std::endl; 
    return EXIT_FAILURE; 
    }  

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typename ImageType::Pointer image = ImageType::New();

  if( type == 1 )
    {
    int direction = atoi( argv[5] );
    if ( direction < 0 || direction > ImageDimension-1 )
      {
      std::cerr << "Invalid direction." << std::endl;
      return EXIT_FAILURE; 
      }
     
    typedef itk::GradientRecursiveGaussianImageFilter<ImageType>
      GradientType;
    typename GradientType::Pointer gradient = GradientType::New();
    gradient->SetInput( reader->GetOutput() );
    gradient->SetUseImageDirection( false );
    gradient->SetSigma( atof( argv[6] ) );
    gradient->SetNormalizeAcrossScale( true );
    gradient->Update();

    image->SetOrigin( reader->GetOutput()->GetOrigin() );  
    image->SetSpacing( reader->GetOutput()->GetSpacing() );  
    image->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );  
    image->Allocate();
    image->FillBuffer( 0 );
    
    itk::ImageRegionIterator<typename GradientType::OutputImageType>
      ItG( gradient->GetOutput(), gradient->GetOutput()->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<ImageType> 
      ItI( image, image->GetLargestPossibleRegion() );
    for( ItG.GoToBegin(), ItI.GoToBegin(); !ItG.IsAtEnd(); ++ItG, ++ItI )
      {
      ItI.Set( ItG.Get()[direction] );
      }

    typedef itk::ZeroCrossingImageFilter<ImageType, LabelImageType> ZeroCrossingsType;
    typename ZeroCrossingsType::Pointer zerocrossings = ZeroCrossingsType::New();
    zerocrossings->SetInput( image );
    zerocrossings->SetBackgroundValue( 0 );
    zerocrossings->SetForegroundValue( 1 );
    zerocrossings->Update();
  
    itk::ImageRegionIterator<LabelImageType> ItZ( zerocrossings->GetOutput(),
      zerocrossings->GetOutput()->GetLargestPossibleRegion() );
    typedef itk::NeighborhoodIterator<ImageType> NeighborhoodIteratorType; 
    typename NeighborhoodIteratorType::RadiusType radius;
    radius.Fill( 1 );
    NeighborhoodIteratorType ItN( radius, image, 
      image->GetLargestPossibleRegion() );
    
    for( ItN.GoToBegin(), ItZ.GoToBegin(); !ItN.IsAtEnd(); ++ItN, ++ItZ )
      {
      if( ItZ.Get() == zerocrossings->GetForegroundValue() )
        {
        if( ItN.GetCenterPixel() > ItN.GetNext( direction ) &&
            ItN.GetCenterPixel() < ItN.GetPrevious( direction ) )
          {
          ItZ.Set( zerocrossings->GetForegroundValue() ); 
          }
        else
          {
          ItZ.Set( zerocrossings->GetBackgroundValue() ); 
          }  
        } 
      }
    
  
    typedef itk::ImageFileWriter<LabelImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( zerocrossings->GetOutput() );
    writer->SetFileName( argv[3] );                                          
    writer->Update();

    }
  else
    {
    typedef itk::LaplacianRecursiveGaussianImageFilter<ImageType, ImageType>
      LaplacianType;
    typename LaplacianType::Pointer laplacian = LaplacianType::New();
    laplacian->SetInput( reader->GetOutput() );
    laplacian->SetSigma( atof( argv[5] ) );
    laplacian->SetNormalizeAcrossScale( true );
    laplacian->Update();
    
    image = laplacian->GetOutput();

    typedef itk::ZeroCrossingImageFilter<ImageType, LabelImageType> ZeroCrossingsType;
    typename ZeroCrossingsType::Pointer zerocrossings = ZeroCrossingsType::New();
    zerocrossings->SetInput( image );
    zerocrossings->SetBackgroundValue( 0 );
    zerocrossings->SetForegroundValue( 1 );
    zerocrossings->Update();
  
    typedef itk::ImageFileWriter<LabelImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( zerocrossings->GetOutput() );
    writer->SetFileName( argv[3] );                                          
    writer->Update();

    }    

    
  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension inputImage outputImage " 
      << "type parameters" << std::endl;
    std::cerr << "   type, parameters" << std::endl;
    std::cerr << "     1. firstDerivative, direction sigma" << std::endl;
    std::cerr << "     2. laplacian, sigma" << std::endl;
    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) ) 
   {
   case 2:
     ZeroCrossingImage<2>( argc, argv );
     break;
   case 3:
     ZeroCrossingImage<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


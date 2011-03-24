#include "itkArray.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileWriter.h"

#include <string>
#include <vector>

#define SCHEME 1

template <unsigned int ImageDimension>
int CreateADCImage( int argc, char * argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType>  ReaderType;

  typedef unsigned int LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;
  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;

  /**
   * list the files
   */

  std::vector<typename ImageType::Pointer> images;
  std::vector<float> bvalues;

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[3] );
  reader->Update();
  images.push_back( reader->GetOutput() );
  bvalues.push_back( 0 );

  for( unsigned int n = 4; n < static_cast<unsigned int>( argc ) - 1 ; n+=2 )
    {
    bvalues.push_back( atof( argv[n] ) );
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[n+1] );
    reader->Update();
    images.push_back( reader->GetOutput() );
    }

  typename ImageType::Pointer output = ImageType::New();
  output->SetOrigin( images[0]->GetOrigin() );
  output->SetSpacing( images[0]->GetSpacing() );
  output->SetRegions( images[0]->GetLargestPossibleRegion() );
  output->SetDirection( images[0]->GetDirection() );
  output->Allocate();
  output->FillBuffer( 0 );

  /**
   * Get the actual labels---assume that the first image read has all
   * the actual labels.
   */
  itk::ImageRegionIteratorWithIndex<ImageType> It( images[0],
    images[0]->GetLargestPossibleRegion() );


  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    float D = 1.0;

    float So = It.Get();
    if( images.size() == 2 )
      {
      float S = images[1]->GetPixel( It.GetIndex() );
      D = vcl_log( S / So ) / ( -bvalues[1] );
      }
    else
      {
      float sumLnY = 0.0;
      float sumXLnY = 0.0;
      float sumX = 0.0;
      float sumX2 = 0.0;

      float sumY = 0.0;
      float sumXY = 0.0;
      float sumXYLnY = 0.0;
      float sumYLnY = 0.0;
      float sumX2Y = 0.0;

      for( unsigned int n = 1; n < images.size(); n++ )
        {
        float S = images[n]->GetPixel( It.GetIndex() );

        // scheme 1  This fit gives greater weights to small values so,
        //  in order to weight the points equally, it is often better to
        //  use scheme 2
        //  http://mathworld.wolfram.com/LeastSquaresFittingExponential.html
        sumLnY += vcl_log( S / So );
        sumXLnY += ( bvalues[n] * vcl_log( S / So ) );
        sumX += bvalues[n];
        sumX2 += vnl_math_sqr( bvalues[n] );

        // scheme 2
        sumY += ( S / So );
        sumXY += ( ( S / So ) * bvalues[n] );
        sumXYLnY += ( ( S / So ) * bvalues[n] ) * vcl_log( S / So );
        sumXYLnY += ( S / So ) * vcl_log( S / So );
        sumX2Y += ( ( S / So ) * vnl_math_sqr( bvalues[n] ) );
        }
      float n = static_cast<float>( images.size() - 1 );

      if( SCHEME == 1 )
        {
        D = ( n * sumXLnY - sumX * sumLnY ) /
          ( n * sumX2 - vnl_math_sqr( sumX ) );
        }
      else
        {
        // scheme 2
        D = ( sumY * sumXYLnY - sumXY * sumYLnY ) /
          ( sumX2Y - vnl_math_sqr( sumXY ) );
        }
      }
    output->SetPixel( It.GetIndex(), vnl_math_max( 0.0f, D ) );
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( output );
  writer->SetFileName( argv[2] );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if( argc < 6 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " imageDimension outputImage B0_image "
      << "B1value B1image B2value B2image ... Bnvalue Bnimage" << std::endl;
    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     CreateADCImage<2>( argc, argv );
     break;
   case 3:
     CreateADCImage<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}




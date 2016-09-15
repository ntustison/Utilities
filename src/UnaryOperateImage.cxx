#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileWriter.h"
#include "itkGaussianInterpolateImageFunction.h"
#include "itkNumericSeriesFileNames.h"

#include "Common.h"

template <unsigned int ImageDimension>
int UnaryOperateImage( int argc, char * argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType>  ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  if( argv[3][0] == 'p' )
    {
    PixelType constant = static_cast<PixelType>( atof( argv[4] ) );

    for( unsigned int n = 6; n < static_cast<unsigned int>( argc ); n++ )
      {
      typename ImageType::IndexType index;

      std::vector<int> idx = ConvertVector<int>( std::string( argv[n] ) );
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        index[d] = idx[d];
        }
      reader->GetOutput()->SetPixel( index, constant );
      }
    }
  else if( argv[3][0] == 'i' )
    {
    typename ImageType::Pointer inputImage = reader->GetOutput();
    inputImage->DisconnectPipeline();

    itk::NumericSeriesFileNames::Pointer fileNamesCreator = itk::NumericSeriesFileNames::New();
    fileNamesCreator->SetStartIndex( 0 );
    fileNamesCreator->SetEndIndex( ImageDimension );
    fileNamesCreator->SetSeriesFormat( argv[5] );
    const std::vector<std::string> & imageNames = fileNamesCreator->GetFileNames();

    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      typename ImageType::Pointer indexImage = ImageType::New();
      indexImage->CopyInformation( inputImage );
      indexImage->SetRegions( inputImage->GetRequestedRegion() );
      indexImage->Allocate();
      indexImage->FillBuffer( 0 );

      itk::ImageRegionIteratorWithIndex<ImageType> It( inputImage,
        inputImage->GetRequestedRegion() );
      for( It.GoToBegin(); !It.IsAtEnd(); ++It )
        {
        typename ImageType::IndexType index = It.GetIndex();
        indexImage->SetPixel( index, static_cast<PixelType>( index[d] ) );
        }
      typedef itk::ImageFileWriter<ImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( indexImage );
      writer->SetFileName( imageNames[d] );
      writer->Update();
      }
    return EXIT_SUCCESS;
    }
  else if( argv[3][0] == 'q' )
    {
    PixelType constant = static_cast<PixelType>( atof( argv[4] ) );

    for( unsigned int n = 6; n < static_cast<unsigned int>( argc ); n++ )
      {
      typename ImageType::IndexType index;
      typename ImageType::PointType point;

      std::vector<float> pt = ConvertVector<float>( std::string( argv[n] ) );
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        point[d] = pt[d];
        }
      reader->GetOutput()->TransformPhysicalPointToIndex( point, index );

      reader->GetOutput()->SetPixel( index, constant );
      }
    }
  else if( argv[3][0] == 'g' )
    {
    typedef itk::GaussianInterpolateImageFunction<ImageType, double>
      GaussianInterpolatorType;
    typename GaussianInterpolatorType::Pointer g_interpolator
      = GaussianInterpolatorType::New();
    g_interpolator->SetInputImage( reader->GetOutput() );

    double sigma[ImageDimension];
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      sigma[d] = reader->GetOutput()->GetSpacing()[d];
      }
    double alpha = 1.0;
    g_interpolator->SetParameters( sigma, alpha );

    typename GaussianInterpolatorType::PointType point;
    point.Fill( 0.0 );

    std::vector<double> pt = ConvertVector<double>( std::string( argv[4] ) );
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      point[d] = pt[d];
      }

    if( g_interpolator->IsInsideBuffer( point ) )
      {
      std::cout << g_interpolator->Evaluate( point ) << std::endl;
      }
    else
      {
      std::cout << "outside image buffer" << std::endl;
      }
    }
  else
    {
    itk::ImageRegionIteratorWithIndex<ImageType> It( reader->GetOutput(),
      reader->GetOutput()->GetLargestPossibleRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      switch( argv[3][0] )
        {
        case '+':
          {
          It.Set( It.Get() + atof( argv[4] ) );
          break;
          }
        case '-':
          {
          It.Set( It.Get() - atof( argv[4] ) );
          break;
          }
        case 'x':
          {
          It.Set( It.Get() * atof( argv[4] ) );
          break;
          }
        case '/':
          {
          It.Set( It.Get() / atof( argv[4] ) );
          break;
          }
        case '^':
          {
          It.Set( static_cast<PixelType>( vcl_pow( static_cast<double>( It.Get() ),
            static_cast<double>( atof( argv[4] ) ) ) ) );
          break;
          }
        case 'f':
          {
          It.Set( vnl_math_max( static_cast<float>( It.Get() ), static_cast<float>( atof( argv[4] ) ) ) );
          break;
          }
        case 'e':
          {
          It.Set( static_cast<PixelType>( vcl_exp( static_cast<double>( It.Get() ) ) ) );
          break;
          }
        case 'l':
          {
          It.Set( static_cast<PixelType>( vcl_log( static_cast<double>( It.Get() ) ) ) );
          break;
          }
        case 'b':
          {
          It.Set( 1.0 / ( 1.0 + It.Get() ) );
          break;
          }
        case 's':
          {
          It.Set( 1.0 / ( 1.0 + vcl_exp( -( It.Get() - atof( argv[7] ) ) / ( atof( argv[6] ) ) ) ) );
          break;
          }
        case 'r':
          {
          if( strcmp( "nan", argv[6] ) == 0 )
            {
            if( vnl_math_isnan( It.Get() ) )
              {
              It.Set( atof( argv[7] ) );
              }
            }
          else if( strcmp( "inf", argv[6] ) == 0 )
            {
            if( vnl_math_isinf( It.Get() ) )
              {
              It.Set( atof( argv[7] ) );
              }
            }
          else if( It.Get() == atof( argv[6] ) )
            {
            It.Set( atof( argv[7] ) );
            }
          break;
          }
        case 't':
          {
          if( It.Get() >= atof( argv[6] ) && It.Get() <= atof( argv[7] ) )
            {
            It.Set( atof( argv[8] ) );
            }
          break;
          }
        default:
          {
          std::cerr << "Error: Unknown operation." << std::endl;
          exit( 1 );
          break;
          }
        }
      }
    }

  if( argc > 5 )
    {
    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( reader->GetOutput() );
    writer->SetFileName( argv[5] );
    writer->Update();
    }

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if( argc < 5 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " imageDimension inputImage "
              << "operation constant [outputImage]" << std::endl;
    std::cerr << "  operation: " << std::endl;
    std::cerr << "    +:   Add" << std::endl;
    std::cerr << "    -:   Subtract" << std::endl;
    std::cerr << "    x:   Multiply" << std::endl;
    std::cerr << "    /:   Divide" << std::endl;
    std::cerr << "    ^:   pow" << std::endl;
    std::cerr << "    f:   set voxel to floor value, i.e. max( voxel, constant )" <<  std::endl;
    std::cerr << "    p:   set pixel to constant value [index1] [index2] ... index[n]" <<  std::endl;
    std::cerr << "    q:   set pixel at physical point to constant value [point1] [point2] ... point[n]" <<  std::endl;
    std::cerr << "    g:   get pixel value at physical point (gaussian interpolation)" << std::endl;
    std::cerr << "  The following operations ignore the \'constant\' argument." << std::endl;
    std::cerr << "    i:   index images (outputImage is specified in c-style format, i.e. foo%d.nii.gz)" << std::endl;
    std::cerr << "    e:   exp" << std::endl;
    std::cerr << "    l:   ln" << std::endl;
    std::cerr << "    b:   bounded reciprocal" << std::endl;
    std::cerr << "    s:   sigmoid function [alpha] [beta]" << std::endl;
    std::cerr << "    r:   replace [oldPixel] [newPixel]" << std::endl;
    std::cerr << "    t:   threshold/replace [lowPixel] [highPixel] [newPixel]" << std::endl;
    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     UnaryOperateImage<2>( argc, argv );
     break;
   case 3:
     UnaryOperateImage<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}




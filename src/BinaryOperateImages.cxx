#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkImageFileWriter.h"

#include <string>
#include <vector>
#include <sstream>

template<class TValue>
TValue Convert( std::string optionString )
{
  TValue value;
  std::istringstream iss( optionString );
  iss >> value;
  return value;
}

template<class TValue>
std::vector<TValue> ConvertVector( std::string optionString )
{
  std::vector<TValue> values;
  std::string::size_type crosspos = optionString.find( 'x', 0 );

  if ( crosspos == std::string::npos )
    {
    values.push_back( Convert<TValue>( optionString ) );
    }
  else
    {
    std::string element = optionString.substr( 0, crosspos ) ;
    TValue value;
    std::istringstream iss( element );
    iss >> value;
    values.push_back( value );
    while ( crosspos != std::string::npos )
      {
      std::string::size_type crossposfrom = crosspos;
      crosspos = optionString.find( 'x', crossposfrom + 1 );
      if ( crosspos == std::string::npos )
        {
        element = optionString.substr( crossposfrom + 1, optionString.length() );
        }
      else
        {
        element = optionString.substr( crossposfrom + 1, crosspos ) ;
        }
      std::istringstream iss2( element );
      iss2 >> value;
      values.push_back( value );
      }
    }
  return values;
}

typedef float RealType;

template <unsigned int ImageDimension>
int BinaryOperateImages( int argc, char * argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType>  ReaderType;

  typename ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[2] );
  reader1->Update();

  typename ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[4] );
  reader2->Update();

  typename ReaderType::Pointer reader3 = ReaderType::New();
  if( argc > 6 )
    {
    reader3->SetFileName( argv[6] );
    reader3->Update();
    }

  typename ImageType::Pointer output = ImageType::New();
  output->SetOrigin( reader1->GetOutput()->GetOrigin() );
  output->SetSpacing( reader1->GetOutput()->GetSpacing() );
  output->SetRegions( reader1->GetOutput()->GetLargestPossibleRegion() );
  output->SetDirection( reader1->GetOutput()->GetDirection() );
  output->Allocate();
  output->FillBuffer( 0 );

  std::string op = std::string( argv[3] );

  typename ImageType::SizeType radius;
  radius.Fill( 1 );

  PixelType delta = -1.0;

  itk::NeighborhoodIterator<ImageType> It( radius, output,
    output->GetLargestPossibleRegion() );
  itk::NeighborhoodIterator<ImageType> It1( radius, reader1->GetOutput(),
    reader1->GetOutput()->GetLargestPossibleRegion() );
  itk::NeighborhoodIterator<ImageType> It2( radius, reader2->GetOutput(),
    reader2->GetOutput()->GetLargestPossibleRegion() );

  for( It.GoToBegin(), It1.GoToBegin(), It2.GoToBegin();
    !It.IsAtEnd(); ++It, ++It1, ++It2 )
    {
    if( op.compare( "+" ) == 0 )
      {
      It.SetCenterPixel( It1.GetCenterPixel() + It2.GetCenterPixel() );
      }
    else if( op.compare( "-" ) == 0 )
      {
      It.SetCenterPixel( It1.GetCenterPixel() - It2.GetCenterPixel() );
      }
    else if( op.compare( "x" ) == 0 )
      {
      It.SetCenterPixel( It1.GetCenterPixel() * It2.GetCenterPixel() );
      }
    else if( op.compare( "/" ) == 0 )
      {
      It.SetCenterPixel( It1.GetCenterPixel() / It2.GetCenterPixel() );
      }
    else if( op.compare( "or" ) == 0 )
      {
      if( static_cast<bool>( It1.GetCenterPixel() ) ||
        static_cast<bool>( It2.GetCenterPixel() ) )
        {
        It.SetCenterPixel( itk::NumericTraits<PixelType>::One );
        }
      else
        {
        It.SetCenterPixel( itk::NumericTraits<PixelType>::Zero );
        }
      }
    else if( op.compare( "xor" ) == 0 )
      {
      if( ( static_cast<bool>( It1.GetCenterPixel() ) ||
        static_cast<bool>( It2.GetCenterPixel() ) ) &&
        ( static_cast<bool>( It1.GetCenterPixel() ) !=
        static_cast<bool>( It2.GetCenterPixel() ) ) )
        {
        It.SetCenterPixel( itk::NumericTraits<PixelType>::One );
        }
      else
        {
        It.SetCenterPixel( itk::NumericTraits<PixelType>::Zero );
        }
      }
    else if( op.compare( "and" ) == 0 )
      {
      if( static_cast<bool>( It1.GetCenterPixel() ) &&
        static_cast<bool>( It2.GetCenterPixel() ) )
        {
        It.SetCenterPixel( itk::NumericTraits<PixelType>::One );
        }
      else
        {
        It.SetCenterPixel( itk::NumericTraits<PixelType>::Zero );
        }
      }
    else if( op.compare( "max" ) == 0 )
      {
      It.SetCenterPixel( vnl_math_max( It1.GetCenterPixel(), It2.GetCenterPixel() ) );
      }
    else if( op.compare( "replace" ) == 0 )
      {
      if( It2.GetCenterPixel() != 0 )
        {
        It1.SetCenterPixel( It2.GetCenterPixel() );
        }
      }
    else if( op.compare( "isgreaterthan" ) == 0 )
      {
      It.SetCenterPixel( ( It1.GetCenterPixel() > It2.GetCenterPixel() ) ? 1 : 0 );
      }
    else if( op.compare( "islessthan" ) == 0 )
      {
      It.SetCenterPixel( ( It1.GetCenterPixel() < It2.GetCenterPixel() ) ? 1 : 0 );
      }
    else if( op.compare( "min" ) == 0 )
      {
      It.SetCenterPixel( vnl_math_min( It1.GetCenterPixel(), It2.GetCenterPixel() ) );
      }
    else if( op.compare( "zscore" ) == 0 )
      {
      if( reader3->GetOutput() == NULL )
        {
        std::cerr << "Need to specify third image." << std::endl;
        return EXIT_FAILURE;
        }

      RealType mean = It1.GetCenterPixel();
      RealType std = vcl_sqrt( It2.GetCenterPixel() );
      RealType pixel = reader3->GetOutput()->GetPixel( It1.GetIndex() );
      if( std > 0 )
        {
        It.SetCenterPixel( ( pixel - mean ) / std );
        }
      }
    else if( op.compare( "misc" ) == 0 )
      {
      if( It1.GetCenterPixel() == 0 )
        {
        continue;
        }

      bool isSurroundedIn2 = true;
      bool isIsolated = true;

      for( unsigned int n = 0; n < It1.GetNeighborhood().Size(); n++ )
        {
        isSurroundedIn2 = isSurroundedIn2 && static_cast<bool>( It2.GetPixel( n ) );
        if( n == static_cast<unsigned int>( 0.5 * It1.Size() ) )
          {
          continue;
          }
        if( It1.GetPixel( n ) != 0 )
          {
          isIsolated = false;
          break;
          }
        }
      if( isIsolated /*&& isSurroundedIn2*/ )
        {
        It.SetCenterPixel( It1.GetCenterPixel() );
        }
      }
    else
      {
      std::cerr << "Error: Unknown operation." << std::endl;
      return EXIT_FAILURE;
      }
    }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( output );
  writer->SetFileName( argv[5] );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if( argc < 6 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " imageDimension inputImage1 operation "
              << " inputImage2 outputImage [inputImage3]" << std::endl;
    std::cerr << "  operation: " << std::endl;
    std::cerr << "    +:   Add" << std::endl;
    std::cerr << "    -:   Subtract" << std::endl;
    std::cerr << "    x:   Multiply" << std::endl;
    std::cerr << "    /:   Divide" << std::endl;
    std::cerr << "    replace:  replace pixels in image1 with the non-zero pixels in image2" << std::endl;
    std::cerr << "    min: voxel-wise minimum" << std::endl;
    std::cerr << "    max: voxel-wise maximum" << std::endl;
    std::cerr << "    zscore:  inputImage1 = meanImage, inputImage2 = varianceImage, z = (pixel - mean)/std" << std::endl;
    std::cerr << "    isgreaterthan: mask (1 if image 1 is greather than image 2, 0 o.w.)" << std::endl;
    std::cerr << "    islessthan: mask (1 if image 1 is greather than image 2, 0 o.w.)" << std::endl;
    std::cerr << "    or:  logical \'or\'" << std::endl;
    std::cerr << "    xor: logical \'xor\'" << std::endl;
    std::cerr << "    and: logical \'and\'" << std::endl;
    std::cerr << "    misc: programmer defined" << std::endl;
    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     BinaryOperateImages<2>( argc, argv );
     break;
   case 3:
     BinaryOperateImages<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}




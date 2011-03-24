#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkImageFileWriter.h"

#include <string>

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
    else if( op.compare( "min" ) == 0 )
      {
      It.SetCenterPixel( vnl_math_min( It1.GetCenterPixel(), It2.GetCenterPixel() ) );
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
              << " inputImage2 outputImage " << std::endl;
    std::cerr << "  operation: " << std::endl;
    std::cerr << "    +:   Add" << std::endl;
    std::cerr << "    -:   Subtract" << std::endl;
    std::cerr << "    x:   Multiply" << std::endl;
    std::cerr << "    /:   Divide" << std::endl;
    std::cerr << "    min: voxel-wise minimum" << std::endl;
    std::cerr << "    max: voxel-wise maximum'" << std::endl;
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




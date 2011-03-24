#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h" 
#include "itkBinaryBoxStructuringElement.h" 
#include "itkBinaryDiamondStructuringElement.h" 

#include "vnl/vnl_math.h"

template <unsigned int ImageDimension>
int ErodeImage( int argc, char * argv[] )
{

  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType>  ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );

  PixelType label = itk::NumericTraits<PixelType>::One;
  if ( argc > 6 )
    {
    label = static_cast<PixelType>( atof( argv[6] ) );
    } 

  if ( argc < 5 || atoi( argv[5] ) == 1 )
    {
    typedef itk::BinaryBallStructuringElement< 
                        PixelType,
                        ImageDimension>             StructuringElementType;
  
    typedef itk::BinaryErodeImageFilter<
                              ImageType, 
                              ImageType,
                              StructuringElementType >  ErodeFilterType;
    typename ErodeFilterType::Pointer  binaryErode = ErodeFilterType::New();
    StructuringElementType  structuringElement;
    structuringElement.SetRadius( atoi( argv[4] ) ); 
    structuringElement.CreateStructuringElement();
  
    binaryErode->SetKernel( structuringElement );
    binaryErode->SetInput( reader->GetOutput() );
    binaryErode->SetForegroundValue( label );
    binaryErode->Update();
  
    typedef itk::ImageFileWriter<ImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( binaryErode->GetOutput() );
    writer->SetFileName( argv[3] );
    writer->Update();
    }
  else if ( atoi( argv[5] ) == 0 )
    {
    typedef itk::BinaryBoxStructuringElement< 
                        PixelType,
                        ImageDimension>             StructuringElementType;
  
    typedef itk::BinaryErodeImageFilter<
                              ImageType, 
                              ImageType,
                              StructuringElementType >  ErodeFilterType;
    typename ErodeFilterType::Pointer  binaryErode = ErodeFilterType::New();
    StructuringElementType  structuringElement;
    structuringElement.SetRadius( atoi( argv[4] ) ); 
    structuringElement.CreateStructuringElement();
  
    binaryErode->SetKernel( structuringElement );
    binaryErode->SetInput( reader->GetOutput() );
    binaryErode->SetForegroundValue( label );
    binaryErode->Update();
  
    typedef itk::ImageFileWriter<ImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( binaryErode->GetOutput() );
    writer->SetFileName( argv[3] );
    writer->Update();
    }
  else
    {
    typedef itk::BinaryDiamondStructuringElement< 
                        PixelType,
                        ImageDimension>             StructuringElementType;
  
    typedef itk::BinaryErodeImageFilter<
                              ImageType, 
                              ImageType,
                              StructuringElementType >  ErodeFilterType;
    typename ErodeFilterType::Pointer  binaryErode = ErodeFilterType::New();
    StructuringElementType  structuringElement;
    structuringElement.SetRadius( atoi( argv[4] ) ); 
    structuringElement.CreateStructuringElement();
  
    binaryErode->SetKernel( structuringElement );
    binaryErode->SetInput( reader->GetOutput() );
    binaryErode->SetForegroundValue( label );
    binaryErode->Update();
  
    typedef itk::ImageFileWriter<ImageType>  WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( binaryErode->GetOutput() );
    writer->SetFileName( argv[3] );
    writer->Update();
    }
  

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if( argc < 5 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " imageDimension inputImage outputImage radius  [type: box == 0, ball = 1, diamond = 2] [label]" << std::endl;
    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) ) 
   {
   case 2:
     ErodeImage<2>( argc, argv );
     break;
   case 3:
     ErodeImage<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}



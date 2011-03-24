#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

#include "itkBinaryBoundedSpaceDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h" 
#include "itkBinaryBoxStructuringElement.h" 
#include "itkBinaryDiamondStructuringElement.h" 

#include "vnl/vnl_math.h"
#include "global.h"

int main( int argc, char * argv[] )
{
  if( argc < 5 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputImageFile boundedSpaceImage outputImageFile radius " 
              << "[type: box == 0, ball = 1, diamond = 2] [label] [boundedSpaceValue] [scaling]" << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType>  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::ImageFileReader<ImageType>  LabelReaderType;
  LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName( argv[2] );
  labelReader->Update();

  PixelType label = itk::NumericTraits<PixelType>::One;
  if ( argc > 6 )
    {
    label = static_cast<PixelType>( atof( argv[6] ) );
    } 
  PixelType bsvalue = itk::NumericTraits<PixelType>::One;
  if ( argc > 7 )
    {
    bsvalue = static_cast<PixelType>( atof( argv[7] ) );
    } 
  int scaling = 1;
  if ( argc > 8 )
    {
    scaling = atoi( argv[8] );
    } 
 

  if ( argc < 5 || atoi( argv[5] ) == 1 )
    {
    typedef itk::BinaryBallStructuringElement< 
                        PixelType,
                        ImageDimension>             StructuringElementType;
  
    typedef itk::BinaryBoundedSpaceDilateImageFilter<
                              ImageType, 
                              ImageType,
                              StructuringElementType >  DilateFilterType;
    DilateFilterType::Pointer  binaryDilate = DilateFilterType::New();
    StructuringElementType  structuringElement;
    structuringElement.SetRadius( atoi( argv[4] ) ); 
    structuringElement.CreateStructuringElement();
  
    binaryDilate->SetKernel( structuringElement );
    binaryDilate->SetInput( reader->GetOutput() );
    binaryDilate->SetScaling( scaling );
    binaryDilate->SetForegroundValue( label );
    binaryDilate->SetBoundedSpaceImage( labelReader->GetOutput() );
    binaryDilate->SetBoundedSpaceValue( bsvalue );
    binaryDilate->Update();
  
    typedef itk::ImageFileWriter<ImageType>  WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput( binaryDilate->GetOutput() );
    writer->SetFileName( argv[3] );
    writer->Update();
    }
  else if ( atoi( argv[5] ) == 2 )
    {
    typedef itk::BinaryBoxStructuringElement< 
                        PixelType,
                        ImageDimension>             StructuringElementType;
  
    typedef itk::BinaryBoundedSpaceDilateImageFilter<
                              ImageType, 
                              ImageType,
                              StructuringElementType >  DilateFilterType;
    DilateFilterType::Pointer  binaryDilate = DilateFilterType::New();
    StructuringElementType  structuringElement;
    structuringElement.SetRadius( atoi( argv[4] ) ); 
    structuringElement.CreateStructuringElement();
  
    binaryDilate->SetKernel( structuringElement );
    binaryDilate->SetInput( reader->GetOutput() );
    binaryDilate->SetScaling( scaling );
    binaryDilate->SetForegroundValue( label );
    binaryDilate->SetBoundedSpaceImage( labelReader->GetOutput() );
    binaryDilate->SetBoundedSpaceValue( bsvalue );
    binaryDilate->Update();
  
    typedef itk::ImageFileWriter<ImageType>  WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput( binaryDilate->GetOutput() );
    writer->SetFileName( argv[3] );
    writer->Update();
    }
  else 
    {
    typedef itk::BinaryDiamondStructuringElement< 
                        PixelType,
                        ImageDimension>             StructuringElementType;
  
    typedef itk::BinaryBoundedSpaceDilateImageFilter<
                              ImageType, 
                              ImageType,
                              StructuringElementType >  DilateFilterType;
    DilateFilterType::Pointer  binaryDilate = DilateFilterType::New();
    StructuringElementType  structuringElement;
    structuringElement.SetRadius( atoi( argv[4] ) ); 
    structuringElement.CreateStructuringElement();
  
    binaryDilate->SetKernel( structuringElement );
    binaryDilate->SetInput( reader->GetOutput() );
    binaryDilate->SetScaling( scaling );
    binaryDilate->SetForegroundValue( label );
    binaryDilate->SetBoundedSpaceImage( labelReader->GetOutput() );
    binaryDilate->SetBoundedSpaceValue( bsvalue );
    binaryDilate->Update();
  
    typedef itk::ImageFileWriter<ImageType>  WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput( binaryDilate->GetOutput() );
    writer->SetFileName( argv[3] );
    writer->Update();
    }
  

  return EXIT_SUCCESS;
}



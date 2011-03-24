#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

#include "itkGrayscaleDilateImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkGrayscaleMorphologicalClosingImageFilter.h"
#include "itkGrayscaleMorphologicalOpeningImageFilter.h"
#include "itkGrayscaleGrindPeakImageFilter.h"
#include "itkBinaryBallStructuringElement.h" 
#include "itkBinaryBoxStructuringElement.h" 
#include "itkBinaryDiamondStructuringElement.h" 

#include "vnl/vnl_math.h"

template <unsigned int ImageDimension>
int GrayscaleMorphology( int argc, char * argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType>  ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );

  unsigned int radius = 1;
  if ( argc > 5 )
    {
    radius = atoi( argv[5] ); 
    }  

  
  unsigned int operation = atoi( argv[4] );
  if ( operation < 4 )
    {
    if ( argc < 6 || atoi( argv[6] ) == 1 )
      {
      typedef itk::BinaryBallStructuringElement< 
                          PixelType,
                          ImageDimension> StructuringElementType;
      StructuringElementType  element;
      element.SetRadius( radius ); 
      element.CreateStructuringElement();
     
      switch ( operation )
        {
        case 0:
          {
          typedef itk::GrayscaleDilateImageFilter<ImageType, ImageType,
            StructuringElementType >  FilterType;
          typename FilterType::Pointer  filter = FilterType::New();
          filter->SetKernel( element );
          filter->SetInput( reader->GetOutput() );
          filter->Update();
          typedef itk::ImageFileWriter<ImageType>  WriterType;
          typename WriterType::Pointer writer = WriterType::New();
          writer->SetInput( filter->GetOutput() );
          writer->SetFileName( argv[3] );
          writer->Update();
          break;
          } 
        case 1:
          {
          typedef itk::GrayscaleErodeImageFilter<ImageType, ImageType,
            StructuringElementType >  FilterType;
          typename FilterType::Pointer  filter = FilterType::New();
          filter->SetKernel( element );
          filter->SetInput( reader->GetOutput() );
          filter->Update();
          typedef itk::ImageFileWriter<ImageType>  WriterType;
          typename WriterType::Pointer writer = WriterType::New();
          writer->SetInput( filter->GetOutput() );
          writer->SetFileName( argv[3] );
          writer->Update();
          break;
          } 
        case 2:
          {
          typedef itk::GrayscaleMorphologicalClosingImageFilter<ImageType, ImageType,
            StructuringElementType >  FilterType;
          typename FilterType::Pointer  filter = FilterType::New();
          filter->SetKernel( element );
          filter->SetInput( reader->GetOutput() );
          filter->Update();
          typedef itk::ImageFileWriter<ImageType>  WriterType;
          typename WriterType::Pointer writer = WriterType::New();
          writer->SetInput( filter->GetOutput() );
          writer->SetFileName( argv[3] );
          writer->Update();
          break;
          } 
        case 3:
          {
          typedef itk::GrayscaleMorphologicalOpeningImageFilter<ImageType, ImageType,
            StructuringElementType >  FilterType;
          typename FilterType::Pointer  filter = FilterType::New();
          filter->SetKernel( element );
          filter->SetInput( reader->GetOutput() );
          filter->Update();
          typedef itk::ImageFileWriter<ImageType>  WriterType;
          typename WriterType::Pointer writer = WriterType::New();
          writer->SetInput( filter->GetOutput() );
          writer->SetFileName( argv[3] );
          writer->Update();
          break;
          } 
        default:
          {
          std::cerr << "Invalid operation choice." << std::endl; 
          return EXIT_FAILURE;
          }
        }
      }
    else if ( atoi( argv[6] ) == 0 )
      {
      typedef itk::BinaryBoxStructuringElement< 
                          PixelType,
                          ImageDimension>  StructuringElementType;
      StructuringElementType element;
      element.SetRadius( radius ); 
      element.CreateStructuringElement();
      
      switch ( operation )
        {
        case 0:
          {
          typedef itk::GrayscaleDilateImageFilter<ImageType, ImageType,
            StructuringElementType >  FilterType;
          typename FilterType::Pointer  filter = FilterType::New();
          filter->SetKernel( element );
          filter->SetInput( reader->GetOutput() );
          filter->Update();
          typedef itk::ImageFileWriter<ImageType>  WriterType;
          typename WriterType::Pointer writer = WriterType::New();
          writer->SetInput( filter->GetOutput() );
          writer->SetFileName( argv[3] );
          writer->Update();
          break;
          } 
        case 1:
          {
          typedef itk::GrayscaleErodeImageFilter<ImageType, ImageType,
            StructuringElementType >  FilterType;
          typename FilterType::Pointer  filter = FilterType::New();
          filter->SetKernel( element );
          filter->SetInput( reader->GetOutput() );
          filter->Update();
          typedef itk::ImageFileWriter<ImageType>  WriterType;
          typename WriterType::Pointer writer = WriterType::New();
          writer->SetInput( filter->GetOutput() );
          writer->SetFileName( argv[3] );
          writer->Update();
          break;
          } 
        case 2:
          {
          typedef itk::GrayscaleMorphologicalClosingImageFilter<ImageType, ImageType,
            StructuringElementType >  FilterType;
          typename FilterType::Pointer  filter = FilterType::New();
          filter->SetKernel( element );
          filter->SetInput( reader->GetOutput() );
          filter->Update();
          typedef itk::ImageFileWriter<ImageType>  WriterType;
          typename WriterType::Pointer writer = WriterType::New();
          writer->SetInput( filter->GetOutput() );
          writer->SetFileName( argv[3] );
          writer->Update();
          break;
          } 
        case 3:
          {
          typedef itk::GrayscaleMorphologicalOpeningImageFilter<ImageType, ImageType,
            StructuringElementType >  FilterType;
          typename FilterType::Pointer  filter = FilterType::New();
          filter->SetKernel( element );
          filter->SetInput( reader->GetOutput() );
          filter->Update();
          typedef itk::ImageFileWriter<ImageType>  WriterType;
          typename WriterType::Pointer writer = WriterType::New();
          writer->SetInput( filter->GetOutput() );
          writer->SetFileName( argv[3] );
          writer->Update();
          break;
          } 
        default:
          {
          std::cerr << "Invalid operation choice." << std::endl; 
          return EXIT_FAILURE;
          }
        }
      }  
    else
      {
      typedef itk::BinaryDiamondStructuringElement< 
                          PixelType,
                          ImageDimension> StructuringElementType;
      StructuringElementType  element;
      element.SetRadius( radius ); 
      element.CreateStructuringElement();
      
      switch ( operation )
        {
        case 0:
          {
          typedef itk::GrayscaleDilateImageFilter<ImageType, ImageType,
            StructuringElementType >  FilterType;
          typename FilterType::Pointer  filter = FilterType::New();
          filter->SetKernel( element );
          filter->SetInput( reader->GetOutput() );
          filter->Update();
          typedef itk::ImageFileWriter<ImageType>  WriterType;
          typename WriterType::Pointer writer = WriterType::New();
          writer->SetInput( filter->GetOutput() );
          writer->SetFileName( argv[3] );
          writer->Update();
          break;
          } 
        case 1:
          {
          typedef itk::GrayscaleErodeImageFilter<ImageType, ImageType,
            StructuringElementType >  FilterType;
          typename FilterType::Pointer  filter = FilterType::New();
          filter->SetKernel( element );
          filter->SetInput( reader->GetOutput() );
          filter->Update();
          typedef itk::ImageFileWriter<ImageType>  WriterType;
          typename WriterType::Pointer writer = WriterType::New();
          writer->SetInput( filter->GetOutput() );
          writer->SetFileName( argv[3] );
          writer->Update();
          break;
          } 
        case 2:
          {
          typedef itk::GrayscaleMorphologicalClosingImageFilter<ImageType, ImageType,
            StructuringElementType >  FilterType;
          typename FilterType::Pointer  filter = FilterType::New();
          filter->SetKernel( element );
          filter->SetInput( reader->GetOutput() );
          filter->Update();
          typedef itk::ImageFileWriter<ImageType>  WriterType;
          typename WriterType::Pointer writer = WriterType::New();
          writer->SetInput( filter->GetOutput() );
          writer->SetFileName( argv[3] );
          writer->Update();
          break;
          } 
        case 3:
          {
          typedef itk::GrayscaleMorphologicalOpeningImageFilter<ImageType, ImageType,
            StructuringElementType >  FilterType;
          typename FilterType::Pointer  filter = FilterType::New();
          filter->SetKernel( element );
          filter->SetInput( reader->GetOutput() );
          filter->Update();
          typedef itk::ImageFileWriter<ImageType>  WriterType;
          typename WriterType::Pointer writer = WriterType::New();
          writer->SetInput( filter->GetOutput() );
          writer->SetFileName( argv[3] );
          writer->Update();
          break;
          } 
        default:
          {
          std::cerr << "Invalid operation choice." << std::endl; 
          return EXIT_FAILURE;
          }
        }
      }  
    }
  else
    {
    switch( operation )
      {
      case 4:
        {
        typedef itk::GrayscaleGrindPeakImageFilter<ImageType, ImageType> FilterType;
        typename FilterType::Pointer filter = FilterType::New();
        filter->SetInput( reader->GetOutput() ); 
        filter->Update();
        typedef itk::ImageFileWriter<ImageType>  WriterType;
        typename WriterType::Pointer writer = WriterType::New();
        writer->SetInput( filter->GetOutput() );
        writer->SetFileName( argv[3] );
        writer->Update();
        break;
        } 
      default:
        {
        std::cerr << "Invalid operation choice." << std::endl; 
        return EXIT_FAILURE;
        }
      } 
    }
   

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if( argc < 5 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " imageDimension inputImage outputImage operation "
      << "[radius] [type: box == 0, ball = 1, diamond = 2]" << std::endl;
    std::cerr << "  operation: " << std::endl;
    std::cerr << "    0. dilate" << std::endl;
    std::cerr << "    1. erode " << std::endl;
    std::cerr << "    2. close " << std::endl;
    std::cerr << "    3. open " << std::endl;
    std::cerr << "    4. grind peak " << std::endl;
    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) ) 
   {
   case 2:
     GrayscaleMorphology<2>( argc, argv );
     break;
   case 3:
     GrayscaleMorphology<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


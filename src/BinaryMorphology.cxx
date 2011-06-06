#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryBallStructuringElement.h" 
#include "itkBinaryBoxStructuringElement.h" 
#include "itkBinaryDiamondStructuringElement.h" 
#include "itkBinaryThinning3DImageFilter.h"
#include "itkBinaryThinningImageFilter.h"


#include "vnl/vnl_math.h"

template <unsigned int ImageDimension>
int BinaryMorphology( int argc, char * argv[] )
{
  typedef short PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType>  ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  unsigned int radius = 1;
  if ( argc > 5 )
    {
    radius = atoi( argv[5] ); 
    }  

  PixelType foreground = itk::NumericTraits<PixelType>::One;
  PixelType background = itk::NumericTraits<PixelType>::Zero;
  if ( argc > 7 )
    {
    foreground = static_cast<PixelType>( atof( argv[7] ) ); 
    }
  if ( argc > 8 )
    {
    background = static_cast<PixelType>( atof( argv[8] ) ); 
    }

  
  unsigned int operation = static_cast<unsigned int>( atoi( argv[4] ) );

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
        typedef itk::BinaryDilateImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );
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
        typedef itk::BinaryErodeImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );
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
        typedef itk::BinaryMorphologicalClosingImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetForegroundValue( foreground );
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
        typedef itk::BinaryMorphologicalOpeningImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );
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
        typedef itk::BinaryDilateImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );
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
        typedef itk::BinaryErodeImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );
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
        typedef itk::BinaryMorphologicalClosingImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetForegroundValue( foreground );
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
        typedef itk::BinaryMorphologicalOpeningImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );
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
        typedef itk::BinaryDilateImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );
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
        typedef itk::BinaryErodeImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );
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
        typedef itk::BinaryMorphologicalClosingImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetForegroundValue( foreground );
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
        typedef itk::BinaryMorphologicalOpeningImageFilter<ImageType, ImageType,
          StructuringElementType >  FilterType;
        typename FilterType::Pointer  filter = FilterType::New();
        filter->SetKernel( element );
        filter->SetInput( reader->GetOutput() );
        filter->SetBackgroundValue( background );
        filter->SetForegroundValue( foreground );
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

int BinaryThin2D( int argc, char * argv[] )
{
  typedef short PixelType;
  typedef itk::Image<PixelType, 2> ImageType;

  typedef itk::ImageFileReader<ImageType>  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();
 
  typedef itk::BinaryThinningImageFilter<ImageType, ImageType> FilterType;
  FilterType::Pointer  filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->Update();
  
  typedef itk::ImageFileWriter<ImageType>  WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( argv[3] );
  writer->Update();

  return EXIT_SUCCESS;
}

int BinaryThin3D( int argc, char * argv[] )
{
  typedef short PixelType;
  typedef itk::Image<PixelType, 3> ImageType;

  typedef itk::ImageFileReader<ImageType>  ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();
 
  typedef itk::BinaryThinning3DImageFilter<ImageType, ImageType> FilterType;
  FilterType::Pointer  filter = FilterType::New();
  filter->SetInput( reader->GetOutput() );
  filter->Update();
  
  typedef itk::ImageFileWriter<ImageType>  WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( argv[3] );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if( argc < 5 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " imageDimension inputImage outputImage operation "
      << "[radius] [type: box == 0, ball = 1, diamond = 2] [label]" << std::endl;
    std::cerr << "  operation: " << std::endl;
    std::cerr << "    0. dilate" << std::endl;
    std::cerr << "    1. erode " << std::endl;
    std::cerr << "    2. close " << std::endl;
    std::cerr << "    3. open " << std::endl;
    std::cerr << "    4. thin " << std::endl;
    return EXIT_FAILURE;
    }

  if( atoi( argv[4] ) == 4 )
    {
    switch( atoi( argv[1] ) ) 
     {
     case 2:
       BinaryThin2D( argc, argv );
       break;
     case 3:
       BinaryThin3D( argc, argv );
       break;
     default:
        std::cerr << "Unsupported dimension" << std::endl;
        exit( EXIT_FAILURE );
     }
    }
  else
    {
    switch( atoi( argv[1] ) ) 
     {
     case 2:
       BinaryMorphology<2>( argc, argv );
       break;
     case 3:
       BinaryMorphology<3>( argc, argv );
       break;
     default:
        std::cerr << "Unsupported dimension" << std::endl;
        exit( EXIT_FAILURE );
     }
   }  
}


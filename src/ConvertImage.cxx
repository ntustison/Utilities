#include <stdio.h>

#include "itkBSplineControlPointImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNumericTraits.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkFlipImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

#include "itkFDFImageIOFactory.h"
#include "itkFDFImageIO.h"

#include <string>

template <class TPixel, unsigned int ImageDimension>
int ConvertImage( int argc, char *argv[] )
{
  typedef TPixel PixelType;

  // Register FDF Factory
  itk::FDFImageIOFactory::RegisterOneFactory();

  if( atoi( argv[4] ) == 7 )
    {
    typedef itk::Vector<PixelType, ImageDimension> VectorType;
    typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;
    typedef itk::Image<PixelType, ImageDimension> ComponentImageType;

    typename DeformationFieldType::Pointer deformationField =
      DeformationFieldType::New();

    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      std::string filename = std::string( argv[2] );
      if( d == 0 )
        {
        filename += std::string( "xvec.nii.gz" );
        }
      else if( d == 1 )
        {
        filename += std::string( "yvec.nii.gz" );
        }
      else if( d == 2 )
        {
        filename += std::string( "zvec.nii.gz" );
        }

      typedef itk::ImageFileReader<ComponentImageType>
        ReaderType;
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( filename.c_str() );
      reader->Update();

      if( d == 0 )
        {
        deformationField->SetOrigin( reader->GetOutput()->GetOrigin() );
        deformationField->SetSpacing( reader->GetOutput()->GetSpacing() );
        deformationField->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
        deformationField->SetDirection( reader->GetOutput()->GetDirection() );

        deformationField->Allocate();
        VectorType V;
        V.Fill( 0.0 );

        deformationField->FillBuffer( V );
        }


      itk::ImageRegionConstIterator<ComponentImageType> It( reader->GetOutput(),
        reader->GetOutput()->GetLargestPossibleRegion() );
      itk::ImageRegionIterator<DeformationFieldType> ItD( deformationField,
        deformationField->GetLargestPossibleRegion() );
      for( It.GoToBegin(), ItD.GoToBegin(); !It.IsAtEnd(); ++ItD, ++It )
        {
        VectorType V = ItD.Get();
        V[d] = It.Get();
        ItD.Set( V );
        }
      }

    typedef itk::ImageFileWriter<DeformationFieldType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( deformationField );
    writer->SetFileName( argv[3] );
    writer->Update();
    }
  else if( atoi( argv[4] ) == 8 )
    {
    typedef itk::Vector<PixelType, ImageDimension> VectorType;
    typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;
    typedef itk::Image<PixelType, ImageDimension> ComponentImageType;

    typedef itk::ImageFileReader<DeformationFieldType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[2] );
    reader->Update();


    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      typedef itk::VectorIndexSelectionCastImageFilter<DeformationFieldType, ComponentImageType> SelectorType;
      typename SelectorType::Pointer selector = SelectorType::New();
      selector->SetInput( reader->GetOutput() );
      selector->SetIndex( d );
      selector->Update();

      std::string filename = std::string( argv[3] );
      if( d == 0 )
        {
        filename += std::string( "xvec.nii.gz" );
        }
      else if( d == 1 )
        {
        filename += std::string( "yvec.nii.gz" );
        }
      else if( d == 2 )
        {
        filename += std::string( "zvec.nii.gz" );
        }

      typedef itk::ImageFileWriter<ComponentImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( selector->GetOutput() );

      writer->SetFileName( filename.c_str() );
      writer->Update();
      }

    }
  else if( atoi( argv[4] ) == 9 )
    {
    typedef itk::Vector<PixelType, ImageDimension> VectorType;
    typedef itk::Image<VectorType, ImageDimension+1> VelocityFieldType;
    typedef itk::Image<PixelType, ImageDimension+1> ComponentImageType;

    typedef itk::ImageFileReader<VelocityFieldType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[2] );
    reader->Update();

    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      typedef itk::VectorIndexSelectionCastImageFilter<VelocityFieldType, ComponentImageType> SelectorType;
      typename SelectorType::Pointer selector = SelectorType::New();
      selector->SetInput( reader->GetOutput() );
      selector->SetIndex( d );
      selector->Update();

      std::string filename = std::string( argv[3] );
      if( d == 0 )
        {
        filename += std::string( "xvec.nii.gz" );
        }
      else if( d == 1 )
        {
        filename += std::string( "yvec.nii.gz" );
        }
      else if( d == 2 )
        {
        filename += std::string( "zvec.nii.gz" );
        }

      typedef itk::ImageFileWriter<ComponentImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( selector->GetOutput() );

      writer->SetFileName( filename.c_str() );
      writer->Update();
      }
    }
  else if( atoi( argv[4] ) == 10 )
    {
    typedef itk::Vector<PixelType, ImageDimension> VectorType;
    typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;

    typedef itk::ImageFileReader<DeformationFieldType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[2] );
    reader->Update();

    typedef itk::ImageFileWriter<DeformationFieldType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( reader->GetOutput() );

    writer->SetFileName( argv[3] );
    writer->Update();
    }
  else if( atoi( argv[4] ) < 7 )
    {
    typedef itk::Image<PixelType, ImageDimension> ImageType;
    typedef itk::ImageFileReader<ImageType> ReaderType;

    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[2] );
    reader->Update();

    // If the requested output image is a .png or .jpg, optimize for visualization
    if ( strstr( argv[3], ".png" ) != NULL || strstr( argv[3], ".jpg" ) != NULL )
      {
      typedef itk::Image<unsigned char, ImageDimension> ShortImageType;

      typedef itk::RescaleIntensityImageFilter<ImageType, ShortImageType> FilterType;
      typename FilterType::Pointer filter = FilterType::New();
      filter->SetInput( reader->GetOutput() );
      filter->SetOutputMinimum( 0 );
      filter->SetOutputMaximum( 255 );

      typedef itk::FlipImageFilter<ShortImageType> FlipperType;
      typename FlipperType::Pointer flipper = FlipperType::New();
      flipper->SetInput( filter->GetOutput() );
      typename FlipperType::FlipAxesArrayType array;
      array.Fill( false );
      array[1] = true;
      flipper->SetFlipAxes( array );

      typedef itk::ImageFileWriter<ShortImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( flipper->GetOutput() );
      writer->SetFileName( argv[3] );
      writer->Update();
      }
    else
      {
      typedef itk::ImageFileWriter<ImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( reader->GetOutput() );
      writer->SetFileName( argv[3] );
      writer->Update();
      }
    }

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension "
      << "inputImage outputImage pixelType" << std::endl;
    std::cout << "pixelType:  0 -> float (default)" << std::endl
              << "            1 -> unsigned short" << std::endl
              << "            2 -> unsigned int" << std::endl
              << "            3 -> unsigned long" << std::endl
              << "            4 -> short" << std::endl
              << "            5 -> int" << std::endl
              << "            6 -> long" << std::endl
              << "            7 -> component images to a float vector image" << std::endl
              << "            8 -> vector image to component images" << std::endl
              << "            9 -> time-varying velocity field image to component images (ImageDimension is the dimensionality of the displacement vector)" << std::endl
              << "           10 -> float vector image" << std::endl;

    exit( 0 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     if( argc > 4 && atoi( argv[4] ) == 1 )
       {
       ConvertImage<unsigned short, 2>( argc, argv );
       }
     else if( argc > 4 && atoi( argv[4] ) == 2 )
       {
       ConvertImage<unsigned int, 2>( argc, argv );
       }
     else if( argc > 4 && atoi( argv[4] ) == 3 )
       {
       ConvertImage<unsigned long, 2>( argc, argv );
       }
     else if( argc > 4 && atoi( argv[4] ) == 4 )
       {
       ConvertImage<short, 2>( argc, argv );
       }
     else if( argc > 4 && atoi( argv[4] ) == 5 )
       {
       ConvertImage<int, 2>( argc, argv );
       }
     else if( argc > 4 && atoi( argv[4] ) == 6 )
       {
       ConvertImage<long, 2>( argc, argv );
       }
     else
       {
       ConvertImage<float, 2>( argc, argv );
       }
     break;
   case 3:
     if( argc > 4 && atoi( argv[4] ) == 1 )
       {
       ConvertImage<unsigned short, 3>( argc, argv );
       }
     else if( argc > 4 && atoi( argv[4] ) == 2 )
       {
       ConvertImage<unsigned int, 3>( argc, argv );
       }
     else if( argc > 4 && atoi( argv[4] ) == 3 )
       {
       ConvertImage<unsigned long, 3>( argc, argv );
       }
     else if( argc > 4 && atoi( argv[4] ) == 4 )
       {
       ConvertImage<short, 3>( argc, argv );
       }
     else if( argc > 4 && atoi( argv[4] ) == 5 )
       {
       ConvertImage<int, 3>( argc, argv );
       }
     else if( argc > 4 && atoi( argv[4] ) == 6 )
       {
       ConvertImage<long, 3>( argc, argv );
       }
     else
       {
       ConvertImage<float, 3>( argc, argv );
       }
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

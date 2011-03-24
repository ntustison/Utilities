#include <stdio.h>

#include "itkGridImageSource.h"
#include "itkImageFileWriter.h"
#include "itkTimeProbe.h"
#include "itkCoxDeBoorBSplineKernelFunction.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc != 4 + ImageDimension*7 )
    {
    std::cout << "Usage: " << argv[0] << " outputImageFilename scale whichKernel";
    std::cout << "origin[0] ... origin[n-1] ";  
    std::cout << "size[0] ... size[n-1] ";  
    std::cout << "spacing[0] ... spacing[n-1] ";  
    std::cout << "grid_spacing[0] ... grid_spacing[n-1] ";  
    std::cout << "grid_offset[0] ... grid_offset[n-1] ";  
    std::cout << "sigma[0] ... sigma[n-1] ";
    std::cout << "whichDimension[0] ... whichDimension[n-1]" << std::endl;
    std::cout << "whichKernel = 0 -> BSpline Kernel of order 0" << std::endl; 
    std::cout << "whichKernel = 1 -> BSpline Kernel of order 1" << std::endl; 
    std::cout << "whichKernel = 2 -> BSpline Kernel of order 2" << std::endl; 
    std::cout << "whichKernel = 3 -> BSpline Kernel of order 3" << std::endl; 
    std::cout << "whichKernel = 4 -> BSpline Kernel of order 4" << std::endl; 
    exit( 1 );
    }
  typedef itk::Image<PixelType, ImageDimension> ImageType;   

  itk::TimeProbe timer;
  timer.Start();

  typedef itk::GridImageSource<ImageType> GridSourceType;
  GridSourceType::Pointer gridder = GridSourceType::New();

  double scale = atof( argv[2] );
  ImageType::SizeType size;
  ImageType::PointType origin;
  ImageType::SpacingType spacing;
  GridSourceType::ArrayType gridSpacing; 
  GridSourceType::ArrayType gridOffset; 
  GridSourceType::ArrayType sigma; 
  GridSourceType::BoolArrayType which; 
  unsigned int kernelType = atoi( argv[3] );
  
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    origin[i] = atof( argv[4+i] );
    size[i] = atoi( argv[4+ImageDimension+i] );
    spacing[i] = atof( argv[4+2*ImageDimension+i] );
    gridSpacing[i] = atof( argv[4+3*ImageDimension+i] );
    gridOffset[i] = atof( argv[4+4*ImageDimension+i] );
    sigma[i] = atof( argv[4+5*ImageDimension+i] );
    which[i] = atoi( argv[4+6*ImageDimension+i] );
    }

  typedef itk::CoxDeBoorBSplineKernelFunction<1> KernelType;
  KernelType::Pointer kernel = KernelType::New();

  if ( kernelType < 4 )
    {
    kernel->SetSplineOrder( kernelType );
    gridder->SetKernelFunction( kernel );
    }
  gridder->SetSpacing( spacing );
  gridder->SetOrigin( origin );
  gridder->SetSize( size );
  gridder->SetGridSpacing( gridSpacing );
  gridder->SetGridOffset( gridOffset );
  gridder->SetWhichDimensions( which );
  gridder->SetSigma( sigma );
  gridder->SetScale( scale );
  gridder->SetNumberOfThreads( 3 );
  gridder->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[1] );
  writer->SetInput( gridder->GetOutput() );
  writer->Update();
  
  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage outputImage" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) ) 
   {
   case 2:
     CreateGridImage<2>( argc, argv );
     break;
   case 3:
     CreateGridImage<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


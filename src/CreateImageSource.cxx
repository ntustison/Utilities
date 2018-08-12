#include <stdio.h>

#include "itkBinaryThresholdImageFilter.h"
#include "itkCoxDeBoorBSplineKernelFunction.h"
#include "itkGaussianImageSource.h"
#include "itkGridImageSource.h"
#include "itkGaborImageSource.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "Common.h"

template <unsigned int ImageDimension>
int CreateImageSource( int argc, char *argv[] )
{

  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  unsigned int type = atoi( argv[4] );

  switch( type )
    {
    case 0:  //Gaussian
      {
      typedef itk::GaussianImageSource<ImageType> GaussianSourceType;
      typename GaussianSourceType::Pointer gaussian
        = GaussianSourceType::New();

      double scale = 1.0;
      if( argc > 7 )
        {
        scale = static_cast<double>( atof( argv[7] ) );
        }
      typename GaussianSourceType::ArrayType mean;
      typename GaussianSourceType::ArrayType sigma;

      std::vector<double> m = ConvertVector<double>( argv[5] );
      std::vector<double> s = ConvertVector<double>( argv[6] );
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        mean[d] = m[d];
        sigma[d] = s[d];
        }
      gaussian->SetSpacing( reader->GetOutput()->GetSpacing() );
      gaussian->SetOrigin( reader->GetOutput()->GetOrigin() );
      gaussian->SetSize(
        reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
      gaussian->SetDirection( reader->GetOutput()->GetDirection() );
      gaussian->SetSigma( sigma );
      gaussian->SetScale( scale );
      gaussian->SetMean( mean );
      gaussian->SetNormalized( false );
      gaussian->Update();

      typedef itk::ImageFileWriter<ImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( argv[3] );
      writer->SetInput( gaussian->GetOutput() );
      writer->Update();

      break;
      }
    case 1:  //Circle
      {
      typedef itk::GaussianImageSource<ImageType> GaussianSourceType;
      typename GaussianSourceType::Pointer gaussian
        = GaussianSourceType::New();

      typename GaussianSourceType::ArrayType mean;
      typename GaussianSourceType::ArrayType sigma;

      std::vector<double> m = ConvertVector<double>( argv[5] );
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        mean[d] = m[d];
        sigma[d] = 1000.0;
        }

      gaussian->SetSpacing( reader->GetOutput()->GetSpacing() );
      gaussian->SetOrigin( reader->GetOutput()->GetOrigin() );
      gaussian->SetSize(
        reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
      gaussian->SetDirection( reader->GetOutput()->GetDirection() );
      gaussian->SetSigma( sigma );
      gaussian->SetScale( 1000.0 );
      gaussian->SetMean( mean );
      gaussian->SetNormalized( false );
      gaussian->Update();

      typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> FilterType;

      typename FilterType::Pointer filter = FilterType::New();

      filter->SetInput( gaussian->GetOutput() );

      filter->SetLowerThreshold( 1000.0 *
        vcl_exp( -0.5 * vnl_math_sqr( atof( argv[6] ) / sigma[0] ) ) );


      filter->SetUpperThreshold( itk::NumericTraits<PixelType>::max() );

      filter->SetInsideValue( 1.0 );

      filter->SetOutsideValue( 0.0 );

      filter->Update();

      typedef itk::ImageFileWriter<ImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( argv[3] );
      writer->SetInput( filter->GetOutput() );
      writer->Update();

      break;
      }
    case 2:  //Square
      {
      typename ImageType::Pointer square = ImageType::New();
      square->SetOrigin( reader->GetOutput()->GetOrigin() );
      square->SetSpacing( reader->GetOutput()->GetSpacing() );
      square->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
      square->SetDirection( reader->GetOutput()->GetDirection() );
      square->Allocate();
      square->FillBuffer( 0 );

      std::vector<double> mn = ConvertVector<double>( argv[5] );
      std::vector<double> dm = ConvertVector<double>( argv[6] );

      itk::ImageRegionIteratorWithIndex<ImageType> It( square,
        square->GetLargestPossibleRegion() );
      for( It.GoToBegin(); !It.IsAtEnd(); ++It )
        {
        typename ImageType::PointType point;
        square->TransformIndexToPhysicalPoint(
          It.GetIndex(), point );

        bool isInside = true;
        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          if( point[d] < mn[d] - 0.5*dm[d] || point[d] > mn[d] + 0.5*dm[d] )
            {
            isInside = false;
            break;
            }
          }
        if( isInside )
          {
          It.Set( 1.0 );
          }
        }
      typedef itk::ImageFileWriter<ImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( argv[3] );
      writer->SetInput( square );
      writer->Update();

      break;
      }
    case 3:  //Grid
      {
      typedef itk::GridImageSource<ImageType> GridSourceType;
      typename GridSourceType::Pointer gridder = GridSourceType::New();

      double scale = atof( argv[5] );
      typename GridSourceType::ArrayType gridSpacing;
      typename GridSourceType::ArrayType gridOffset;
      typename GridSourceType::ArrayType sigma;
      typename GridSourceType::BoolArrayType which;
      unsigned int kernelType = atoi( argv[10] );


      std::vector<double> sp = ConvertVector<double>( argv[6] );
      std::vector<double> off = ConvertVector<double>( argv[7] );
      std::vector<double> sg = ConvertVector<double>( argv[8] );
      std::vector<bool> w = ConvertVector<bool>( argv[9] );

      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        gridSpacing[d] = sp[d];
        gridOffset[d] = off[d];
        sigma[d] = sg[d];
        which[d] = w[d];
        }

      typedef itk::CoxDeBoorBSplineKernelFunction<1> KernelType;
      typename KernelType::Pointer kernel = KernelType::New();

      if ( kernelType < 5 )
        {
        kernel->SetSplineOrder( kernelType );
        gridder->SetKernelFunction( kernel );
        }
      gridder->SetSpacing( reader->GetOutput()->GetSpacing() );
      gridder->SetOrigin( reader->GetOutput()->GetOrigin() );
      gridder->SetSize(
        reader->GetOutput()->GetLargestPossibleRegion().GetSize()  );
      gridder->SetGridSpacing( gridSpacing );
      gridder->SetGridOffset( gridOffset );
      gridder->SetWhichDimensions( which );
      gridder->SetSigma( sigma );
      gridder->SetScale( scale );
      gridder->Update();

      typedef itk::ImageFileWriter<ImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( argv[3] );
      writer->SetInput( gridder->GetOutput() );
      writer->Update();

      break;
      }
    case 4:  // Gabor
      {
      typedef itk::GaborImageSource<ImageType> GaborSourceType;
      typename GaborSourceType::Pointer gabor = GaborSourceType::New();

      double frequency = atof( argv[5] );
      typename GaborSourceType::ArrayType mean;
      typename GaborSourceType::ArrayType sigma;
      bool calculateImag = static_cast<bool>( atoi( argv[6] ) );

      std::vector<double> m = ConvertVector<double>( argv[7] );
      std::vector<double> s = ConvertVector<double>( argv[8] );

      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        mean[d] = m[d];
        sigma[d] = s[d];
        }

      gabor->SetSpacing( reader->GetOutput()->GetSpacing() );
      gabor->SetOrigin( reader->GetOutput()->GetOrigin() );
      gabor->SetSize(
        reader->GetOutput()->GetLargestPossibleRegion().GetSize()  );
      gabor->SetFrequency( frequency );
      gabor->SetCalculateImaginaryPart( calculateImag );
      gabor->SetSigma( sigma );
      gabor->SetMean( mean );
      gabor->Update();

      typedef itk::ImageFileWriter<ImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( argv[3] );
      writer->SetInput( gabor->GetOutput() );
      writer->Update();

      break;
      }
    }

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << "Usage: " << argv[0] << " ImageDimension domainImage outputImage type";
    std::cout << "  type: " << std::endl;
    std::cout << "   {0=Gaussian} mean sigma [scale] " << std::endl;
    std::cout << "   {1=Circle} mean radius" << std::endl;
    std::cout << "   {2=Square} mean dimensions(in physical space)" << std::endl;
    std::cout << "   {3=Grid} scale gridSpacing gridOffset sigma whichDimension whichKernel" << std::endl;
    std::cout << "     whichKernel = 0 -> BSpline Kernel of order 0" << std::endl;
    std::cout << "     whichKernel = 1 -> BSpline Kernel of order 1" << std::endl;
    std::cout << "     whichKernel = 2 -> BSpline Kernel of order 2" << std::endl;
    std::cout << "     whichKernel = 3 -> BSpline Kernel of order 3" << std::endl;
    std::cout << "     whichKernel = 4 -> BSpline Kernel of order 4" << std::endl;
    std::cout << "     whichKernel = 5 -> Gaussian Kernel" << std::endl;
    std::cout << "   {3=Gabor} frequency calculateImaginaryPart mean sigma" << std::endl;

    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     CreateImageSource<2>( argc, argv );
     break;
   case 3:
     CreateImageSource<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}



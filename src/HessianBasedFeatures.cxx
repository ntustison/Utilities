#include <stdio.h>

#include "itkDecomposeTensorFunction.h"
#include "itkDiscreteHessianGaussianImageFunction.h"
#include "itkFixedArray.h"
#include "itkHessianToObjectnessMeasureImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkMultiScaleHessianBasedMeasureImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkVariableSizeMatrix.h"

template <unsigned int ImageDimension>
int HessianBasedFeatures( int argc, char *argv[] )
{

  /**
   * Based on the paper by Westin et al., "Geometrical Diffusion
   *    Measures for MRI from Tensor Basis Analysis"
   * and Luca Antiga's Insight Journal paper
   * http://hdl.handle.net/1926/576
   *
   */

  typedef float PixelType;
  typedef float RealType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::RescaleIntensityImageFilter<ImageType, ImageType> RescalerType;
  typename RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetInput( reader->GetOutput() );
  rescaler->SetOutputMinimum( 0.0 );
  rescaler->SetOutputMaximum( 1.0 );
  rescaler->Update();

  float sigmaMin = rescaler->GetOutput()->GetSpacing()[0];
  if( argc > 5 )
    {
    sigmaMin = atof( argv[5] );
    }
  float sigmaMax = sigmaMin;
  if( argc > 6 )
    {
    sigmaMax = atof( argv[6] );
    }
  unsigned int numberOfSigmaSteps = 1;
  if( argc > 7 )
    {
    numberOfSigmaSteps = atoi( argv[7] );
    }

  int type = atoi( argv[4] );
  float alpha = 0.5;
  float beta = 0.5;
  float gamma = 5.0;
  bool brightObject = true;

  if( argc > 8 )
    {
    brightObject = static_cast<bool>( atoi( argv[11] ) );
    }
  if( argc > 9 )
    {
    alpha = atof( argv[9] );
    }
  if( argc > 10 )
    {
    beta = atof( argv[10] );
    }
  if( argc > 11 )
    {
    gamma = atof( argv[11] );
    }

  if( type < 0 || type > 2 )
    {
    std::cerr << "Invalid feature type." << std::endl;
    return EXIT_FAILURE;
    }

    // Define the dimension of the images
  typedef typename itk::NumericTraits<PixelType>::RealType RealPixelType;

  typedef itk::SymmetricSecondRankTensor<RealPixelType, ImageDimension>
    HessianPixelType;
  typedef itk::Image<HessianPixelType, ImageDimension> HessianImageType;
  typedef itk::HessianToObjectnessMeasureImageFilter
    <HessianImageType, ImageType> ObjectnessFilterType;
  typedef itk::MultiScaleHessianBasedMeasureImageFilter
    <ImageType, HessianImageType, ImageType> MultiScaleEnhancementFilterType;

  typename MultiScaleEnhancementFilterType::Pointer multiScaleEnhancementFilter =
    MultiScaleEnhancementFilterType::New();
  multiScaleEnhancementFilter->SetInput( rescaler->GetOutput() );
  multiScaleEnhancementFilter->SetSigmaStepMethodToLogarithmic();
  multiScaleEnhancementFilter->SetSigmaMinimum( sigmaMin  );
  multiScaleEnhancementFilter->SetSigmaMaximum( sigmaMax );
  multiScaleEnhancementFilter->SetNumberOfSigmaSteps( numberOfSigmaSteps );

  typename ObjectnessFilterType::Pointer objectnessFilter = ObjectnessFilterType::New();
  objectnessFilter->SetScaleObjectnessMeasure( false );
  objectnessFilter->SetBrightObject( brightObject );
  objectnessFilter->SetAlpha( alpha );
  objectnessFilter->SetBeta( beta );
  objectnessFilter->SetGamma( gamma );
  objectnessFilter->SetObjectDimension( type );

  multiScaleEnhancementFilter->SetHessianToMeasureFilter( objectnessFilter );

  try
    {
    multiScaleEnhancementFilter->Update();
    }
  catch ( itk::ExceptionObject &e )
    {
    std::cerr << e << std::endl;
    }

  typename RescalerType::Pointer rescaler2 = RescalerType::New();
  rescaler2->SetInput( multiScaleEnhancementFilter->GetOutput() );
  rescaler2->SetOutputMinimum( 0.0 );
  rescaler2->SetOutputMaximum( 1.0 );

  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( rescaler2->GetOutput() );
  writer->SetFileName( argv[3] );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cerr << "Usage: " << argv[0]
              << " imageDimension inputImage outputImage type [sigmaMin] [sigmaMax]"
              << " [numberOfSigmaSteps]"
              << "  [brightObject] [alpha] [beta] [gamma]" << std::endl;
    std::cerr << "   Type:  0. sphere" << std::endl;
    std::cerr << "          1. line" << std::endl;
    std::cerr << "          2. plane" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     HessianBasedFeatures<2>( argc, argv );
     break;
   case 3:
     HessianBasedFeatures<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}



//
//   int type = atoi( argv[4] );
//   float alpha = 0.5;
//   float beta = 0.5;
//   float gamma = 5.0;
//   bool brightObject = true;
//
//   if( argc > 8 )
//     {
//     alpha = atof( argv[8] );
//     }
//   if( argc > 9 )
//     {
//     beta = atof( argv[9] );
//     }
//   if( argc > 10 )
//     {
//     gamma = atof( argv[10] );
//     }
//   if( argc > 11 )
//     {
//     brightObject = static_cast<bool>( atoi( argv[11] ) );
//     }
//   typedef float PixelType;
//   typedef float RealType;
//
//   typedef itk::Image<PixelType, ImageDimension> ImageType;
//   typedef itk::ImageFileReader<ImageType> ReaderType;
//   typedef itk::Image<int, ImageDimension> LabelImageType;
//
//   typedef itk::FixedArray<PixelType, ImageDimension> EigenValueArrayType;
//   typedef itk::Image<EigenValueArrayType, ImageDimension> EigenValueImageType;
//
//   typename ReaderType::Pointer reader = ReaderType::New();
//   reader->SetFileName( argv[2] );
//   reader->Update();
//
//   typename LabelImageType::Pointer maskImage = LabelImageType::New();
//   maskImage = NULL;
//   if( argc > 12 )
//     {
//     typedef itk::ImageFileReader<LabelImageType> ReaderType;
//
//     typename ReaderType::Pointer reader = ReaderType::New();
//     reader->SetFileName( argv[12] );
//     reader->Update();
//     maskImage = reader->GetOutput();
//     }
//   typename ImageType::Pointer output = ImageType::New();
//   output->SetOrigin( reader->GetOutput()->GetOrigin() );
//   output->SetDirection( reader->GetOutput()->GetDirection() );
//   output->SetSpacing( reader->GetOutput()->GetSpacing() );
//   output->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
//   output->Allocate();
//   output->FillBuffer( 0 );
//
//   float sigmaMin = reader->GetOutput()->GetSpacing()[0];
//   if( argc > 5 )
//     {
//     sigmaMin = atof( argv[5] );
//     }
//   float sigmaMax = sigmaMin;
//   if( argc > 6 )
//     {
//     sigmaMax = atof( argv[6] );
//     }
//   unsigned int numberOfSigmaSteps = 1;
//   if( argc > 7 )
//     {
//     numberOfSigmaSteps = atoi( argv[7] );
//     }
//
//   typedef itk::VariableSizeMatrix<RealType> MatrixType;
//   typedef itk::DecomposeTensorFunction<MatrixType, RealType> DecomposerType;
//   typename DecomposerType::Pointer decomposer = DecomposerType::New();
//
//   float ds = ( sigmaMax - sigmaMin ) / ( numberOfSigmaSteps - 1 );
//   if( sigmaMin == sigmaMax || numberOfSigmaSteps <= 1 )
//     {
//     ds = 1e10;
//     }
//
//   for( float sigma = sigmaMin; sigma <= sigmaMax; sigma+= ds )
//     {
//     typedef itk::DiscreteHessianGaussianImageFunction<ImageType>
//       HessianFunctionType;
//     typename HessianFunctionType::Pointer hessian = HessianFunctionType::New();
//     hessian->SetInputImage( reader->GetOutput() );
//     hessian->SetNormalizeAcrossScale( false );
//     hessian->SetUseImageSpacing( true );
//   //  hessian->SetMaximumKernelWidth( 5 );
//   //  hessian->SetMaximumError( 0.0001 );
//   //  hessian->SetInterpolationMode();
//     hessian->SetSigma( sigma );
//     hessian->Initialize();
//
//     itk::ImageRegionConstIteratorWithIndex<ImageType> It( hessian->GetInputImage(),
//       hessian->GetInputImage()->GetLargestPossibleRegion() );
//     itk::ImageRegionIterator<ImageType> ItO( output,
//       output->GetLargestPossibleRegion() );
//
//     RealType totalNumberOfPixels = static_cast<RealType>(
//       output->GetLargestPossibleRegion().GetNumberOfPixels() );
//     RealType count = 0.0;
//     std::cout << "sigma = " << sigma << ": /" << std::flush;
//     RealType oldProgress = 0.0;
//     RealType newProgress = 0.0;
//     unsigned int progressCount = 0;
//     for( It.GoToBegin(), ItO.GoToBegin(); !It.IsAtEnd(); ++It, ++ItO )
//       {
//       count++;
//       newProgress = count / totalNumberOfPixels;
//       if( newProgress - oldProgress > 0.01 )
//         {
//         std::cout << "*" << std::flush;
//         oldProgress = newProgress;
//         progressCount++;
//         if( progressCount % 10 == 0 )
//           {
//           std::cout << progressCount << std::flush;
//           }
//         }
//       if( maskImage && maskImage->GetPixel( It.GetIndex() ) == 0 )
//         {
//         continue;
//         }
//
//       typename HessianFunctionType::OutputType H
//         = hessian->EvaluateAtIndex( It.GetIndex() );
//
//       MatrixType M;
//
//      bool isValidHessian = true;
//
//       M.SetSize( ImageDimension, ImageDimension );
//       for( unsigned int m = 0; m < ImageDimension; m++ )
//         {
//         for( unsigned int n = m; n < ImageDimension; n++ )
//           {
//           M(m, n) = M(n, m) = H(m, n);
//           if( std::isnan( H(m, n) ) || std::isinf( H(m, n) )
//             || std::isinf( -H(m, n) ) )
//             {
//             isValidHessian = false;
//             break;
//             }
//           }
//         }
//
//       if( !isValidHessian )
//         {
//         continue;
//         }
//
//       MatrixType D;
//       MatrixType V;
//
//       decomposer->EvaluateSymmetricEigenDecomposition( M, D, V );
//
//       /**
//        *
//        */
//
//       itk::FixedArray<RealType, ImageDimension> lambda;
//
//       if( ImageDimension == 2 )
//         {
//         if( vnl_math_abs( D(0, 0) ) > vnl_math_abs( D(1, 1) ))
//           {
//           lambda[0] = D(0, 0);
//           lambda[1] = D(1, 1);
//           }
//         else
//           {
//           lambda[0] = D(1, 1);
//           lambda[1] = D(0, 0);
//           }
//
//         }
//       if( ImageDimension == 3 )
//         {
//         if( vnl_math_abs( D(0, 0) ) > vnl_math_abs( D(1, 1) )
//           && vnl_math_abs( D(0, 0) ) > vnl_math_abs( D(2, 2) ) )
//           {
//           lambda[0] = D(0, 0);
//           if( vnl_math_abs( D(1, 1) ) > vnl_math_abs( D(2, 2) ) )
//             {
//             lambda[1] = D(1, 1);
//             lambda[2] = D(2, 2);
//             }
//           else
//             {
//             lambda[1] = D(2, 2);
//             lambda[2] = D(1, 1);
//             }
//           }
//         if( vnl_math_abs( D(1, 1) ) > vnl_math_abs( D(2, 2) )
//           && vnl_math_abs( D(1, 1) ) > vnl_math_abs( D(0, 0) ) )
//           {
//           lambda[0] = D(1, 1);
//           if( vnl_math_abs( D(2, 2) ) > vnl_math_abs( D(0, 0) ) )
//             {
//             lambda[1] = D(2, 2);
//             lambda[2] = D(0, 0);
//             }
//           else
//             {
//             lambda[1] = D(0, 0);
//             lambda[2] = D(2, 2);
//             }
//           }
//         else
//           {
//           lambda[0] = D(2, 2);
//           if( vnl_math_abs( D(0, 0) ) > vnl_math_abs( D(1, 1) ) )
//             {
//             lambda[1] = D(0, 0);
//             lambda[2] = D(1, 1);
//             }
//           else
//             {
//             lambda[1] = D(1, 1);
//             lambda[2] = D(0, 0);
//             }
//           }
//         }
//
//       RealType S = 0;
//       for( unsigned int d = 0; d < ImageDimension; d++ )
//         {
//         S += vnl_math_sqr( lambda[d] );
//         }
//       S = vcl_sqrt( S );
//
//       RealType Ra = 0;
//       if( ImageDimension == type - 1 )
//         {
//         Ra = itk::NumericTraits<RealType>::max();
//         }
//       else
//         {
//         int M = type;
//
//         RealType numerator = vnl_math_abs( lambda[M+1-1] );
//
//         RealType denominator = 1;
//         for( unsigned int d = M+2; d <= ImageDimension; d++ )
//           {
//           denominator *= vcl_pow( vnl_math_abs( lambda[d-1] ),
//             static_cast<RealType>( 1.0 / ( ImageDimension - M ) ) );
//           }
//
//         Ra = numerator / denominator;
//         }
//
//       RealType Rb = 0;
//       if( type > 0 )
//         {
//         int M = type;
//
//         RealType numerator = vnl_math_abs( lambda[M] );
//
//         RealType denominator = 1;
//         for( unsigned int d = M+1; d <= ImageDimension; d++ )
//           {
//           denominator *= vcl_pow( vnl_math_abs( lambda[d-1] ),
//             static_cast<RealType>( 1.0 / ( ImageDimension - M ) ) );
//           }
//
//         Rb = numerator / denominator;
//         }
//
//       bool validObject = true;
//       for( unsigned int d = type; d < ImageDimension; d++ )
//         {
//         if( lambda[d] > 0 && !brightObject )
//           {
//           validObject = false;
//           break;
//           }
//         else if( lambda[d] < 0 && brightObject )
//           {
//           validObject = false;
//           break;
//           }
//         }
//
//       RealType measure = 0;
//       if( validObject == true )
//         {
//         measure = ( 1.0 - vcl_exp( -0.5 * vnl_math_sqr( Ra / alpha ) ) )
//           * vcl_exp( -0.5 * vnl_math_sqr( Rb / beta ) )
//           * ( 1.0 - vcl_exp( -0.5 * vnl_math_sqr( S / gamma ) ) );
//         }
//
//       ItO.Set( vnl_math_max( ItO.Get(), measure ) );
//       }
//     std::cout << "/" << std::endl;
//     }

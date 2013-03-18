#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMeasurementVectorTraits.h"
#include "itkGaussianRandomSpatialNeighborSubsampler.h"
#include "itkPatchBasedDenoisingImageFilter.h"

#include <string>
#include <vector>

void ConvertToLowerCase( std::string& str )
{
  std::transform( str.begin(), str.end(), str.begin(), tolower );
  // You may need to cast the above line to (int(*)(int))
  // tolower - this works as is on VC 7.1 but may not work on
  // other compilers
}

template <unsigned int ImageDimension>
int DenoiseImage( unsigned int argc, char *argv[] )
{
//   typedef float PixelType;
//   typedef itk::Image<PixelType, ImageDimension> ImageType;
//
//   typedef itk::ImageFileReader<ImageType> ReaderType;
//
//   typedef itk::PatchBasedDenoisingImageFilter<ImageType, ImageType> FilterType;
//
//   typedef typename FilterType::PatchWeightsType PatchType;
//
//   typedef itk::Statistics::GaussianRandomSpatialNeighborSubsampler<
//     typename FilterType::PatchSampleType, typename ImageType::RegionType> SamplerType;
//
//   // read the noisy image to be denoised
//   typename ReaderType::Pointer reader = ReaderType::New();
//   reader->SetFileName( argv[2] );
//   try
//     {
//     reader->Update();
//     }
//   catch( itk::ExceptionObject & excp )
//     {
//     std::cerr << "Problem encountered while reading image file : " << argv[2] << std::endl;
//     std::cerr << excp << std::endl;
//     return EXIT_FAILURE;
//     }
//
//   // create filter and initialize
//   // give image to filter and run it
//   // get filter output and write to file
//
//   typename FilterType::Pointer filter = FilterType::New();
//   filter->SetInput( reader->GetOutput() );
//
//   // patch radius is same for all dimensions of the image
//   unsigned int patchRadius = 4;
//   if( argc > 4 )
//     {
//     patchRadius = atoi( argv[4] );
//     }
//   filter->SetPatchRadius( patchRadius );
//
//   // instead of directly setting the weights, could also specify type
//   filter->UseSmoothDiscPatchWeightsOn();
//   filter->UseFastTensorComputationsOn();
//
//   std::string noiseModel = "gaussian";
//   if( argc > 5 )
//     {
//     noiseModel = std::string( argv[5] );
//     ConvertToLowerCase( noiseModel );
//     }
//   // noise model to use
//   if( noiseModel == "gaussian" )
//     {
//     filter->SetNoiseModel( FilterType::GAUSSIAN );
//     }
//   else if( noiseModel == "rician" )
//     {
//     filter->SetNoiseModel( FilterType::RICIAN );
//     }
//   else if( noiseModel == "poisson" )
//     {
//     filter->SetNoiseModel( FilterType::POISSON );
//     }
//
//   // stepsize or weight for smoothing term
//   // Large stepsizes may cause instabilities.
//   filter->SetSmoothingWeight( 1 );
//
//   float fidelityWeight = 0.1;
//   if( argc > 6 )
//     {
//     fidelityWeight = atof( argv[6] );
//     }
//   // stepsize or weight for fidelity term
//   // use a positive weight to prevent oversmoothing
//   // (penalizes deviations from noisy data based on a noise model)
//   filter->SetFidelityWeight( fidelityWeight );
//
//   // number of iterations over the image of denoising
//   unsigned int numIterations = 1;
//   if( argc > 7 )
//     {
//     numIterations = atoi( argv[7] );
//     }
//   filter->SetNumberOfIterations( numIterations );
//
//   // sampling the image to find similar patches
//   typename SamplerType::Pointer sampler = SamplerType::New();
//   // variance (in physical units) for semi-local Gaussian sampling
//   sampler->SetVariance( 400 );
//   // rectangular window restricting the Gaussian sampling
//   sampler->SetRadius( 50 ); // 2.5 * standard deviation
//   // number of random sample "patches" to use for computations
//   sampler->SetNumberOfResultsRequested( 1000 );
//   // Sampler can be complete neighborhood sampler, random neighborhood sampler,
//   // Gaussian sampler, etc.
//   filter->SetSampler( sampler );
//
//   // automatic estimation of the kernel bandwidth
//   filter->DoKernelBandwidthEstimationOn();
//   // update bandwidth every 'n' iterations
//   filter->SetKernelBandwidthUpdateFrequency( 3 );
//   // use 33% of the pixels for the sigma update calculation
//   filter->SetFractionPixelsForSigmaUpdate( 0.20 );
//
//   // multiplication factor modifying the automatically-estimated kernel sigma
//   float sigmaMultiplicationFactor = 1.0;
//   if( argc > 8 )
//     {
//     sigmaMultiplicationFactor = atof( argv[8] );
//     }
//   filter->SetSigmaMultiplicationFactor( sigmaMultiplicationFactor );
//
//
//   // manually-selected Gaussian kernel sigma
//   // filter->DoKernelBandwidthEstimationOff();
//   // typename FilterType::RealArrayType gaussianKernelSigma;
//   // gaussianKernelSigma.SetSize(reader->GetOutput()->GetNumberOfComponentsPerPixel());
//   // gaussianKernelSigma.Fill(11);
//   // filter->SetGaussianKernelSigma (gaussianKernelSigma);
//
//   // denoise the image
// //   std::cout << "Filter prior to update:\n";
// //   filter->Print( std::cout );
//   try
//     {
//     filter->Update();
//     }
//   catch( itk::ExceptionObject & excp )
//     {
//     std::cout << "Error: In " __FILE__ ", line " << __LINE__ << "\n"
//            << "Caught exception <" << excp
//            << "> while running patch-based denoising image filter."
//            << "\n\n";
//     return EXIT_FAILURE;
//     }
//
//   filter->GetOutput()->Print( std::cout, 3 );
//
//   // write the denoised image to file
//   typedef typename FilterType::OutputImageType OutputImageType;
//
//   typedef itk::ImageFileWriter<OutputImageType> WriterType;
//   typename WriterType::Pointer writer = WriterType::New();
//   writer->SetFileName( argv[3] );
//   writer->SetInput( filter->GetOutput() );
// //   try
// //   {
//     writer->Update();
// //   }
// //   catch( itk::ExceptionObject & excp )
// //   {
// //     std::cerr << excp << std::endl;
// //     return EXIT_FAILURE;
// //   }

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cerr << "Usage :  " << argv[0] << " imageDimension"
              << " inputImageFileName outputImageFileName"
              << " [patchRadius=4] [noiseModel=gaussian]"
              << " [fidelityWeight=0.1]"
              << " [numIterations=1] [sigmaMultiplicationFactor=1.0]"
              << std::endl;
    std::cout << "Notes:  " << std::endl;
    std::cout << "Noise models include poisson, rician, and gaussian. " << std::endl;
    std::cout << "Fidelity weight: step-size or weight for fidelity term use a positive weight to prevent " << std::endl;
    std::cout << "oversmoothing (penalizes deviations from noisy data based on a noise model)." << std::endl;

    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     DenoiseImage<2>( argc, argv );
     break;
   case 3:
     DenoiseImage<3>( argc, argv );
     break;
   case 4:
     DenoiseImage<4>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


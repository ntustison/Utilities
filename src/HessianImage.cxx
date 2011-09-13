#include <stdio.h>

#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkVector.h"

template <unsigned int ImageDimension>
int Tensor( int argc, char *argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::Image<unsigned int, ImageDimension> MaskImageType;

  typedef itk::Vector<PixelType, ImageDimension> CovariantVectorType;
  typedef itk::Image<CovariantVectorType, ImageDimension> CovariantVectorImageType;

  typedef itk::SymmetricSecondRankTensor<float, ImageDimension> TensorType;
  typedef itk::Image<TensorType, ImageDimension> TensorImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  bool sigma = 0.0;
  if( argc > 5 )
    {
    sigma = atof( argv[5] );
    }

  typename MaskImageType::Pointer maskImage = NULL;

  if( argc > 5 )
    {
    typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
    typename MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader->SetFileName( argv[5] );
    maskImage = maskReader->GetOutput();
    try
      {
      maskImage->Update();
      maskImage->DisconnectPipeline();
      }
    catch( ... )
      {
      maskImage = NULL;
      std::cout << "Mask not read." << std::endl;
      }
    }

  typedef itk::HessianRecursiveGaussianImageFilter<ImageType, TensorImageType> FilterType;
  typename FilterType::Pointer hessian = FilterType::New();
  hessian->SetInput( reader->GetOutput() );
  hessian->SetSigma( sigma );
  hessian->Update();

  itk::ImageRegionIteratorWithIndex<TensorImageType>
    ItH( hessian->GetOutput(),
    hessian->GetOutput()->GetLargestPossibleRegion() );

  typedef itk::ImageFileWriter<TensorImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( hessian->GetOutput() );
  writer->Update();

  if( argc > 6 )
    {
    CovariantVectorType zeroVector;
    zeroVector.Fill( 0.0 );

    std::vector<typename CovariantVectorImageType::Pointer> eigenImages;
    eigenImages.resize( ImageDimension );
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      eigenImages[d] = CovariantVectorImageType::New();
      eigenImages[d]->CopyInformation( reader->GetOutput() );
      eigenImages[d]->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
      eigenImages[d]->Allocate();
      eigenImages[d]->FillBuffer( zeroVector );
      }

    for( ItH.GoToBegin(); !ItH.IsAtEnd(); ++ItH )
      {
      if( !maskImage || maskImage->GetPixel( ItH.GetIndex() ) != 0 )
        {
        TensorType tensor = ItH.Get();

        typename TensorType::EigenValuesArrayType eigenvalues;
        typename TensorType::EigenVectorsMatrixType eigenvectors;
        tensor.ComputeEigenAnalysis( eigenvalues, eigenvectors );

        for( unsigned int d = 0; d < ImageDimension; d++ )
          {
          CovariantVectorType eigenVector;
          for( unsigned int e = 0; e < ImageDimension; e++ )
            {
            eigenVector[e] = eigenvectors( d, e );
            }
          eigenVector *= eigenvalues[d];

          eigenImages[d]->SetPixel( ItH.GetIndex(), eigenVector );
          }
        }
      }

    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      std::stringstream ss;
      ss << d;

      std::string filename = std::string( argv[6] ) + ss.str() +
        std::string( ".nii.gz" );

      typedef itk::ImageFileWriter<CovariantVectorImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( filename.c_str() );
      writer->SetInput( eigenImages[d] );
      writer->Update();
      }
    }

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " imageDimension inputImage outputTensorImage [maskImage] [sigma=0.0] [eigenImagesPrefix]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     Tensor<2>( argc, argv );
     break;
   case 3:
     Tensor<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

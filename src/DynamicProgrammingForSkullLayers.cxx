#include "itkBinaryThresholdImageFilter.h"
#include "itkBresenhamLine.h"
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLabelContourImageFilter.h"
#include "itkNeighborhoodIterator.h"
#include "itkVariableSizeMatrix.h"
#include "itkVector.h"

#include <string>
#include <vector>


template<class PixelType>
void DynamicProgramming( std::vector<PixelType> &lineProfile, std::vector<PixelType> mu, std::vector<PixelType> var  )
{

  itk::VariableSizeMatrix<PixelType> D;
  D.SetSize( mu.size(), lineProfile.size() );

  // Formulate the cost table
  for( unsigned int m = 0; m < D.Rows(); m++ )
    {
    for( unsigned int n = 0; n < D.Cols(); n++ )
      {
      if( m == 0 && n == 0 )
        {
        D(m, n) = 0.0;
        }
      else if( n == 0 )
        {
        D(m, n) = itk::NumericTraits<PixelType>::max();
        }
      else
        {
        PixelType d_l = vnl_math_sqr( lineProfile[n] - mu[m] );
//        PixelType d_l = ( 1.0 - 1.0 / vcl_sqrt( 2.0 * vnl_math::pi * var[m] ) *
//            vcl_exp( -vnl_math_sqr( lineProfile[n] - mu[m] ) / ( 2.0 * var[m] ) ) );

//        PixelType d_l = 0.0;
//        for( unsigned int l = 1; l < m; l++ )
//          {
//          d_l += vnl_math_sqr( lineProfile[l] - mu[m] );
//          }
//        d_l /= static_cast<PixelType>( m - 1 );

        PixelType cost1 = D(m, n-1) + d_l;

        PixelType cost2 = itk::NumericTraits<PixelType>::max();

        if( m > 0 )
          {
          cost2 = D(m-1, n-1) + vnl_math_sqr( lineProfile[n] - mu[m] );
//          cost2 = D(m-1, n-1) + ( 1.0 - 1.0 / vcl_sqrt( 2.0 * vnl_math::pi * var[m] ) *
//            vcl_exp( -vnl_math_sqr( lineProfile[n] - mu[m] ) / ( 2.0 * var[m] ) ) );
          }
        D(m, n) = vnl_math_min( cost1, cost2 );
        }
      }
    }

  // Find the optimal histogram mapping by tracing back through the
  // cost table.
  int m = D.Rows() - 1;
  int n = D.Cols() - 1;
  while( m > 0 || n > 0 )
    {
    lineProfile[n] = m+1;

    unsigned int minM = vnl_math_max(0, m);
    unsigned int minN = vnl_math_max(0, n - 1);
    if( D(vnl_math_max(0, m-1), vnl_math_max(0, n-1)) < D(minM, minN) )
      {
      minM = vnl_math_max(0, m - 1);
      minN = vnl_math_max(0, n - 1);
      }
    m = minM;
    n = minN;
    }
  lineProfile[n] = m+1;

}

template <unsigned int ImageDimension>
int DrawLayers( int argc, char *argv[] )
{
  const unsigned int numberOfLayers = 3;

  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef unsigned int LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
  typename LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName( argv[3] );
  labelReader->Update();

  typedef itk::Vector<unsigned int, numberOfLayers> VectorType;
  typedef itk::Image<VectorType, ImageDimension> VectorImageType;

  typename ImageType::IndexType targetIndex;

  vnl_vector<float> centerOfMass( ImageDimension );
  centerOfMass.fill( 0.0 );
  float N = 0.0;

  itk::ImageRegionIteratorWithIndex<LabelImageType> ItL(
    labelReader->GetOutput(), labelReader->GetOutput()->GetLargestPossibleRegion() );
  for( ItL.GoToBegin(); !ItL.IsAtEnd(); ++ItL )
    {
    if( ItL.Get() == 1 || ItL.Get() == 3 )
      {
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        centerOfMass[d] += ItL.GetIndex()[d];
        }
      N++;
      }
    }
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    centerOfMass[d] /= N;
    }
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    targetIndex[d] = static_cast<int>( centerOfMass[d] );
    }

  typedef itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( labelReader->GetOutput() );
  filter->SetLowerThreshold( 1 );
  filter->SetUpperThreshold( 3 );
  filter->SetInsideValue( 1 );
  filter->SetOutsideValue( 0 );
  filter->Update();

  typedef itk::LabelContourImageFilter<LabelImageType, LabelImageType> ContourFilterType;
  typename ContourFilterType::Pointer contours = ContourFilterType::New();
  contours->SetInput( labelReader->GetOutput() );
  contours->SetFullyConnected( true );
  contours->SetBackgroundValue( 0 );
  contours->Update();

  typename ContourFilterType::Pointer contours2 = ContourFilterType::New();
  contours2->SetInput( filter->GetOutput() );
  contours2->SetFullyConnected( true );
  contours2->SetBackgroundValue( 0 );
  contours2->Update();

  itk::ImageRegionIteratorWithIndex<LabelImageType> ItF(
    contours2->GetOutput(), contours2->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<LabelImageType> ItC(
    contours->GetOutput(), contours->GetOutput()->GetLargestPossibleRegion() );
  for( ItF.GoToBegin(), ItC.GoToBegin(); !ItF.IsAtEnd(); ++ItF, ++ItC )
    {
    if( ItC.Get() != 1 )
      {
      ItF.Set( 0 );
      }
    }

  typedef itk::BresenhamLine<ImageDimension> LinerType;
  LinerType liner;

  //
  // Get an initial guess of the mean and variance of the different layers
  //

  std::vector<PixelType> mu;  mu.resize( numberOfLayers );
  std::vector<PixelType> var; var.resize( numberOfLayers );
  std::vector<PixelType> pix; pix.resize( numberOfLayers );

  N = 0.0;

  for( ItF.GoToBegin(); !ItF.IsAtEnd(); ++ItF )
    {
    if( ItF.Get() == 1 )
      {
      typename ImageType::IndexType startIndex = ItF.GetIndex();

      typename LinerType::IndexArray indices = liner.BuildLine( startIndex, targetIndex );

      std::vector<PixelType> lineProfile;

      lineProfile.push_back( reader->GetOutput()->GetPixel( startIndex ) );

      typename LinerType::IndexArray::const_iterator it;
      for( it = indices.begin(); it != indices.end(); it++ )
        {
        if( labelReader->GetOutput()->GetPixel( *it ) == 2 )
          {
          break;
          }
        lineProfile.push_back( reader->GetOutput()->GetPixel( *it ) );
        }

      if( lineProfile.size() > 2 )
        {
        pix[0] = 0.5 * ( lineProfile[0] + lineProfile[lineProfile.size()-1] );
        pix[1] = lineProfile[static_cast<unsigned int>( 0.5*(lineProfile.size()))];
        pix[2] = pix[0];

        N += 1.0;
        for( unsigned int d = 0; d < numberOfLayers; d++ )
          {
          mu[d] = mu[d] * ( N - 1.0 ) / N + pix[d] / N;
          }
        if ( N > 1.0 )
          {
          for( unsigned int d = 0; d < numberOfLayers; d++ )
            {
            var[d] = var[d] * ( N - 1.0 ) / N + vnl_math_sqr( pix[d] - mu[d] ) / ( N - 1.0 );
            }
          }
        }
      }
    }

  for( unsigned int d = 0; d < numberOfLayers; d++ )
    {
    std::cout << mu[d] << " " << var[d] << std::endl;
    }

  VectorType zeroVector( 0.0 );

  typename VectorImageType::Pointer countImage = VectorImageType::New();
  countImage->CopyInformation( reader->GetOutput() );
  countImage->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  countImage->Allocate();
  countImage->FillBuffer( zeroVector );

  for( ItF.GoToBegin(); !ItF.IsAtEnd(); ++ItF )
    {
    if( ItF.Get() == 1 )
      {
      typename ImageType::IndexType startIndex = ItF.GetIndex();

      typename LinerType::IndexArray indices = liner.BuildLine( startIndex, targetIndex );

      std::vector<PixelType> lineProfile;

      typename LinerType::IndexArray::const_iterator it;
      for( it = indices.begin(); it != indices.end(); it++ )
        {
        if( labelReader->GetOutput()->GetPixel( *it ) == 2 )
          {
          break;
          }
        lineProfile.push_back( reader->GetOutput()->GetPixel( *it ) );
        }

      if( lineProfile.size() >= numberOfLayers )
        {
        DynamicProgramming<PixelType>( lineProfile, mu, var );

        typename LinerType::IndexArray::const_iterator it;
        for( it = indices.begin(); it != indices.end(); it++ )
          {
          if( labelReader->GetOutput()->GetPixel( *it ) == 2 )
            {
            break;
            }
          VectorType count = countImage->GetPixel( *it );
          count[static_cast<unsigned int>( lineProfile[it-indices.begin()] )-1] += 1;
          countImage->SetPixel( *it, count );
          }
        }
      }
    }

  contours2->GetOutput()->FillBuffer( 0 );

  itk::ImageRegionIteratorWithIndex<VectorImageType> ItV(
    countImage, countImage->GetLargestPossibleRegion() );
  for( ItF.GoToBegin(), ItV.GoToBegin(); !ItF.IsAtEnd(); ++ItF, ++ItV )
    {
    if( labelReader->GetOutput()->GetPixel( ItF.GetIndex() ) == 1 ||
      labelReader->GetOutput()->GetPixel( ItF.GetIndex() ) == 3 )
      {
      VectorType countVector = ItV.Get();

      unsigned int maxCount = countVector[0];
      LabelType maxLabel = 1;

      for( unsigned int d = 1; d < numberOfLayers; d++ )
        {
        if( maxCount < countVector[d] )
          {
          maxLabel = static_cast<LabelType>( d + 1 );
          }
        }
      ItF.Set( maxLabel );
      }
    }

  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[4] );
  writer->SetInput( contours2->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << argv[0] << " imageDimension inputImage labeledMask outputImage" << std::endl;
    exit( 0 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     DrawLayers<2>( argc, argv );
     break;
   case 3:
     DrawLayers<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


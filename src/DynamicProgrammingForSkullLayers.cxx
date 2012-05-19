#include "itkBinaryThresholdImageFilter.h"
#include "itkBresenhamLine.h"
#include "itkGradientImageFilter.h"
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLabelContourImageFilter.h"
#include "itkNeighborhoodIterator.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkVariableSizeMatrix.h"
#include "itkVector.h"

#include <string>
#include <vector>


template<class PixelType>
void DynamicProgramming( std::vector<PixelType> &lineProfile, std::vector<PixelType> mu, std::vector<PixelType> var  )
{
  unsigned int lineProfileSize = lineProfile.size();

  if( lineProfileSize < 3 && lineProfileSize > 20 )
    {
    lineProfile.resize( 0 );
    return;
    }

  std::vector<int> lineProfileLabels;
  lineProfileLabels.resize( lineProfileSize );

  for( int n = 0; n < lineProfileSize; n++ )
    {
    if( lineProfile[n] == 0 )
      {
      lineProfile[n] = 1;
      }
    else
      {
      break;
      }
    }

  for( int n = 0; n < lineProfileSize; n++ )
    {
    if( lineProfile[lineProfileSize - n - 1] == 0 )
      {
      lineProfile[lineProfileSize - n - 1] = 1;
      }
    else
      {
      break;
      }
    }

  unsigned int edgeCount = 0;
  for( int n = 0; n < lineProfileSize; n++ )
    {
    if( lineProfile[n] == 0 && lineProfile[n+1] == 1 )
      {
      edgeCount++;
      }
    }

  if( edgeCount != 3 )
    {
    lineProfile.resize( 0 );
    return;
    }
  else
    {
    edgeCount = 0;
    for( unsigned int n = 0; n < lineProfileSize - 1; n++ )
      {
      lineProfileLabels[n] = vnl_math_min( 3, static_cast<int>( edgeCount + 1 ) );
      if( lineProfile[n] == 0 && lineProfile[n+1] == 1 )
        {
        edgeCount++;
        }
      }
    lineProfileLabels[lineProfileSize-1] = 3;
    }

  for( unsigned int n = 0; n < lineProfileSize; n++ )
    {
    lineProfile[n] = lineProfileLabels[n];
    }

//   for( unsigned int n = 0; n < lineProfile.size(); n++ )
//     {
//     std::cout << lineProfile[n] << std::flush;
//     }
//   std::cout << std::endl;


/*
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
//         PixelType d_l = vnl_math_sqr( lineProfile[n] - mu[m] );
       PixelType d_l = ( 1.0 - 1.0 / vcl_sqrt( 2.0 * vnl_math::pi * var[m] ) *
           vcl_exp( -vnl_math_sqr( lineProfile[n] - mu[m] ) / ( 2.0 * var[m] ) ) );

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
//           cost2 = D(m-1, n-1) + vnl_math_sqr( lineProfile[n] - mu[m] );
         cost2 = D(m-1, n-1) + ( 1.0 - 1.0 / vcl_sqrt( 2.0 * vnl_math::pi * var[m] ) *
           vcl_exp( -vnl_math_sqr( lineProfile[n] - mu[m] ) / ( 2.0 * var[m] ) ) );
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
      minM = vnl_math_max( 0, m - 1 );
      minN = vnl_math_max( 0, n - 1 );
      }
    m = minM;
    n = minN;
    }
  lineProfile[n] = m + 1;
*/
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

  typedef itk::BinaryThresholdImageFilter<LabelImageType, LabelImageType> ThresholderType;
  typename ThresholderType::Pointer thresholder = ThresholderType::New();
  thresholder->SetInput( labelReader->GetOutput() );
  thresholder->SetLowerThreshold( 3 );
  thresholder->SetUpperThreshold( 3 );
  thresholder->SetInsideValue( 1 );
  thresholder->SetOutsideValue( 0 );

  typedef itk::SignedMaurerDistanceMapImageFilter<LabelImageType, ImageType> DistancerType;
  typename DistancerType::Pointer distancer = DistancerType::New();
  distancer->SetInput( thresholder->GetOutput() );
  distancer->SetSquaredDistance( false );
  distancer->SetUseImageSpacing( true );
  distancer->SetInsideIsPositive( false );
  distancer->Update();

  typedef itk::GradientImageFilter<ImageType> GradientType;
  typename GradientType::Pointer grad = GradientType::New();
  grad->SetInput( distancer->GetOutput() );
  grad->UseImageDirectionOff();
  grad->Update();

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


  VectorType zeroVector( 0.0 );

  typename VectorImageType::Pointer countImage = VectorImageType::New();
  countImage->CopyInformation( reader->GetOutput() );
  countImage->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  countImage->Allocate();
  countImage->FillBuffer( zeroVector );

  typename LabelImageType::Pointer output = LabelImageType::New();
  output->CopyInformation( contours2->GetOutput() );
  output->SetRegions( contours2->GetOutput()->GetRequestedRegion() );
  output->Allocate();
  output->FillBuffer( 0 );

  for( ItF.GoToBegin(); !ItF.IsAtEnd(); ++ItF )
    {
    if( ItF.Get() == 1 )
      {
      typename ImageType::IndexType startIndex = ItF.GetIndex();

      typename GradientType::OutputPixelType vector = grad->GetOutput()->GetPixel( startIndex );

      typename LinerType::LType direction;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        direction[d] = -vector[d];
        }
      typename LinerType::OffsetArray offsets = liner.BuildLine( direction, 100 );

      std::vector<PixelType> lineProfile;

      lineProfile.push_back( reader->GetOutput()->GetPixel( startIndex ) );

      typename LinerType::OffsetArray::const_iterator it;
      for( it = offsets.begin(); it != offsets.end(); it++ )
        {
        if( labelReader->GetOutput()->GetPixel( startIndex + *it ) == 2 ||  labelReader->GetOutput()->GetPixel( startIndex + *it ) == 0 )
          {
          break;
          }
        lineProfile.push_back( reader->GetOutput()->GetPixel( startIndex + *it ) );
        }

      if( lineProfile.size() >= numberOfLayers )
        {
        DynamicProgramming<PixelType>( lineProfile, mu, var );

        typename LinerType::OffsetArray::const_iterator it;
        for( it = offsets.begin(); it != offsets.end(); it++ )
          {
          if( it-offsets.begin() == lineProfile.size() )
            {
            break;
            }
          output->SetPixel( startIndex + *it, lineProfile[it-offsets.begin()] );
          }
        }
      }
    }


//   typedef itk::NeighborhoodIterator<LabelImageType> NIteratorType;
//   typename NIteratorType::RadiusType radius;
//   radius.Fill( 1 );
//
//   NIteratorType ItN( radius, output, output->GetRequestedRegion() );
//
//   for( ItN.GoToBegin(); !ItN.IsAtEnd(); ++ItN )
//     {
//     if( ItN.GetCenterPixel() == 0 && labelReader->GetOutput()->GetPixel( ItN.GetIndex() ) == 1  )
//       {
//       vnl_vector<unsigned int> count;
//       count.set_size( 3 );
//
//       count.fill( 0.0 );
//
//       for( unsigned int n = 0; n < ( ItN.GetNeighborhood() ).Size(); n++ )
//         {
//         if( ItN.GetPixel( n ) > 0 )
//           {
//           count[ItN.GetPixel( n ) - 1]++;
//           }
//         }
//       if( count[0] != 0 && count[1] == 0 && count[2] == 0 )
//         {
//         ItN.SetCenterPixel( 1 );
//         }
//       if( count[1] != 0 && count[2] == 0 && count[0] == 0 )
//         {
//         ItN.SetCenterPixel( 2 );
//         }
//       if( count[2] != 0 && count[0] == 0 && count[1] == 0 )
//         {
//         ItN.SetCenterPixel( 3 );
//         }
//       }
//     }

  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[4] );
  writer->SetInput( output );
  writer->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 5 )
    {
    std::cout << argv[0] << " imageDimension inputImage labeledMask outputImage [outputSampleFile]" << std::endl;
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


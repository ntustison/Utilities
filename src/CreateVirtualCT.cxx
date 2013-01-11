#include "itkBresenhamLine.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLabelContourImageFilter.h"
#include "itkNeighborhoodIterator.h"
#include "itkVector.h"

#include <string>
#include <vector>

template<class TValue>
TValue Convert( std::string optionString )
			{
			TValue value;
			std::istringstream iss( optionString );
			iss >> value;
			return value;
			}

template<class TValue>
std::vector<TValue> ConvertVector( std::string optionString )
			{
			std::vector<TValue> values;
			std::string::size_type crosspos = optionString.find( 'x', 0 );

			if ( crosspos == std::string::npos )
					{
					values.push_back( Convert<TValue>( optionString ) );
					}
			else
					{
					std::string element = optionString.substr( 0, crosspos ) ;
					TValue value;
					std::istringstream iss( element );
					iss >> value;
					values.push_back( value );
					while ( crosspos != std::string::npos )
							{
							std::string::size_type crossposfrom = crosspos;
							crosspos = optionString.find( 'x', crossposfrom + 1 );
							if ( crosspos == std::string::npos )
									{
									element = optionString.substr( crossposfrom + 1, optionString.length() );
									}
							else
									{
									element = optionString.substr( crossposfrom + 1, crosspos ) ;
									}
							std::istringstream iss( element );
							iss >> value;
							values.push_back( value );
							}
					}
			return values;
			}


template <unsigned int ImageDimension>
int DrawLines( int argc, char *argv[] )
{
  typedef float PixelType;
  typedef unsigned int LabelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );

  typename ImageType::Pointer mri = reader->GetOutput();
  mri->Update();
  mri->DisconnectPipeline();

  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
  typename LabelReaderType::Pointer labelReader = LabelReaderType::New();
  labelReader->SetFileName( argv[3] );

  typename LabelImageType::Pointer seg = labelReader->GetOutput();
  seg->Update();
  seg->DisconnectPipeline();

  typename ImageType::IndexType targetIndex;

  vnl_vector<float> centerOfMass( ImageDimension );
  centerOfMass.fill( 0.0 );
  float N = 0.0;

  itk::ImageRegionIteratorWithIndex<ImageType> It(
    mri, mri->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    typename ImageType::IndexType index = It.GetIndex();
    PixelType weight = It.Get();

    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      centerOfMass[d] += ( weight * index[d] );
      }
    N += weight;
    }
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    centerOfMass[d] /= N;
    }
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    targetIndex[d] = static_cast<int>( centerOfMass[d] );
    }
  std::cout << "Target index = " << targetIndex << std::endl;

  typename ImageType::Pointer output = ImageType::New();
  output->CopyInformation( mri );
  output->SetRegions( mri->GetLargestPossibleRegion() );
  output->Allocate();
  output->FillBuffer( 0 );

  typename ImageType::Pointer outputCount = ImageType::New();
  outputCount->CopyInformation( mri );
  outputCount->SetRegions( mri->GetLargestPossibleRegion() );
  outputCount->Allocate();
  outputCount->FillBuffer( 0 );

  typedef itk::BresenhamLine<ImageDimension> LinerType;
  LinerType liner;
  typedef typename LinerType::LType VectorType;
  typedef typename LinerType::OffsetType OffsetType;
  typedef typename LinerType::IndexType IndexType;

  typename ImageType::SizeType size = mri->GetLargestPossibleRegion().GetSize();

  unsigned long maxLength = 0;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    maxLength += size[d] * size[d];
    }
  maxLength = static_cast<unsigned long>( vcl_sqrt( maxLength ) );

  itk::ImageRegionIteratorWithIndex<LabelImageType> ItL(
    seg, seg->GetLargestPossibleRegion() );
  for( ItL.GoToBegin(); !ItL.IsAtEnd(); ++ItL )
    {
    if( ItL.Get() == 1 )
      {
      typename ImageType::IndexType startIndex = ItL.GetIndex();

      VectorType direction;
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        direction[d] = startIndex[d] - targetIndex[d];
        }
      direction.Normalize();

      typename LinerType::OffsetArray offsets = liner.BuildLine( direction, maxLength );

      IndexType currentIndex = targetIndex;

      bool isFound = false;
      unsigned int begOffsetIndex = 0;
      unsigned int endOffsetIndex = 0;
      IndexType begIndex;
      IndexType endIndex;

      typename LinerType::OffsetArray::const_iterator it;

      for( it = offsets.begin(); it != offsets.end(); ++it )
        {
        if( begOffsetIndex == 0 && seg->GetPixel( currentIndex ) == 3 )
          {
          begOffsetIndex = it - offsets.begin();
          begIndex = currentIndex;
          }
        if( !seg->GetLargestPossibleRegion().IsInside( currentIndex ) )
          {
          break;
          }
        if( startIndex == currentIndex && begOffsetIndex != 0 )
          {
          if( seg->GetPixel( targetIndex + offsets[it - offsets.begin() + 1] ) == 0 )
            {
            isFound = true;
            endOffsetIndex = it - offsets.begin();
            endIndex = currentIndex;
            break;
            }
          }
        currentIndex = targetIndex + *it;
        }

      if( !isFound )
        {
        continue;
        }

      // Calculate virtual total thickness, Y
      // Y_CT = 1.62396 + 0.55682 * X_T1

      typename ImageType::PointType begPoint;
      typename ImageType::PointType endPoint;

      output->TransformIndexToPhysicalPoint( begIndex, begPoint );
      output->TransformIndexToPhysicalPoint( endIndex, endPoint );

      float X_T1 = begPoint.EuclideanDistanceTo( endPoint );
      float Y_CT = 1.62396 + 0.55682 * X_T1;

      unsigned int minBegOffsetIndex = begOffsetIndex;
      unsigned int minEndOffsetIndex = endOffsetIndex;
      float minDistanceDifference = vnl_math_abs( X_T1 - Y_CT );
      for( unsigned int t = 1; t <= 10; t++ )
        {
        typename ImageType::PointType begPointPre;
        output->TransformIndexToPhysicalPoint( targetIndex + offsets[begOffsetIndex-t], begPointPre );
        typename ImageType::PointType endPointPost;
        output->TransformIndexToPhysicalPoint( targetIndex + offsets[endOffsetIndex+t-1], endPointPost );

        float distanceDifference = vnl_math_abs( Y_CT - begPointPre.EuclideanDistanceTo( endPointPost ) );
        if( distanceDifference < minDistanceDifference )
          {
          minDistanceDifference = distanceDifference;
          minBegOffsetIndex = begOffsetIndex - t;
          minEndOffsetIndex = endOffsetIndex + t - 1;
          }

        output->TransformIndexToPhysicalPoint( targetIndex + offsets[endOffsetIndex+t], endPointPost );
        distanceDifference = vnl_math_abs( Y_CT - begPointPre.EuclideanDistanceTo( endPointPost ) );
        if( distanceDifference < minDistanceDifference )
          {
          minDistanceDifference = distanceDifference;
          minBegOffsetIndex = begOffsetIndex - t;
          minEndOffsetIndex = endOffsetIndex + t;
          }
        }

      float Y_CT_INNER = 0.00010 + 0.30800 * X_T1;
      float Y_CT_MIDDLE = 0.73595 + 0.18319 * X_T1;
      float Y_CT_OUTER = 0.87618 + 0.06697 * X_T1;

      float Y_totalDistance = Y_CT_INNER + Y_CT_MIDDLE + Y_CT_OUTER;

      float Y_PORTION_INNER = Y_CT_INNER / Y_totalDistance;
      float Y_PORTION_MIDDLE = Y_CT_MIDDLE / Y_totalDistance;
      float Y_PORTION_OUTER = Y_CT_OUTER / Y_totalDistance;

      unsigned int numberOfOffsets = minEndOffsetIndex - minBegOffsetIndex + 1;

      unsigned int numberOfInnerOffsets = static_cast<unsigned int>( Y_PORTION_INNER * ( numberOfOffsets ) + 0.5 );
      unsigned int numberOfOuterOffsets = static_cast<unsigned int>( Y_PORTION_OUTER * ( numberOfOffsets ) + 0.5 );

      unsigned int numberOfMiddleOffsets = 0;
      if( numberOfOffsets > ( numberOfInnerOffsets + numberOfOuterOffsets ) )
        {
        numberOfMiddleOffsets = numberOfOffsets - ( numberOfInnerOffsets + numberOfOuterOffsets );
        }

      float averageMR = 0.0;
      for( unsigned int n = 0; n < numberOfOffsets; n++ )
        {
        averageMR += mri->GetPixel( targetIndex + offsets[minBegOffsetIndex + n] );
        }
      averageMR /= static_cast<float>( numberOfOffsets );

      float densityCT = 1528.719 - 29.072 * vcl_sqrt( averageMR );

      float densityCT_INNER = ( 1.0 - Y_PORTION_INNER ) * densityCT;
      float densityCT_MIDDLE = ( 1.0 - Y_PORTION_MIDDLE ) * densityCT;
      float densityCT_OUTER = ( 1.0 - Y_PORTION_OUTER ) * densityCT;

      for( unsigned int n = 0; n < numberOfOffsets; n++ )
        {
        typename ImageType::IndexType index1 = targetIndex + offsets[minBegOffsetIndex + n];

        float outputValue = output->GetPixel( index1 );
        float outputCountValue = outputCount->GetPixel( index1 );
        outputCount->SetPixel( index1, outputCountValue + 1.0 );
        if( n < numberOfInnerOffsets )
          {
          output->SetPixel( index1, outputValue + densityCT_INNER );
          }
        else if( n >= numberOfInnerOffsets && n < numberOfInnerOffsets + numberOfMiddleOffsets )
          {
          output->SetPixel( index1, outputValue + densityCT_MIDDLE );
          }
        else // if( n >= numberOfInnerOffsets + numberOfMiddleOffsets && n < numberOfOuterOffsets )
          {
          output->SetPixel( index1, outputValue + densityCT_OUTER );
          }
        }

      // Debug
//       if( isFound )
//         {
//         output->SetPixel( targetIndex + offsets[minBegOffsetIndex], 1 );
//         output->SetPixel( targetIndex + offsets[minEndOffsetIndex], 2 );
//         }
      }
    }

  itk::ImageRegionIterator<ImageType> ItO( output, output->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> ItC( outputCount, outputCount->GetLargestPossibleRegion() );
  for( ItO.GoToBegin(), ItC.GoToBegin(); !ItO.IsAtEnd(); ++ItO, ++ItC )
    {
    if( ItC.Get() > 0 )
      {
      ItO.Set( ItO.Get() / ItC.Get() );
      }
    }

  // Use Poisson diffusion to fill in holes

  typename ImageType::Pointer outputFilled;

  typedef itk::ImageDuplicator<ImageType> DuplicatorType;
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage( output );

  outputFilled = duplicator->GetOutput();
  outputFilled->Update();
  outputFilled->DisconnectPipeline();

  itk::ImageRegionIterator<ImageType> ItF( outputFilled, outputFilled->GetLargestPossibleRegion() );

  float delta = itk::NumericTraits<float>::max();

  unsigned int maxIterations = 100;
  unsigned int iteration = 0;

  while( iteration++ < maxIterations && delta >= 1e-4 )
    {
    std::cout << "iteration " << iteration << ": delta = " << delta << std::endl;
    typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> SmoothingFilterType;
    typename SmoothingFilterType::Pointer smoothingFilter = SmoothingFilterType::New();

    smoothingFilter->SetUseImageSpacing( true );
    smoothingFilter->SetVariance( 1.0 );
    smoothingFilter->SetMaximumError( 0.01f );
    smoothingFilter->SetInput( outputFilled );
    smoothingFilter->Update();

    typename ImageType::Pointer smoothImage = smoothingFilter->GetOutput();

    itk::ImageRegionIterator<ImageType> ItS( smoothImage, smoothImage->GetLargestPossibleRegion() );

    delta = 0;
    N = 0.0;
    for( ItO.GoToBegin(), ItS.GoToBegin(), ItL.GoToBegin(), ItF.GoToBegin(); !ItO.IsAtEnd(); ++ItO, ++ItS, ++ItL, ++ItF )
      {
      if( ItL.Get() != 0 && ItO.Get() == 0.0 )
        {
        delta += vnl_math_abs( ItS.Get() - ItF.Get() );
        ItF.Set( ItS.Get() );
        N++;
        }
      else
        {
        ItF.Set( ItO.Get() );
        }
      }
    delta /= N;
    }


  typedef itk::ImageFileWriter<ImageType> WriterType;
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
    std::cout << argv[0] << " imageDimension inputT1 labelMask outputVirtualCT" << std::endl;
    exit( 0 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     DrawLines<2>( argc, argv );
     break;
   case 3:
     DrawLines<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


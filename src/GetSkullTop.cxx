#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBresenhamLine.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLabelContourImageFilter.h"
#include "itkRobustAutomaticThresholdImageFilter.h"

int main( int argc, char *argv[] )
{
  if( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " inputImage outputImage" << std::endl;
    }

  const unsigned int ImageDimension = 3;

  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  ImageType::Pointer input = reader->GetOutput();

  ImageType::IndexType startIndex = input->GetLargestPossibleRegion().GetIndex();

  ImageType::IndexType endIndex = startIndex;
  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    endIndex[d] += ( input->GetLargestPossibleRegion().GetSize()[d] - 1 );
    }

  vnl_vector<float> avgZDistances( 2 * ImageDimension, 0 );
  for( unsigned int d = 0 ; d < ImageDimension; d++ )
    {
    for( unsigned int s = 0; s <= 1; s++ )
      {
      for( unsigned int c = 0; c < ( 2 << ( ImageDimension - 2 ) ); c++ )
        {
        ImageType::IndexType cornerIndex;

        unsigned int nonDimensionCount = 0;
        for( unsigned int dd = 0; dd < ImageDimension; dd++ )
          {
          if( dd == d )
            {
            if( s == 0 )
              {
              cornerIndex[dd] = startIndex[dd];
              }
            else
              {
              cornerIndex[dd] = endIndex[dd];
              }
            }
          else
            {
            if( c & ( 1 << nonDimensionCount ) )
              {
              cornerIndex[dd] = endIndex[dd];
              }
            else
              {
              cornerIndex[dd] = startIndex[dd];
              }
            nonDimensionCount++;
            }
          }
        ImageType::PointType cornerPoint;
        input->TransformIndexToPhysicalPoint( cornerIndex, cornerPoint );

        avgZDistances[2*d + s] += cornerPoint[ImageDimension-1];

//         std::cout << cornerIndex << ", " << cornerPoint << std::endl;
        }
      avgZDistances[2*d + s] /= static_cast<float>( 2 << ( ImageDimension - 2 ) );
      }
    }

unsigned int whichDimension = ( avgZDistances.arg_max() ) / 2;
unsigned int offset = ( avgZDistances.arg_max() ) - 2 * whichDimension;

ImageType::PointType topCenterPoint;
topCenterPoint.Fill( 0.0 );

for( unsigned int c = 0; c < ( 2 << ( ImageDimension - 2 ) ); c++ )
  {
  ImageType::IndexType cornerIndex;

  unsigned int nonDimensionCount = 0;
  for( unsigned int dd = 0; dd < ImageDimension; dd++ )
    {
    if( dd == whichDimension )
      {
      if( offset == 0 )
        {
        cornerIndex[dd] = startIndex[dd];
        }
      else
        {
        cornerIndex[dd] = endIndex[dd];
        }
      }
    else
      {
      if( c & ( 1 << nonDimensionCount ) )
        {
        cornerIndex[dd] = endIndex[dd];
        }
      else
        {
        cornerIndex[dd] = startIndex[dd];
        }
      nonDimensionCount++;
      }
    }
  ImageType::PointType cornerPoint;
  input->TransformIndexToPhysicalPoint( cornerIndex, cornerPoint );

  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    topCenterPoint[d] += ( cornerPoint[d] / static_cast<float>( 2 << ( ImageDimension - 2 ) ) );
    }
  }

ImageType::IndexType topCenterIndex;
input->TransformPhysicalPointToIndex( topCenterPoint, topCenterIndex );

/*************************************************************/
/*  Calculate Otsu threshold                                 */
/*************************************************************/

typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> SmoothingFilterType;
SmoothingFilterType::Pointer smoothingFilter = SmoothingFilterType::New();
smoothingFilter->SetInput( input );
smoothingFilter->SetUseImageSpacingOff();
smoothingFilter->SetVariance( 4*4 );
smoothingFilter->SetMaximumError( 0.01f );
smoothingFilter->Update();

typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<ImageType, ImageType> GradientType;
GradientType::Pointer gradFilter = GradientType::New();
gradFilter->SetInput( smoothingFilter->GetOutput() );
gradFilter->SetSigma( 3 );
gradFilter->Update();

typedef itk::RobustAutomaticThresholdImageFilter <ImageType, ImageType> ThresholdImageFilterType;
ThresholdImageFilterType::Pointer histFilter = ThresholdImageFilterType::New();
histFilter->SetInput( smoothingFilter->GetOutput() );
histFilter->SetGradientImage( gradFilter->GetOutput() );
histFilter->SetInsideValue( 1 );
histFilter->SetOutsideValue( 0 );
histFilter->SetPow( 2 );
histFilter->Update();

typedef itk::BinaryBallStructuringElement<PixelType, ImageDimension> StructuringElementType;
StructuringElementType element;
element.SetRadius( 5 );
element.CreateStructuringElement();

typedef itk::BinaryMorphologicalClosingImageFilter<ImageType, ImageType, StructuringElementType> MorphFilterType;
MorphFilterType::Pointer morphFilter = MorphFilterType::New();
morphFilter->SetKernel( element );
morphFilter->SetInput( histFilter->GetOutput() );
morphFilter->SetForegroundValue( 1 );

typedef itk::LabelContourImageFilter<ImageType, ImageType> ContourFilterType;
ContourFilterType::Pointer contourFilter = ContourFilterType::New();
contourFilter->SetInput( morphFilter->GetOutput() );
contourFilter->SetFullyConnected( true );
contourFilter->Update();

ImageType::Pointer output = ImageType::New();
output->CopyInformation( input );
output->SetRegions( input->GetLargestPossibleRegion() );
output->Allocate();
output->FillBuffer( 0 );

typedef itk::BresenhamLine<ImageDimension> LinerType;
LinerType liner;

itk::ImageRegionIteratorWithIndex<ImageType> It( contourFilter->GetOutput(),
  contourFilter->GetOutput()->GetLargestPossibleRegion() );
for( It.GoToBegin(); !It.IsAtEnd(); ++It )
  {
  if( It.Get() == 1 )
    {
    ImageType::IndexType startIndex = It.GetIndex();

    LinerType::IndexArray indices = liner.BuildLine( startIndex, topCenterIndex );
    bool overlap = false;
    LinerType::IndexArray::const_iterator it;
    for( it = indices.begin() + 1; it != indices.end(); it++ )
      {
      if( output->GetLargestPossibleRegion().IsInside( *it ) )
        {
        if( contourFilter->GetOutput()->GetPixel( *it ) == 1 )
          {
          overlap = true;
          break;
          }
        }
      else
        {
        break;
        }
      }

    if( !overlap )
      {
      if( output->GetLargestPossibleRegion().IsInside( indices[1] ) )
        {
        output->SetPixel( indices[1], 1 );
        }
//       for( it = indices.begin(); it != indices.end(); it++ )
//         {
//         if( output->GetLargestPossibleRegion().IsInside( *it ) )
//           {
//           output->SetPixel( *it, 1 );
//           }
//         }
      }
    }
  }

typedef itk::ImageFileWriter<ImageType> WriterType;
WriterType::Pointer writer = WriterType::New();
writer->SetInput( output );
writer->SetFileName( argv[2] );
writer->Update();

return 0;
}


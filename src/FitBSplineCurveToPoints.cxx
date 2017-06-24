#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkPointSet.h"

#include <stdio.h>
#include <vector>
#include <fstream>
#include <string>

template <unsigned int PointDimension>
int FitBSplineCurveToPoints( unsigned int argc, char *argv[] )
{
  typedef double RealType;
  typedef int LabelType;
  typedef itk::Image<LabelType, PointDimension> LabelImageType;

  typedef itk::ImageFileReader<LabelImageType> LabelImageReaderType;
  typename LabelImageReaderType::Pointer reader = LabelImageReaderType::New();
  reader->SetFileName( argv[2] );

  typename LabelImageType::Pointer inputImage = reader->GetOutput();
  inputImage->Update();
  inputImage->DisconnectPipeline();

  typedef itk::Vector<RealType, PointDimension> VectorType;
  typedef itk::Image<VectorType, 1> CurveImageType;

  typedef itk::PointSet<VectorType, 1> PointSetType;
  typename PointSetType::Pointer pointSet = PointSetType::New();
  pointSet->Initialize();

  typedef itk::BSplineScatteredDataPointSetToImageFilter
     <PointSetType, CurveImageType>  FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  typename FilterType::WeightsContainerType::Pointer weights = FilterType::WeightsContainerType::New();


  typename itk::ImageSliceIteratorWithIndex<LabelImageType>
    It( inputImage, inputImage->GetRequestedRegion() );

  LabelType whichLabel = atoi( argv[4] );
  unsigned int whichAxis = atoi( argv[5] );

  if( whichAxis == 0 )
    {
    It.SetFirstDirection( 1 );
    It.SetSecondDirection( 2 );
    }
  else if( whichAxis == 1 )
    {
    It.SetFirstDirection( 2 );
    It.SetSecondDirection( 0 );
    }
  else
    {
    It.SetFirstDirection( 0 );
    It.SetSecondDirection( 1 );
    }

  RealType minIndex = itk::NumericTraits<RealType>::max();
  RealType maxIndex = itk::NumericTraits<RealType>::NonpositiveMin();
  unsigned long count = 0;
  typename LabelImageType::PointType averagePoint( 0.0 );

  // std::cout << "Iterating." << std::endl;

  It.GoToBegin();
  while( !It.IsAtEnd() )
    {
    RealType indexValue = static_cast<RealType>( It.GetIndex()[whichAxis] );

    while( !It.IsAtEndOfSlice() )
      {
      while( !It.IsAtEndOfLine() )
        {
        LabelType label = It.Get();

        if( label == whichLabel )
          {
          if( indexValue < minIndex )
            {
            minIndex = indexValue;
            }
          if( indexValue > maxIndex )
            {
            maxIndex = indexValue;
            }
          typename LabelImageType::PointType imagePoint;
          for( unsigned int d = 0; d < PointDimension; d++ )
            {
            imagePoint[d] = It.GetIndex()[d];
            }
          averagePoint += imagePoint.GetVectorFromOrigin();

          pointSet->SetPointData( count, imagePoint.GetVectorFromOrigin() );

          typename PointSetType::PointType point;
          point[0] = indexValue;
          pointSet->SetPoint( count, point );
          weights->InsertElement( count, 1.0 );
          count++;

          It.Set( 0 );
          }
        ++It;
        }
      It.NextLine();
      }
    It.NextSlice();
    }

  for( unsigned int d = 0; d < PointDimension; d++ )
    {
    averagePoint[d] /= static_cast<RealType>( count );
    }

  for ( unsigned int i = 0; i < pointSet->GetNumberOfPoints(); i++ )
    {
    typename PointSetType::PointType point;
    pointSet->GetPoint( i, &point );
    point[0] = ( point[0] - minIndex ) / ( maxIndex - minIndex );
    pointSet->SetPoint( i, point );

    typename PointSetType::PixelType pointData;
    pointSet->GetPointData( i, &pointData );
    pointData -= averagePoint.GetVectorFromOrigin();
    pointSet->SetPointData( i, pointData );
    }

  // std::cout << "Fit the model." << std::endl;

  filter->SetInput( pointSet );
  filter->SetGenerateOutputImage( true );

  typename CurveImageType::PointType origin;
  origin.Fill( 0.0 );
  filter->SetOrigin( origin );
  typename CurveImageType::SpacingType spacing;
  spacing[0] = 0.0001;
  filter->SetSpacing( spacing );
  typename CurveImageType::SizeType size;
  size[0] = static_cast<unsigned int>( 1.0 / spacing[0] + 1 );
  filter->SetSize( size );
  typename FilterType::ArrayType order;
  order[0] = 3;
  filter->SetSplineOrder( order );
  typename FilterType::ArrayType ncps;
  ncps[0] = order[0] + 1;
  if ( argc > 7 )
    {
    ncps[0] = atoi( argv[7] );
    }
  filter->SetNumberOfControlPoints( ncps );
  typename FilterType::ArrayType nlevels;
  nlevels[0] = 5;
  if ( argc > 6 )
    {
    nlevels[0] = atoi( argv[6] );
    }
  filter->SetNumberOfLevels( nlevels );
  typename FilterType::ArrayType close;
  close[0] = 0;
  filter->SetCloseDimension( close );
  filter->Update();

  // std::cout << "Place curve in image." << std::endl;

  itk::ImageRegionIterator<CurveImageType> ItB(
    filter->GetOutput(), filter->GetOutput()->GetLargestPossibleRegion() );
  for ( ItB.GoToBegin(); !ItB.IsAtEnd(); ++ItB )
    {
    VectorType vector = ItB.Get() + averagePoint.GetVectorFromOrigin();
    typename LabelImageType::IndexType index;
    for( unsigned int d = 0; d < PointDimension; d++ )
      {
      index[d] = vector[d];
      }
    inputImage->SetPixel( index, whichLabel );
    }

  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( inputImage );
  writer->Update();

  return 0;
}


int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputLabelImage outputLabelImage"
      << "whichLabel whichAxis [nlevels=5] [numberOfControlPoints=4]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 3:
     FitBSplineCurveToPoints<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

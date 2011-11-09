#include "itkBresenhamLine.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"

#include <vector>
#include <fstream>

template <unsigned int ImageDimension>
int Convert( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " inputFile outputFile " << std::endl;
    exit( 1 );
    }

  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[3] );
  reader->Update();

  typedef unsigned int LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;

  typename LabelImageType::Pointer labelImage = LabelImageType::New();
  labelImage->CopyInformation( reader->GetOutput() );
  labelImage->SetRegions( reader->GetOutput()->GetBufferedRegion() );
  labelImage->Allocate();
  labelImage->FillBuffer( 0 );

  typedef float RealType;

  typedef typename LabelImageType::PointType PointType;

  std::vector<PointType> points;
  std::vector<LabelType> labels;
  std::vector<LabelType> distinctLabels;

  RealType x[3];
  LabelType l;

  std::fstream str( argv[2] );

  while ( str >> x[0] >> x[1] >> x[2] >> l )
    {
    if ( x[0] == 0 && x[1] == 0 && x[2] == 0 && l == 0 )
      {
      continue;
      }

    PointType point;
    for( unsigned int d = 0; d < ImageDimension; d++ )
      {
      point[d] = x[d];
      }
    points.push_back( point );

    if( std::find( distinctLabels.begin(), distinctLabels.end(), l ) == distinctLabels.end() )
      {
      distinctLabels.push_back( l );
      }
    labels.push_back( l );
    }
  str.close();

  std::cout << "HERE" << std::endl;

  for( l = 0; l < distinctLabels.size(); l++ )
    {
    typename LabelImageType::IndexType index;
    typename LabelImageType::IndexType sourceIndex;
    typename LabelImageType::IndexType targetIndex;

    bool firstOneFound = false;
    for( unsigned int n = 0; n < points.size(); n++ )
      {
      if( labels[n] == distinctLabels[l] )
        {
        labelImage->TransformPhysicalPointToIndex( points[n], index );

        if( firstOneFound )
          {
          targetIndex = sourceIndex;
          sourceIndex = index;

          if( sourceIndex != targetIndex )
            {
            typedef itk::BresenhamLine<ImageDimension> LinerType;
            typedef typename LinerType::IndexType IndexType;

            LinerType liner;

            typename LinerType::IndexArray indices = liner.BuildLine( sourceIndex, targetIndex );

            typename LinerType::IndexArray::const_iterator it;
            for( it = indices.begin(); it != indices.end(); it++ )
              {
              if( labelImage->GetPixel( *it ) == 0 )
                {
                labelImage->SetPixel( *it, labels[n] );
                }
              }
            }
          else
            {
            labelImage->SetPixel( sourceIndex, labels[n] );
            }
          }
        else
          {
          sourceIndex = index;
          firstOneFound = true;
          }
        }
      }
    }

  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[4] );
  writer->SetInput( labelImage );
  writer->Update();

  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " imageDimension inputAvantsFile referenceImage outputImage" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     Convert<2>( argc, argv );
     break;
   case 3:
     Convert<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

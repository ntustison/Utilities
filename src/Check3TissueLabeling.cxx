#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkLabelStatisticsImageFilter.h"

#include <string>

int main( int argc, char* argv[] )
{
  if( argc < 8 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " priorWarped1 priorWarped2 priorWarped3 segmentation posteriorWarped1 posteriorWarped2 posteriorWarped3" << std::endl;
    exit( 1 );
    }

  typedef double PixelType;
  typedef unsigned int LabelType;

  const unsigned int Dimension = 3;

  typedef itk::Image<PixelType, Dimension>           ImageType;
  typedef itk::ImageFileReader<ImageType>            ReaderType;
  typedef itk::Image<LabelType, Dimension>           LabelImageType;
  typedef itk::ImageFileReader<LabelImageType>       LabelReaderType;

		typedef itk::ImageFileReader<LabelImageType> LabelReaderType;
	 LabelReaderType::Pointer labelImageReader = LabelReaderType::New();
		labelImageReader->SetFileName( argv[4] );

		LabelImageType::Pointer labelImage = labelImageReader->GetOutput();
		labelImage->Update();
		labelImage->DisconnectPipeline();

  ImageType::Pointer priors[Dimension];
  ImageType::Pointer posteriors[Dimension];

  for( unsigned int d = 0; d < Dimension; d++ )
    {
    ReaderType::Pointer priorReader = ReaderType::New();
    priorReader->SetFileName( argv[d+1] );

    priors[d] = priorReader->GetOutput();
    priors[d]->Update();
    priors[d]->DisconnectPipeline();

    ReaderType::Pointer posteriorReader = ReaderType::New();
    posteriorReader->SetFileName( argv[d+5] );

    posteriors[d] = posteriorReader->GetOutput();
    posteriors[d]->Update();
    posteriors[d]->DisconnectPipeline();
    }

  LabelType movingLabels[3];
  LabelType fixedLabels[3];
  for( unsigned int d = 0; d < 3; d++ )
    {
    typedef itk::LabelStatisticsImageFilter<ImageType, LabelImageType> HistogramGeneratorType;
    HistogramGeneratorType::Pointer stats = HistogramGeneratorType::New();
    stats->SetInput( priors[d] );
    stats->SetLabelInput( labelImage );
    stats->Update();

    LabelType maxLabel = 1;
    for( LabelType l = 2; l <= 3; l++ )
      {
      if( stats->GetMean( l ) > stats->GetMean( maxLabel ) )
        {
        maxLabel = l;
        }
      }
    movingLabels[d] = maxLabel;
    fixedLabels[d] = d + 1;
    }

  for( unsigned int d = 0; d < 3; d++ )
    {
    std::cout << fixedLabels[d] << " -> " << movingLabels[d] << std::endl;

    bool foundLabel = false;
    for( LabelType l = 1; l <= 3; l++ )
      {
      if( movingLabels[d] == l )
        {
        foundLabel = true;
        }
      }
    if( !foundLabel )
      {
      std::cerr << "Not all labels were found." << std::endl;
      return EXIT_FAILURE;
      }
    }


  bool writeSegmentationImage = false;
  for( unsigned int d = 0; d < 3; d++ )
    {
    LabelType movingLabel = movingLabels[d];
    LabelType fixedLabel = fixedLabels[d];

    if( movingLabel == fixedLabel )
      {
      continue;
      }
    else
      {
      writeSegmentationImage = true;

      std::cout << "Writing posteriors " << movingLabels[d] << " " << argv[5+d] << std::endl;
      std::cout << "Writing posteriors " << movingLabels[movingLabels[d]-1] << " " << argv[5+movingLabels[d]-1] << std::endl;

      typedef itk::ImageFileWriter<ImageType> WriterType;
      WriterType::Pointer writer1 = WriterType::New();
      writer1->SetInput( posteriors[movingLabels[d]-1] );
      writer1->SetFileName( argv[5 + d] );
      writer1->Update();

      WriterType::Pointer writer2 = WriterType::New();
      writer2->SetInput( posteriors[movingLabels[movingLabels[d]-1] - 1] );
      writer2->SetFileName( argv[5+movingLabels[d]-1] );
      writer2->Update();

      LabelType tmp = movingLabels[d];
      movingLabels[d] = movingLabels[tmp-1];
      movingLabels[tmp-1] = tmp;
      }

    itk::ImageRegionIterator<LabelImageType> ItL( labelImage, labelImage->GetRequestedRegion() );
    for( ItL.GoToBegin(); !ItL.IsAtEnd(); ++ItL )
      {
      LabelType currentLabel = ItL.Get();
      if( currentLabel == movingLabel )
        {
        ItL.Set( fixedLabel );
        }
      else if( currentLabel == fixedLabel )
        {
        ItL.Set( movingLabel );
        }
      }
    }

  if( writeSegmentationImage )
    {
    typedef itk::ImageFileWriter<LabelImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput( labelImage );
    writer->SetFileName( argv[4] );
    writer->Update();
    }

  return EXIT_SUCCESS;
}

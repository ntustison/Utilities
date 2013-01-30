#include "itkImageDuplicator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkLabelOverlapMeasuresImageFilter.h"
#include "itkMatrix.h"

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

  LabelImageType::Pointer maxPriorLabelImage = LabelImageType::New();
  maxPriorLabelImage->CopyInformation( labelImage );
  maxPriorLabelImage->SetRegions( labelImage->GetRequestedRegion() );
  maxPriorLabelImage->Allocate();
  maxPriorLabelImage->FillBuffer( 0 );

  itk::ImageRegionIteratorWithIndex<LabelImageType> ItL( labelImage, labelImage->GetRequestedRegion() );
  itk::ImageRegionIterator<LabelImageType> ItM( maxPriorLabelImage, maxPriorLabelImage->GetRequestedRegion() );
  for( ItL.GoToBegin(), ItM.GoToBegin(); !ItL.IsAtEnd(); ++ItM, ++ItL )
    {
    if( ItL.Get() == 0 )
      {
      continue;
      }

    LabelType maxLabel = 1;
    PixelType maxPrior = priors[0]->GetPixel( ItL.GetIndex() );
    for( LabelType d = 2; d <= 3; d++ )
      {
      PixelType prior = priors[d-1]->GetPixel( ItL.GetIndex() );
      if( prior > maxPrior )
        {
        maxPrior = prior;
        maxLabel = d;
        }
      }
    ItM.Set( maxLabel );
    }

  itk::Matrix<LabelType, 6, 3> permutations;

  unsigned int which = 0;

  permutations(which, 0) = 1;
  permutations(which, 1) = 2;
  permutations(which, 2) = 3;

  permutations(++which, 0) = 1;
  permutations(which  , 1) = 3;
  permutations(which  , 2) = 2;

  permutations(++which, 0) = 2;
  permutations(which  , 1) = 1;
  permutations(which  , 2) = 3;

  permutations(++which, 0) = 2;
  permutations(which  , 1) = 3;
  permutations(which  , 2) = 1;

  permutations(++which, 0) = 3;
  permutations(which  , 1) = 1;
  permutations(which  , 2) = 2;

  permutations(++which, 0) = 3;
  permutations(which  , 1) = 2;
  permutations(which  , 2) = 1;

  PixelType maxDice = 0.0;
  int maxPermutationRow = -1;
  for( unsigned r = 0; r < 6; r++ )
    {
    typedef itk::ImageDuplicator<LabelImageType> DuplicatorType;
    DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator->SetInputImage( labelImage );
    duplicator->Update();

    LabelImageType::Pointer permutedLabelImage = duplicator->GetOutput();

    itk::ImageRegionIterator<LabelImageType> ItP( permutedLabelImage, permutedLabelImage->GetRequestedRegion() );
    for( ItP.GoToBegin(); !ItP.IsAtEnd(); ++ItP )
      {
      LabelType permutedLabel = ItP.Get();
      if( permutedLabel != 0 )
        {
        unsigned int whichColumn = permutedLabel - 1;
        ItP.Set( permutations( r, whichColumn ) );
        }
      }

    typedef itk::LabelOverlapMeasuresImageFilter<LabelImageType> FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetSourceImage( permutedLabelImage );
    filter->SetTargetImage( maxPriorLabelImage );
    filter->Update();

    PixelType dice = filter->GetMeanOverlap();
    std::cout << r << ": " << dice << std::endl;
    if( dice > maxDice )
      {
      maxPermutationRow = r;
      maxDice = dice;
      }
    }

  if( maxPermutationRow == -1 )
    {
    std::cerr << "Unexpected problem." << std::endl;
    return EXIT_FAILURE;
    }

  LabelType movingLabels[3];
  LabelType fixedLabels[3];
  for( unsigned int d = 0; d < 3; d++ )
    {
    std::cout << d+1 << " -> " << permutations( maxPermutationRow, d ) << std::endl;
    movingLabels[d] = permutations( maxPermutationRow, d );
    fixedLabels[d] = d + 1;
    }

  if( maxPermutationRow == 0 )
    {
    std::cout << "No need to change labels/posteriors." << std::endl;
    }
  else
    {
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
      if( d == 0 )
        {
        typedef itk::ImageDuplicator<LabelImageType> DuplicatorType;
        DuplicatorType::Pointer duplicator = DuplicatorType::New();
        duplicator->SetInputImage( labelImage );
        duplicator->Update();

        LabelImageType::Pointer permutedLabelImage = duplicator->GetOutput();

        itk::ImageRegionIterator<LabelImageType> ItP( permutedLabelImage, permutedLabelImage->GetRequestedRegion() );
        for( ItP.GoToBegin(); !ItP.IsAtEnd(); ++ItP )
          {
          LabelType permutedLabel = ItP.Get();
          if( permutedLabel != 0 )
            {
            unsigned int whichColumn = permutedLabel - 1;
            ItP.Set( permutations( maxPermutationRow, whichColumn ) );
            }
          }

        typedef itk::ImageFileWriter<LabelImageType> WriterType;
        WriterType::Pointer writer = WriterType::New();
        writer->SetInput( permutedLabelImage );
        writer->SetFileName( argv[4] );
        writer->Update();
        }
      }
    }

  return EXIT_SUCCESS;
}




//   double means[3] = { 0.0 };
//
//   LabelType movingLabels[3];
//   LabelType fixedLabels[3];
//   for( unsigned int d = 0; d < 3; d++ )
//     {
//     typedef itk::LabelStatisticsImageFilter<ImageType, LabelImageType> HistogramGeneratorType;
//     HistogramGeneratorType::Pointer stats = HistogramGeneratorType::New();
//     stats->SetInput( priors[d] );
//     stats->SetLabelInput( labelImage );
//     stats->Update();
//
//     std::cout << "Mean for prior " << d + 1 << std::endl;
//     for( unsigned int e = 0; e < 3; e++ )
//       {
//       std::cout << "  " << e+1 << ": " << stats->GetMean( e+1 ) << std::endl;
//       }
//
//     LabelType maxLabel = 1;
//     for( LabelType l = 2; l <= 3; l++ )
//       {
//
//       if( stats->GetMean( l ) > stats->GetMean( maxLabel ) )
//         {
//         maxLabel = l;
//         }
//       }
//     movingLabels[d] = maxLabel;
//     fixedLabels[d] = d + 1;
//     }
//
//   for( unsigned int d = 0; d < 3; d++ )
//     {
//     std::cout << fixedLabels[d] << " -> " << movingLabels[d] << std::endl;
//     }
//
//   for( LabelType l = 1; l <= 3; l++ )
//     {
//     bool foundLabel = false;
//     for( unsigned int d = 0; d < 3; d++ )
//       {
//       if( movingLabels[d] == l )
//         {
//         foundLabel = true;
//         }
//       }
//     if( !foundLabel )
//       {
//       std::cerr << "Not all labels were found." << std::endl;
//       return EXIT_FAILURE;
//       }
//     }
//
//   bool writeSegmentationImage = false;
//   for( unsigned int d = 0; d < 3; d++ )
//     {
//     LabelType movingLabel = movingLabels[d];
//     LabelType fixedLabel = fixedLabels[d];
//
//     if( movingLabel == fixedLabel )
//       {
//       continue;
//       }
//     else
//       {
//       writeSegmentationImage = true;
//
//       std::cout << "Writing posteriors " << movingLabels[d] << " " << argv[5+d] << std::endl;
//       std::cout << "Writing posteriors " << movingLabels[movingLabels[d]-1] << " " << argv[5+movingLabels[d]-1] << std::endl;
//
//       typedef itk::ImageFileWriter<ImageType> WriterType;
//       WriterType::Pointer writer1 = WriterType::New();
//       writer1->SetInput( posteriors[movingLabels[d]-1] );
//       writer1->SetFileName( argv[5 + d] );
//       writer1->Update();
//
//       WriterType::Pointer writer2 = WriterType::New();
//       writer2->SetInput( posteriors[movingLabels[movingLabels[d]-1] - 1] );
//       writer2->SetFileName( argv[5+movingLabels[d]-1] );
//       writer2->Update();
//
//       LabelType tmp = movingLabels[d];
//       movingLabels[d] = movingLabels[tmp-1];
//       movingLabels[tmp-1] = tmp;
//       }
//
//     itk::ImageRegionIterator<LabelImageType> ItL( labelImage, labelImage->GetRequestedRegion() );
//     for( ItL.GoToBegin(); !ItL.IsAtEnd(); ++ItL )
//       {
//       LabelType currentLabel = ItL.Get();
//       if( currentLabel == movingLabel )
//         {
//         ItL.Set( fixedLabel );
//         }
//       else if( currentLabel == fixedLabel )
//         {
//         ItL.Set( movingLabel );
//         }
//       }
//     }
//
//   if( writeSegmentationImage )
//     {
//     typedef itk::ImageFileWriter<LabelImageType> WriterType;
//     WriterType::Pointer writer = WriterType::New();
//     writer->SetInput( labelImage );
//     writer->SetFileName( argv[4] );
//     writer->Update();
//     }

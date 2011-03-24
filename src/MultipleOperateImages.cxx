#include "itkArray.h"
#include "itkImageDuplicator.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileWriter.h"

#include <string>
#include <vector>

template <unsigned int ImageDimension>
int MultipleOperateImages( int argc, char * argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType>  ReaderType;

  typedef unsigned int LabelType;
  typedef itk::Image<LabelType, ImageDimension> LabelImageType;
  typedef itk::ImageFileReader<LabelImageType> LabelReaderType;

  /**
   * list the files
   */
  std::cout << "Using the following files: " << std::endl;
  for( unsigned int n = 4; n < static_cast<unsigned int>( argc ); n++ )
    {
    std::cout << "   " << n-3 << ": " << argv[n] << std::endl;
    }

  std::string op = std::string( argv[2] );

  if( op.compare( std::string( "mean" ) ) == 0 )
    {
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[4] );
    reader->Update();

    typename ImageType::Pointer output = reader->GetOutput();

    float N = 1.0;

    for( unsigned int n = 5; n < static_cast<unsigned int>( argc ); n++ )
      {
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( argv[n] );
      reader->Update();

      N += 1.0;
      itk::ImageRegionIterator<ImageType> It( reader->GetOutput(),
        reader->GetOutput()->GetLargestPossibleRegion() );
      itk::ImageRegionIterator<ImageType> ItO( output,
        output->GetLargestPossibleRegion() );
      for( It.GoToBegin(), ItO.GoToBegin(); !It.IsAtEnd(); ++It, ++ItO )
        {
        ItO.Set( ItO.Get() * ( N - 1.0 ) / N + It.Get() / N );
        }
      }

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( output );
    writer->SetFileName( argv[3] );
    writer->Update();
    }
  else if( op.compare( std::string( "var" ) ) == 0 )
    {
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[4] );
    reader->Update();

    typename ImageType::Pointer meanImage = reader->GetOutput();

    float N = 1.0;

    typename ImageType::Pointer output = ImageType::New();
    output->SetOrigin( meanImage->GetOrigin() );
    output->SetSpacing( meanImage->GetSpacing() );
    output->SetRegions( meanImage->GetLargestPossibleRegion() );
    output->SetDirection( meanImage->GetDirection() );
    output->Allocate();
    output->FillBuffer( 0 );

    for( unsigned int n = 5; n < static_cast<unsigned int>( argc ); n++ )
      {
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( argv[n] );
      reader->Update();


      N += 1.0;

      itk::ImageRegionIterator<ImageType> It( reader->GetOutput(),
        reader->GetOutput()->GetLargestPossibleRegion() );
      itk::ImageRegionIterator<ImageType> ItO( output,
        output->GetLargestPossibleRegion() );
      itk::ImageRegionIterator<ImageType> ItM( meanImage,
        meanImage->GetLargestPossibleRegion() );
      for( It.GoToBegin(), ItO.GoToBegin(), ItM.GoToBegin(); !It.IsAtEnd(); ++It, ++ItO, ++ItM )
        {
        ItM.Set( ItM.Get() * ( N - 1.0 ) / N + It.Get() / N );
        if ( N > 1.0 )
          {
          ItO.Set( ItO.Get() * ( N - 1.0 )/N +
            ( It.Get() - ItM.Get() )*( It.Get() - ItM.Get() ) / ( N - 1.0 ) );
          }
        }
      }

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( output );
    writer->SetFileName( argv[3] );
    writer->Update();
    }
  else if( op.compare( std::string( "s" ) ) == 0 )
    {
    std::vector<typename LabelImageType::Pointer> labelImages;
    for( unsigned int n = 4; n < static_cast<unsigned int>( argc ); n++ )
      {
      typename LabelReaderType::Pointer reader = LabelReaderType::New();
      reader->SetFileName( argv[n] );
      reader->Update();
      labelImages.push_back( reader->GetOutput() );
      }

    typename ImageType::Pointer output = ImageType::New();
    output->SetOrigin( labelImages[0]->GetOrigin() );
    output->SetSpacing( labelImages[0]->GetSpacing() );
    output->SetRegions( labelImages[0]->GetLargestPossibleRegion() );
    output->SetDirection( labelImages[0]->GetDirection() );
    output->Allocate();
    output->FillBuffer( 0 );

    /**
     * Get the actual labels---assume that the first image read has all
     * the actual labels.
     */
    itk::ImageRegionIteratorWithIndex<LabelImageType> It( labelImages[0],
      labelImages[0]->GetLargestPossibleRegion() );

    std::vector<LabelType> labels;
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      if( std::find( labels.begin(), labels.end(), It.Get() ) == labels.end() )
        {
        labels.push_back( It.Get() );
        }
      }
    itk::Array<unsigned long> labelCount( labels.size() );
    std::vector<float> maxProbabilities;
    std::vector<typename LabelImageType::IndexType> maxProbabilityIndices;
    std::vector<float> maxProbabilityCount;
    maxProbabilities.resize( labels.size() );
    maxProbabilityIndices.resize( labels.size() );
    maxProbabilityCount.resize( labels.size() );

    std::cout << "Labels: ";
    for( unsigned int n = 0; n < labels.size(); n++ )
      {
      maxProbabilities[n] = 0.0;
      maxProbabilityCount[n] = 0.0;
      std::cout << labels[n] << " ";
      }
    std::cout << std::endl;


    std::vector<LabelType>::iterator loc
      = std::find( labels.begin(), labels.end(), 0 );
    if( loc == labels.end() )
      {
      std::cerr << "Background label 0 not found." << std::endl;
      return EXIT_FAILURE;
      }
    unsigned int backgroundIndex = static_cast<unsigned int>(
      std::distance( labels.begin(), loc ) );

    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      labelCount.Fill( 0 );
      for( unsigned int n = 0; n < labelImages.size(); n++ )
        {
        LabelType label = labelImages[n]->GetPixel( It.GetIndex() );
        std::vector<LabelType>::iterator loc
          = std::find( labels.begin(), labels.end(), label );
        if( loc == labels.end() )
          {
          std::cerr << "Label " << label << " not found." << std::endl;
          return EXIT_FAILURE;
          }
        unsigned int index = static_cast<unsigned int>(
          std::distance( labels.begin(), loc ) );
        labelCount[index]++;
        }

      if( labelCount[backgroundIndex] == labelImages.size() )
        {
        continue;
        }

      float totalProbability = 0.0;
      for( unsigned int m = 0; m < labelCount.Size(); m++ )
        {
        if( m == backgroundIndex )
          {
          continue;
          }
        float probability = static_cast<float>( labelCount[m] ) /
          static_cast<float>( labelImages.size() );

        if( probability > maxProbabilities[m] )
          {
          maxProbabilities[m] = probability;
          }
        for( unsigned int n = 0; n < labelCount.Size(); n++ )
          {
          if( m == n || n == backgroundIndex )
            {
            continue;
            }
          probability *= ( 1.0 - static_cast<float>( labelCount[n] ) /
          static_cast<float>( labelImages.size() ) );
          }
        totalProbability += probability;
        }
      output->SetPixel( It.GetIndex(), totalProbability );
      }

    /**
     * Find the central cluster of the max probabilities
     */
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      labelCount.Fill( 0 );
      for( unsigned int n = 0; n < labelImages.size(); n++ )
        {
        LabelType label = labelImages[n]->GetPixel( It.GetIndex() );
        std::vector<LabelType>::iterator loc
          = std::find( labels.begin(), labels.end(), label );
        if( loc == labels.end() )
          {
          std::cerr << "Label " << label << " not found." << std::endl;
          return EXIT_FAILURE;
          }
        unsigned int index = static_cast<unsigned int>(
          std::distance( labels.begin(), loc ) );
        labelCount[index]++;
        }

      if( labelCount[backgroundIndex] == labelImages.size() )
        {
        continue;
        }

      for( unsigned int m = 0; m < labelCount.Size(); m++ )
        {
        if( m == backgroundIndex )
          {
          continue;
          }
        float probability = static_cast<float>( labelCount[m] ) /
          static_cast<float>( labelImages.size() );
        if( probability == maxProbabilities[m] )
          {
          for( unsigned int d = 0; d < ImageDimension; d++ )
            {
            if( maxProbabilityCount[m] == 0 )
              {
              maxProbabilityIndices[m][d] = It.GetIndex()[d];
              }
            else
              {
              maxProbabilityIndices[m][d] += It.GetIndex()[d];
              }
            }
          maxProbabilityCount[m]++;
          }
        }
      }

    for( unsigned int n = 0; n < maxProbabilityIndices.size(); n++ )
      {
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        maxProbabilityIndices[n][d] = vcl_floor( maxProbabilityIndices[n][d] /
          maxProbabilityCount[n] );
        }
      }

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( output );
    writer->SetFileName( argv[3] );
    writer->Update();

    for( unsigned int m = 0; m < labels.size(); m++ )
      {
      if( labels[m] == 0 )
        {
        continue;
        }
      std::cout << labels[m] << " "
        << maxProbabilityIndices[m] << " "
        << maxProbabilities[m] << std::endl;
      }
    }
  else if( op.compare( std::string( "w" ) ) == 0 )
    {
    std::vector<typename ImageType::Pointer> images;
    for( unsigned int n = 4; n < static_cast<unsigned int>( argc ); n++ )
      {
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( argv[n] );
      reader->Update();
      images.push_back( reader->GetOutput() );
      }

    typename ImageType::Pointer output = ImageType::New();
    output->SetOrigin( images[0]->GetOrigin() );
    output->SetSpacing( images[0]->GetSpacing() );
    output->SetRegions( images[0]->GetLargestPossibleRegion() );
    output->SetDirection( images[0]->GetDirection() );
    output->Allocate();
    output->FillBuffer( 0 );

    itk::ImageRegionIteratorWithIndex<ImageType> ItO( output,
      output->GetLargestPossibleRegion() );
    for( ItO.GoToBegin(); !ItO.IsAtEnd(); ++ItO )
      {
      float probability = 0.0;
      for( unsigned int i = 0; i < images.size(); i++ )
        {
        float negation = 1.0;
        for( unsigned int j = 0; j < images.size(); j++ )
          {
          if( i == j )
            {
            continue;
            }
          negation *= ( 1.0 - images[j]->GetPixel( ItO.GetIndex() ) );
          }
        probability += negation * images[i]->GetPixel( ItO.GetIndex() );
        }
      ItO.Set( probability );
      }
    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( output );
    writer->SetFileName( argv[3] );
    writer->Update();
    }
  else if( op.compare( std::string( "ex" ) ) == 0 )
    {
    std::vector<typename ImageType::Pointer> images;
    for( unsigned int n = 4; n < static_cast<unsigned int>( argc ); n++ )
      {
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( argv[n] );
      reader->Update();
      images.push_back( reader->GetOutput() );
      }

    typename ImageType::Pointer output = ImageType::New();
    output->SetOrigin( images[0]->GetOrigin() );
    output->SetSpacing( images[0]->GetSpacing() );
    output->SetRegions( images[0]->GetLargestPossibleRegion() );
    output->SetDirection( images[0]->GetDirection() );
    output->Allocate();
    output->FillBuffer( 0 );

    itk::ImageRegionIteratorWithIndex<ImageType> ItO( output,
      output->GetLargestPossibleRegion() );
    for( ItO.GoToBegin(); !ItO.IsAtEnd(); ++ItO )
      {
      float exVent = 0.0;
      for( unsigned int i = 0; i < images.size(); i++ )
        {
        exVent += ( ( i + 1 ) * images[i]->GetPixel( ItO.GetIndex() ) );
        }
      ItO.Set( exVent );
      }
    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( output );
    writer->SetFileName( argv[3] );
    writer->Update();
    }
  else
    {
    std::cout << "Option not recognized." << std::endl;
    }

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if( argc < 6 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " imageDimension operation outputImage "
              << "image-list-via-wildcard " << std::endl;
    std::cerr << "  operations: " << std::endl;
    std::cerr << "    s:    Create speed image from atlas" << std::endl;
    std::cerr << "    mean: Create mean image" << std::endl;
    std::cerr << "    var:  Create variance image" << std::endl;
    std::cerr << "    w:    create probabilistic weight image from label probability images" << std::endl;
    std::cerr << "    ex:   Create expected ventilation from posterior prob. images" << std::endl;
    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     MultipleOperateImages<2>( argc, argv );
     break;
   case 3:
     MultipleOperateImages<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}




#include "itkArray.h"
#include "itkImageDuplicator.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileWriter.h"

#include <itksys/SystemTools.hxx>

#include <string>
#include <vector>

#include <fstream>

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
  std::string firstFile = std::string( argv[5] );
  std::vector<std::string> filenames;
  std::cout << "Using the following files ";

  if( firstFile.find( ".txt" ) != std::string::npos || firstFile.find( ".csv" ) != std::string::npos )
    {
    std::cout << "(reading images names from text file):" << std::endl;
    std::fstream str( firstFile.c_str(), std::ios::in );
    std::string file;

    while( str >> file )
      {
      bool isFile = itksys::SystemTools::FileExists( file.c_str(), true );
      if( isFile )
        {
        filenames.push_back( file );
        }
      else
        {
        std::cout << "File does not exist---" << file.c_str() << std::endl;
        }
      }
    }
  else
    {
    std::cout << "(reading images names from command line):" << std::endl;
    for( unsigned int n = 5; n < static_cast<unsigned int>( argc ); n++ )
      {
      filenames.push_back( std::string( argv[n] ) );
      }
    }

  for( unsigned int n = 0; n < filenames.size(); n++ )
    {
    std::cout << "   " << n+1 << ": " << filenames[n].c_str() << std::endl;
    }

  typename LabelImageType::Pointer mask = LabelImageType::New();
  mask = NULL;
  try
    {
    typedef itk::ImageFileReader<LabelImageType> ReaderType;
    typename ReaderType::Pointer labelImageReader = ReaderType::New();
    labelImageReader->SetFileName( argv[4] );
    labelImageReader->Update();
    mask = labelImageReader->GetOutput();
//    labelImageReader->DisconnectPipeline();
    }
  catch(...)
    {
    std::cout << "Not using a mask." << std::endl;
    }

  std::string op = std::string( argv[2] );

  if( op.compare( std::string( "mean" ) ) == 0 )
    {
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( filenames[0].c_str() );
    reader->Update();

    typename ImageType::Pointer output = reader->GetOutput();

    float N = 1.0;

    for( unsigned int n = 1; n < filenames.size(); n++ )
      {
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( filenames[n].c_str() );
      reader->Update();

      N += 1.0;
      itk::ImageRegionIteratorWithIndex<ImageType> It( reader->GetOutput(),
        reader->GetOutput()->GetLargestPossibleRegion() );
      itk::ImageRegionIterator<ImageType> ItO( output,
        output->GetLargestPossibleRegion() );
      for( It.GoToBegin(), ItO.GoToBegin(); !It.IsAtEnd(); ++It, ++ItO )
        {
        if( !mask || mask->GetPixel( It.GetIndex() ) != 0 )
          {
          ItO.Set( ItO.Get() * ( N - 1.0 ) / N + It.Get() / N );
          }
        }
      }

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( output );
    writer->SetFileName( argv[3] );
    writer->Update();
    }
  else if( op.compare( std::string( "sum" ) ) == 0 )
    {
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( filenames[0].c_str() );
    reader->Update();

    typename ImageType::Pointer output = reader->GetOutput();

    for( unsigned int n = 1; n < filenames.size(); n++ )
      {
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( filenames[n].c_str() );
      reader->Update();

      itk::ImageRegionIteratorWithIndex<ImageType> It( reader->GetOutput(),
        reader->GetOutput()->GetLargestPossibleRegion() );
      itk::ImageRegionIterator<ImageType> ItO( output,
        output->GetLargestPossibleRegion() );
      for( It.GoToBegin(), ItO.GoToBegin(); !It.IsAtEnd(); ++It, ++ItO )
        {
        if( !mask || mask->GetPixel( It.GetIndex() ) != 0 )
          {
          ItO.Set( ItO.Get() + It.Get() );
          }
        }
      }

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( output );
    writer->SetFileName( argv[3] );
    writer->Update();
    }
  else if( op.compare( std::string( "max" ) ) == 0 )
    {
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( filenames[0].c_str() );
    reader->Update();

    typename ImageType::Pointer output = reader->GetOutput();

    for( unsigned int n = 1; n < filenames.size(); n++ )
      {
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( filenames[n].c_str() );
      reader->Update();

      itk::ImageRegionIteratorWithIndex<ImageType> It( reader->GetOutput(),
        reader->GetOutput()->GetLargestPossibleRegion() );
      itk::ImageRegionIterator<ImageType> ItO( output,
        output->GetLargestPossibleRegion() );
      for( It.GoToBegin(), ItO.GoToBegin(); !It.IsAtEnd(); ++It, ++ItO )
        {
        if( !mask || mask->GetPixel( It.GetIndex() ) != 0 )
          {
          ItO.Set( vnl_math_max( ItO.Get(), It.Get() ) );
          }
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
    reader->SetFileName( filenames[0].c_str() );
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

    for( unsigned int n = 1; n < filenames.size(); n++ )
      {
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( filenames[n].c_str() );
      reader->Update();

      N += 1.0;

      itk::ImageRegionIteratorWithIndex<ImageType> It( reader->GetOutput(),
        reader->GetOutput()->GetLargestPossibleRegion() );
      itk::ImageRegionIterator<ImageType> ItO( output,
        output->GetLargestPossibleRegion() );
      itk::ImageRegionIterator<ImageType> ItM( meanImage,
        meanImage->GetLargestPossibleRegion() );
      for( It.GoToBegin(), ItO.GoToBegin(), ItM.GoToBegin(); !It.IsAtEnd(); ++It, ++ItO, ++ItM )
        {
        if( !mask || mask->GetPixel( It.GetIndex() ) != 0 )
          {
          ItM.Set( ItM.Get() * ( N - 1.0 ) / N + It.Get() / N );
          if ( N > 1.0 )
            {
            ItO.Set( ItO.Get() * ( N - 1.0 )/N +
              ( It.Get() - ItM.Get() )*( It.Get() - ItM.Get() ) / ( N - 1.0 ) );
            }
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
    for( unsigned int n = 0; n < filenames.size(); n++ )
      {
      typename LabelReaderType::Pointer reader = LabelReaderType::New();
      reader->SetFileName( filenames[n].c_str() );
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
      if( !mask || mask->GetPixel( It.GetIndex() ) != 0 )
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
    for( unsigned int n = 0; n < filenames.size(); n++ )
      {
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( filenames[n].c_str() );
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
      if( !mask || mask->GetPixel( ItO.GetIndex() ) != 0 )
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
      }
    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( output );
    writer->SetFileName( argv[3] );
    writer->Update();
    }
  else if( op.compare( std::string( "seg" ) ) == 0 )
    {
    std::vector<typename ImageType::Pointer> images;
    for( unsigned int n = 0; n < filenames.size(); n++ )
      {
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( filenames[n].c_str() );
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
      if( mask && mask->GetPixel( ItO.GetIndex() ) != 0 )
        {
        float maxProbability = 0;
        float maxLabel = 0;
        for( unsigned int i = 0; i < images.size(); i++ )
          {
          if( images[i]->GetPixel( ItO.GetIndex() ) >= maxProbability )
            {
            maxProbability = images[i]->GetPixel( ItO.GetIndex() );
            maxLabel = i + 1;
            }
          }
        ItO.Set( maxLabel );
        }
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
    for( unsigned int n = 0; n < filenames.size(); n++ )
      {
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( filenames[n].c_str() );
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
      if( !mask || mask->GetPixel( ItO.GetIndex() ) != 0 )
        {
        float exVent = 0.0;
        for( unsigned int i = 0; i < images.size(); i++ )
          {
          exVent += ( ( i + 1 ) * images[i]->GetPixel( ItO.GetIndex() ) );
          }
        ItO.Set( exVent );
        }
      }
    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( output );
    writer->SetFileName( argv[3] );
    writer->Update();
    }
  else if( op.compare( std::string( "sample" ) ) == 0 )
    {
    std::vector<typename ImageType::Pointer> images;
    for( unsigned int n = 0; n < filenames.size(); n++ )
      {
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( filenames[n].c_str() );
      reader->Update();
      images.push_back( reader->GetOutput() );
      }

    std::string sampleFilename = std::string( argv[3] ) +
      std::string( "Samples.txt" );

    std::string indexFilename = std::string( argv[3] ) +
      std::string( "Indices.txt" );

    std::ofstream str( sampleFilename.c_str() );
    std::ofstream str2( indexFilename.c_str() );

    for( unsigned int n = 0; n < filenames.size() - 1; n++ )
      {
      str << filenames[n] << ",";
      }
    str << filenames[filenames.size()-1] << std::endl;

    itk::ImageRegionIteratorWithIndex<ImageType> It( images[0],
      images[0]->GetLargestPossibleRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      if( !mask || mask->GetPixel( It.GetIndex() ) != 0 )
        {
        for( unsigned int d = 0; d < ImageDimension-1; d++ )
          {
          str2 << It.GetIndex()[d] << ",";
          }
        str2 << It.GetIndex()[ImageDimension-1] << std::endl;
        for( unsigned int n = 0; n < images.size()-1; n++ )
          {
          str << images[n]->GetPixel( It.GetIndex() ) << ",";
          }
        str << images[images.size()-1]->GetPixel( It.GetIndex() ) << std::endl;
        }
      }
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
    std::cerr << argv[0] << " imageDimension operation outputImage [maskImage] "
              << "image-list-via-wildcard " << std::endl;
    std::cerr << "  operations: " << std::endl;
    std::cerr << "    s:      Create speed image from atlas" << std::endl;
    std::cerr << "    mean:   Create mean image" << std::endl;
    std::cerr << "    sum:    Create sum image" << std::endl;
    std::cerr << "    max:    Create max image" << std::endl;
    std::cerr << "    var:    Create variance image" << std::endl;
    std::cerr << "    w:      create probabilistic weight image from label probability images" << std::endl;
    std::cerr << "    seg:    create labe image from label probability images" << std::endl;
    std::cerr << "    ex:     Create expected ventilation from posterior prob. images" << std::endl;
    std::cerr << "    sample: Print samples to output text/index files (prefix specified in place of outputImage)" << std::endl;
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




#include "itkAddImageFilter.h"
#include "itkArray.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileWriter.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkRescaleIntensityImageFilter.h"

#include <itksys/SystemTools.hxx>

#include "vnl/algo/vnl_fft_1d.h"
#include "vnl/vnl_complex_traits.h"
#include "vcl_complex.h"

#include <string>
#include <vector>
#include <sstream>

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
      std::istringstream iss2( element );
      iss2 >> value;
      values.push_back( value );
      }
    }
  return values;
}

typedef float RealType;

RealType CalculatePearsonCoefficient( std::vector<RealType> X, std::vector<RealType> Y )
{
  if( X.size() != Y.size() )
    {
    std::cerr << "Vectors are not the same size" << std::endl;
    exit( 1 );
    }

  RealType N = X.size();

  std::vector<RealType>::const_iterator itX;
  std::vector<RealType>::const_iterator itY;

  RealType sumX = 0.0;
  RealType sumY = 0.0;
  RealType sumX2 = 0.0;
  RealType sumXY = 0.0;
  RealType sumY2 = 0.0;

  for( itX = X.begin(), itY = Y.begin(); itX != X.end(); ++itX, ++itY )
    {
    sumX  += (*itX);
    sumY  += (*itY);
    sumXY += (*itX) * (*itY);
    sumX2 += (*itX) * (*itX);
    sumY2 += (*itY) * (*itY);
    }

  RealType r = ( N * sumXY - sumX*sumY ) / ( ( vcl_sqrt( N *sumX2 - (sumX*sumX) ) ) * ( vcl_sqrt( N *sumY2 - (sumY*sumY) ) ) );

  return r;
}

std::vector<RealType> FitRegressionLine(
  std::vector<RealType> X, std::vector<RealType> Y )
{
  if( X.size() != Y.size() )
    {
    std::cerr << "Vectors are not the same size" << std::endl;
    exit( 1 );
    }

  RealType N = X.size();

  std::vector<RealType>::const_iterator itX;
  std::vector<RealType>::const_iterator itY;

  RealType sumX = 0.0;
  RealType sumY = 0.0;
  RealType sumX2 = 0.0;
  RealType sumXY = 0.0;

  for( itX = X.begin(), itY = Y.begin(); itX != X.end(); ++itX, ++itY )
    {
    sumX  += (*itX);
    sumY  += (*itY);
    sumXY += (*itX) * (*itY);
    sumX2 += (*itX) * (*itX);
    }

  std::vector<RealType> line( 2 );
  line[0] = ( N * sumXY - sumX*sumY ) / ( N *sumX2 - sumX*sumX );
  if( sumX2 == 0 )
    {
    line[1] = sumY / N;
    }
  else
    {
    line[1] = ( sumY - line[0] * sumX ) / N;
    }

  return line;
}


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
      if( !mask || mask->GetPixel( ItO.GetIndex() ) != 0 )
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
        if( maxProbability > 0 )
          {
          ItO.Set( maxLabel );
          }
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
  else if( op.compare( std::string( "labelAvg" ) ) == 0 )
    {
    typedef unsigned int LabelType;

    std::vector<LabelType> labels;
    typedef itk::Image<LabelType, ImageDimension> LabelImageType;

    std::vector<typename LabelImageType::Pointer> images;
    std::vector<typename ImageType::Pointer> outputImages;
    std::vector<std::string> outputFilenames;
    for( unsigned int n = 0; n < filenames.size(); n++ )
      {
      typename LabelReaderType::Pointer reader = LabelReaderType::New();
      reader->SetFileName( filenames[n].c_str() );
      reader->Update();
      images.push_back( reader->GetOutput() );

      typedef itk::LabelGeometryImageFilter<LabelImageType, ImageType>  FilterType;
      typename FilterType::Pointer filter = FilterType::New();
      filter->SetInput( images[n] );
      filter->Update();

      std::vector<LabelType> imageLabels = filter->GetLabels();

      std::vector<LabelType>::const_iterator it;
      for( it = imageLabels.begin(); it != imageLabels.end(); ++it )
        {
        if( std::find( labels.begin(), labels.end(), *it ) == labels.end() && *it != 0 )
          {
          labels.push_back( *it );
          }
        }
      }

    for( unsigned int n = 0; n < labels.size(); n++ )
      {
      std::ostringstream str;
      str << n;

      std::string outname = std::string( argv[3] ) + std::string( "" ) + str.str() + std::string( ".nii.gz" );
      outputFilenames.push_back( outname );

      typename ImageType::Pointer outImage = ImageType::New();
      outImage->CopyInformation( images[0] );
      outImage->SetRegions( images[0]->GetLargestPossibleRegion() );
      outImage->Allocate();
      outImage->FillBuffer( 0 );

      outputImages.push_back( outImage );
      }

    RealType numberOfImages = static_cast<RealType>( filenames.size() );

    std::vector<LabelType>::const_iterator it;
    for( it = labels.begin(); it != labels.end(); ++it )
      {
      for( unsigned int n = 0; n < filenames.size(); n++ )
        {
        typedef itk::BinaryThresholdImageFilter<LabelImageType, ImageType> ThresholderType;
        typename ThresholderType::Pointer thresholder = ThresholderType::New();
        thresholder->SetInput( images[n] );
        thresholder->SetLowerThreshold( *it );
        thresholder->SetUpperThreshold( *it );
        thresholder->SetInsideValue( 1.0 );
        thresholder->SetOutsideValue( 0.0 );
        thresholder->Update();

        typedef itk::ImageDuplicator<ImageType> DuplicatorType;
        typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
        duplicator->SetInputImage( outputImages[it - labels.begin()] );
        duplicator->Update();

        typedef itk::AddImageFilter<ImageType, ImageType, ImageType> AdderType;
        typename AdderType::Pointer adder = AdderType::New();
        adder->SetInput1( thresholder->GetOutput() );
        adder->SetInput2( duplicator->GetOutput() );
        adder->Update();

        outputImages[it - labels.begin()] = adder->GetOutput();
        outputImages[it - labels.begin()]->DisconnectPipeline();
        }
      itk::ImageRegionIterator<ImageType> It( outputImages[it - labels.begin()], outputImages[it - labels.begin()]->GetLargestPossibleRegion() );
      for( It.GoToBegin(); !It.IsAtEnd(); ++It )
        {
        It.Set( It.Get() / numberOfImages );
        }
      }

    for( unsigned int n = 0; n < labels.size(); n++ )
      {
      typedef itk::ImageFileWriter<ImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( outputImages[n] );
      writer->SetFileName( outputFilenames[n] );
      writer->Update();
      }
    }
  else if( op.compare( std::string( "fft" ) ) == 0 )
    {
    std::vector<typename ImageType::Pointer> images;
    std::vector<typename ImageType::Pointer> outputImages;
    std::vector<std::string> outputFilenames;
    for( unsigned int n = 0; n < filenames.size(); n++ )
      {
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( filenames[n].c_str() );
      reader->Update();
      images.push_back( reader->GetOutput() );
      }

    unsigned int numberOfImages = filenames.size();

    RealType exponent = vcl_ceil( vcl_log( static_cast<RealType>( numberOfImages ) ) / vcl_log( 2.0 ) );
    unsigned int paddedSize = static_cast<unsigned int>( vcl_pow( static_cast<RealType>( 2.0 ), exponent ) + 0.5 );

    for( unsigned int n = 0; n < paddedSize; n++ )
      {
      std::string leadingZeros = std::string( 4, '0' );

      if( n > 0 )
        {
        std::ostringstream str;
        str << n;
        leadingZeros = std::string( 4 - static_cast<unsigned int>( vcl_log( n )/vcl_log( 10 ) + 1 ), '0' ).append( str.str() );
        }

      std::string outname = std::string( argv[3] ) + std::string( "FT" ) + leadingZeros + std::string( ".nii.gz" );
      outputFilenames.push_back( outname );

      typename ImageType::Pointer outImage = ImageType::New();
      outImage->CopyInformation( images[0] );
      outImage->SetRegions( images[0]->GetLargestPossibleRegion() );
      outImage->Allocate();
      outImage->FillBuffer( 0 );

      outputImages.push_back( outImage );
      }

    vnl_fft_1d<RealType> fft( paddedSize );

    itk::ImageRegionIteratorWithIndex<ImageType> It( images[0],
      images[0]->GetLargestPossibleRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      if( !mask || mask->GetPixel( It.GetIndex() ) != 0 )
        {
        vnl_vector< vcl_complex<RealType> > V( paddedSize, vcl_complex<RealType>( 0.0, 0.0 ) );

        for( unsigned int n = 0; n < numberOfImages; n++ )
          {
          // Multiply intensity by Hann window
          V[n] = images[n]->GetPixel( It.GetIndex() ) * 0.5 * ( 1 - vcl_cos( 2 * vnl_math::pi * n / ( paddedSize - 1 ) ) );
          }
        fft.fwd_transform( V );

        for( unsigned int n = 0; n < paddedSize; n++ )
          {
          outputImages[n]->SetPixel( It.GetIndex(),  vcl_norm( V[n] ) );
          }
        }
      else
        {
        It.Set( 0 );
        }
      }

    for( unsigned int n = 0; n < paddedSize; n++ )
      {
      typedef itk::ImageFileWriter<ImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( outputImages[n] );
      writer->SetFileName( outputFilenames[n] );
      writer->Update();
      }
    }
  else if( op.compare( 0, 4, std::string( "corr", 0, 4 ) ) == 0 )
    {
    std::vector<typename ImageType::Pointer> images;
    for( unsigned int n = 0; n < filenames.size(); n++ )
      {
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( filenames[n].c_str() );
      reader->Update();
      images.push_back( reader->GetOutput() );
      }

    std::string vectorString = op.substr( 5 );

    std::vector<RealType> corrVector = ConvertVector<RealType>( vectorString );

    if( corrVector.size() != filenames.size() )
      {
      std::cerr << "Error: the size of the specified correlation vector does not equal the number of images." << std::endl;
      return EXIT_FAILURE;
      }

    itk::ImageRegionIteratorWithIndex<ImageType> It( images[0],
      images[0]->GetLargestPossibleRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      if( !mask || mask->GetPixel( It.GetIndex() ) != 0 )
        {
        std::vector<RealType> intensities;
        for( unsigned int n = 0; n < filenames.size(); n++ )
          {
          intensities.push_back( images[n]->GetPixel( It.GetIndex() ) );
          }

        RealType r = CalculatePearsonCoefficient( corrVector, intensities );
        It.Set( r );
        }
      else
        {
        It.Set( 0 );
        }
      }

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( images[0] );
    writer->SetFileName( argv[3] );
    writer->Update();
    }
  else if( op.compare( 0, 5, std::string( "slope", 0, 5 ) ) == 0 )
    {
    std::vector<typename ImageType::Pointer> images;
    for( unsigned int n = 0; n < filenames.size(); n++ )
      {
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( filenames[n].c_str() );
      reader->Update();
      images.push_back( reader->GetOutput() );
      }

    std::string vectorString = op.substr( 6 );

    std::vector<RealType> corrVector = ConvertVector<RealType>( vectorString );

    std::vector<RealType> line;

    if( corrVector.size() != filenames.size() )
      {
      std::cerr << "Error: the size of the specified correlation vector does not equal the number of images." << std::endl;
      return EXIT_FAILURE;
      }

    itk::ImageRegionIteratorWithIndex<ImageType> It( images[0],
      images[0]->GetLargestPossibleRegion() );
    for( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      if( !mask || mask->GetPixel( It.GetIndex() ) != 0 )
        {
        std::vector<RealType> intensities;
        for( unsigned int n = 0; n < filenames.size(); n++ )
          {
          intensities.push_back( images[n]->GetPixel( It.GetIndex() ) );
          }

        line = FitRegressionLine( corrVector, intensities );
        It.Set( line[0] );
        }
      else
        {
        It.Set( 0 );
        }
      }

    typedef itk::ImageFileWriter<ImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetInput( images[0] );
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
      std::string( "Samples.csv" );

    std::string indexFilename = std::string( argv[3] ) +
      std::string( "Indices.csv" );

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
  else if( op.compare( 0, 6, std::string( "cohort", 0, 6 ) ) == 0 )
    {
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( filenames[0].c_str() );
    reader->Update();

    std::string numberString = op.substr( 7 );

    unsigned int numberOfSubjects = Convert<unsigned int>( numberString );

    typename ImageType::Pointer meanImage = reader->GetOutput();

    float N = 1.0;

    typename ImageType::Pointer variance = ImageType::New();
    variance->SetOrigin( meanImage->GetOrigin() );
    variance->SetSpacing( meanImage->GetSpacing() );
    variance->SetRegions( meanImage->GetLargestPossibleRegion() );
    variance->SetDirection( meanImage->GetDirection() );
    variance->Allocate();
    variance->FillBuffer( 0 );

    for( unsigned int n = 1; n < filenames.size(); n++ )
      {
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( filenames[n].c_str() );
      reader->Update();

      N += 1.0;

      itk::ImageRegionIteratorWithIndex<ImageType> It( reader->GetOutput(),
        reader->GetOutput()->GetLargestPossibleRegion() );
      itk::ImageRegionIterator<ImageType> ItO( variance,
        variance->GetLargestPossibleRegion() );
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

    typedef typename itk::Statistics::MersenneTwisterRandomVariateGenerator RandomizerType;
    typename RandomizerType::Pointer randomizer = RandomizerType::New();
    randomizer->Initialize();

    for( unsigned int n = 0; n < numberOfSubjects; n++ )
      {
      typename ImageType::Pointer output = ImageType::New();
      output->SetOrigin( meanImage->GetOrigin() );
      output->SetSpacing( meanImage->GetSpacing() );
      output->SetRegions( meanImage->GetLargestPossibleRegion() );
      output->SetDirection( meanImage->GetDirection() );
      output->Allocate();
      output->FillBuffer( 0 );

      std::cout << "Creating subject " << n << std::endl;

      itk::ImageRegionIterator<ImageType> ItV( variance,
        variance->GetLargestPossibleRegion() );
      itk::ImageRegionIterator<ImageType> ItM( meanImage,
        meanImage->GetLargestPossibleRegion() );
      itk::ImageRegionIteratorWithIndex<ImageType> ItO( output,
        output->GetLargestPossibleRegion() );

      for( ItO.GoToBegin(), ItM.GoToBegin(), ItV.GoToBegin(); !ItM.IsAtEnd(); ++ItO, ++ItM, ++ItV )
        {
        if( !mask || mask->GetPixel( ItO.GetIndex() ) != 0 )
          {
          double intensity = randomizer->GetNormalVariate( ItM.Get(), ItV.Get() );
          intensity = vnl_math_min( 1.0, intensity );
          intensity = vnl_math_max( 0.0, intensity );
          ItO.Set( intensity );
          }
        }
      std::ostringstream str;
      str << n;

      std::string subjectFilename = std::string( argv[3] ) + std::string( "Subject" ) +
         str.str() + std::string( ".nii.gz" );

      typedef itk::ImageFileWriter<ImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( output );
      writer->SetFileName( subjectFilename.c_str() );
      writer->Update();
      }
    }
  else if( op.compare( 0, 9, std::string( "normalize", 0, 9 ) ) == 0 )
    {
    std::string paramString = op.substr( 10 );
    std::vector<RealType> params = ConvertVector<RealType>( paramString );

    if( params.size() != 3 )
      {
      std::cerr << "Incorrect number of parameters.  Please see option for \"normalize\"." << std::endl;
      return EXIT_FAILURE;
      }

    RealType p1 = params[0];
    RealType p2 = params[1];
    unsigned int nInnerQuantiles = static_cast<unsigned int>( params[2] );

    if( !mask )
      {
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( filenames[0].c_str() );
      reader->Update();

      mask = LabelImageType::New();
      mask->CopyInformation( reader->GetOutput() );
      mask->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
      mask->Allocate();
      mask->FillBuffer( 1 );
      }

    vnl_vector<RealType> p1_n( filenames.size(), 0.0 );
    vnl_vector<RealType> p2_n( filenames.size(), 0.0 );
    std::vector<vnl_vector<RealType> > IQs_n;
    IQs_n.resize( filenames.size() );

    vnl_vector<RealType> avgIQ( nInnerQuantiles, 0.0 );

    // Calculate the average quantile values from all the images.

    for( unsigned int n = 0; n < filenames.size(); n++ )
      {
      IQs_n[n].set_size( nInnerQuantiles );
      IQs_n[n].fill( 0.0 );

      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( filenames[n].c_str() );
      reader->Update();

      RealType maxValue = itk::NumericTraits<RealType>::min();
      RealType minValue = itk::NumericTraits<RealType>::max();

      typename ImageType::Pointer image = reader->GetOutput();
      itk::ImageRegionConstIteratorWithIndex<ImageType> ItI( image, image->GetRequestedRegion() );
      for( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
        {
        if( mask->GetPixel( ItI.GetIndex() ) == 1 )
          {
          if( ItI.Get() < minValue )
            {
            minValue = ItI.Get();
            }
          else if( ItI.Get() > maxValue )
            {
            maxValue = ItI.Get();
            }
          }
        }
      typedef itk::LabelStatisticsImageFilter<ImageType, LabelImageType> HistogramGeneratorType;
      typename HistogramGeneratorType::Pointer stats = HistogramGeneratorType::New();
      stats->SetInput( reader->GetOutput() );
      stats->SetLabelInput( mask );
      stats->SetUseHistograms( true );
      stats->SetHistogramParameters( 255, minValue, maxValue );
      stats->Update();

      typedef typename HistogramGeneratorType::HistogramType  HistogramType;
      const HistogramType *histogram = stats->GetHistogram( 1 );

      p1_n[n] = histogram->Quantile( 0, p1 );
      p2_n[n] = histogram->Quantile( 0, p2 );

      RealType deltaQ = static_cast<RealType>( ( p2 - p1 ) / ( nInnerQuantiles + 2 - 1 ) );
      for( unsigned int q = 0; q < nInnerQuantiles; q++ )
        {
        double quantileValue = histogram->Quantile( 0, p1 + ( q + 1 ) * deltaQ );
        avgIQ[q] += quantileValue;
        IQs_n[n][q] = quantileValue;
        }
      }

    // Renormalize images based on quantile values

    RealType avgP1 = p1_n.mean();
    RealType avgP2 = p2_n.mean();
    for( unsigned int q = 0; q < nInnerQuantiles; q++ )
      {
      avgIQ[q] /= static_cast<RealType>( filenames.size() );
      }

    for( unsigned int n = 0; n < filenames.size(); n++ )
      {
      typename ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( filenames[n].c_str() );
      reader->Update();

      typename ImageType::Pointer output = ImageType::New();
      output->CopyInformation( reader->GetOutput() );
      output->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
      output->Allocate();
      output->FillBuffer( 0 );

      typename ImageType::Pointer image = reader->GetOutput();
      itk::ImageRegionConstIteratorWithIndex<ImageType> ItI( image, image->GetRequestedRegion() );
      for( ItI.GoToBegin(); !ItI.IsAtEnd(); ++ItI )
        {
        if( mask->GetPixel( ItI.GetIndex() ) == 1 )
          {
          RealType x = vnl_math_max( vnl_math_min( ItI.Get(), p2_n[n] ), p1_n[n] );

          RealType x1 = p1_n[n];
          RealType x2 = IQs_n[n][0];
          RealType y1 = avgP1;
          RealType y2 = avgIQ[0];
          if( x > IQs_n[n][0] )
            {
            if( x >= IQs_n[n][nInnerQuantiles-1] )
              {
              x1 = IQs_n[n][nInnerQuantiles-1];
              x2 = p2_n[n];
              y1 = avgIQ[nInnerQuantiles-1];
              y2 = avgP2;
              }
            else
              {
              for( unsigned int q = 0; q < nInnerQuantiles-1; q++ )
                {
                if( x >= IQs_n[n][q] && x < IQs_n[n][q+1] )
                  {
                  x1 = IQs_n[n][q];
                  x2 = IQs_n[n][q+1];
                  y1 = avgIQ[q];
                  y2 = avgIQ[q+1];
                  }
                }
              }
            }
          RealType y = ( y2 - y1 ) / ( x2 - x1 ) * ( x - x1 ) + y1;
          output->SetPixel( ItI.GetIndex(), y );
          }
        }

      std::string path = itksys::SystemTools::GetFilenamePath( argv[3] );
      std::string filenameRoot = itksys::SystemTools::GetFilenameName( filenames[n] );
      std::string filenamePrefix = itksys::SystemTools::GetFilenameName( argv[3] );

      std::string outputName = path + std::string( "/" ) + filenamePrefix + filenameRoot;

      typedef itk::ImageFileWriter<ImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( output );
      writer->SetFileName( outputName.c_str() );
      writer->Update();
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
    std::cerr << argv[0] << " imageDimension operation outputImage [mask] "
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
    std::cerr << "    fft:    Perform voxelwise fft" << std::endl;
    std::cerr << "    labelAvg:   Perform labelwise averaging." << std::endl;
    std::cerr << "    corr=mxnxoxp...:   Create voxelwise correlation map with vector <m,n,x,o,p>" << std::endl;
    std::cerr << "    slope=mxnxoxp...:   Create voxelwise regression slope map with vector <m,n,x,o,p>" << std::endl;
    std::cerr << "    cohort=n...:   Create random cohort of n subjects from sample (gaussian modeling)" << std::endl;
    std::cerr << "    normalize=p1xp2xn:   Create normalized image set.  0 <= p1 < p2 <= 1.0 percentile min/max input intensity" << std::endl;
    std::cerr << "                         n is number of inner quantiles. " << std::endl;
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




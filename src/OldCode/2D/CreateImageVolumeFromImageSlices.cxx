#include <stdio.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNumericSeriesFileNames.h"
#include "itkMatrix.h"

#include "global.h"
#include "fstream.h"


int main( unsigned int argc, char *argv[] )
{
  if ( argc != 5 )
    {
    std::cout << "Usage: " << argv[0] << " format startIndex endIndex outputImage" << std::endl;
    exit( 1 );
    }
  
//  int numberOfSlices = atoi( argv[3] )-atoi( argv[2] )+1;

  typedef itk::NumericSeriesFileNames    NameGeneratorType;
  NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
  nameGenerator->SetSeriesFormat( argv[1] );
  nameGenerator->SetStartIndex( atoi( argv[2] ) );
  nameGenerator->SetEndIndex( atoi( argv[3] ) );
  nameGenerator->SetIncrementIndex( 1 );
  
  typedef itk::Image<PixelType, 2> Image2DType;
  typedef itk::ImageFileReader<Image2DType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();

  for ( int i = 0; i < atoi( argv[3] )-atoi( argv[2] )+1; i++ )
    {
    std::cout << "i = " << i << ", name = " << nameGenerator->GetFileNames()[i] << std::endl;

    reader->SetFileName( nameGenerator->GetFileNames()[i] );
    try 
      {
      reader->Update(); 
      } 
    catch( itk::ExceptionObject & err ) 
      { 
      std::cout << "ERROR: "<< err << std::endl; 
      return -1; 
      } 
    }
  return 0;
}


/*

int main( int argc, char *argv[] )
{
  if ( argc != 5 )
    {
    std::cout << "Usage: " << argv[0] << " format startIndex endIndex outputImage" << std::endl;
    exit( 1 );
    }
  
  int numberOfSlices = atoi( argv[3] )-atoi( argv[2] )+1;

  typedef itk::Image<PixelType, 3> Image3DType;
  Image3DType::Pointer image = Image3DType::New();
  Image3DType::SizeType size;

  typedef itk::Image<PixelType, 2> Image2DType;

  typedef itk::NumericSeriesFileNames    NameGeneratorType;
  NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();
  nameGenerator->SetSeriesFormat( argv[1] );
  nameGenerator->SetStartIndex( atoi( argv[2] ) );
  nameGenerator->SetEndIndex( atoi( argv[3] ) );
  nameGenerator->SetIncrementIndex( 1 );
  
  typedef itk::Image<PixelType, 2> Image2DType;
  typedef itk::ImageFileReader<Image2DType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( nameGenerator->GetFileNames()[0] );
  reader->Update();
  
  size[0] = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[0];
  size[1] = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[1];
  size[2] = numberOfSlices;
  
  image->SetRegions( size );
  image->Allocate();  
  image->FillBuffer( 0 );
  Image3DType::SpacingType spacing;
  spacing[0] = reader->GetOutput()->GetSpacing()[0];
  spacing[1] = reader->GetOutput()->GetSpacing()[0];
  spacing[2] = reader->GetOutput()->GetSpacing()[0];

  image->SetSpacing( spacing );

  for ( unsigned int i = 0; i < atoi( argv[3] )-atoi( argv[2] )+1; i++ )
    {
    
    std::cout << "i = " << i << ", line = " << nameGenerator->GetFileNames()[i] << std::endl;

    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( nameGenerator->GetFileNames()[i] );
    try 
      {
      reader->Update(); 
      } 
    catch( itk::ExceptionObject & err ) 
      { 
      std::cout << "ERROR: "<< err << std::endl; 
      return -1; 
      } 
    itk::ImageRegionIteratorWithIndex<Image2DType> It
      ( reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );
    for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
      {
      Image3DType::IndexType idx;
      idx[0] = It.GetIndex()[0];
      idx[1] = It.GetIndex()[1];
      idx[2] = i;
      image->SetPixel( idx, It.Get() ); 
      }
    }
  itk::Matrix<double> direction;
  direction.SetIdentity();

  image->SetDirection( direction );

  typedef itk::ImageFileWriter<Image3DType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( image );
  writer->SetFileName( argv[4] );
  writer->Update(); 
  return 0;
}
*/ 
#include "itkBinaryThinning3DImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkPathIterator.h"
#include "itkLiveWireImageFunction.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkTimeProbe.h"

#include "global.h"

int main( unsigned int argc, char *argv[] )
{
  if ( argc < 2*ImageDimension + 3 )
    {
    std::cout << "Usage: " << argv[0] << " inputImage outputImage "
              << "index1[0] .. index1[n] index2[0] .. index2[n] "
              << "[gradientMagnitudeWeight] [zeroCrossingWeight] "
              << "[gradientDirectionWeight] [findStepEdges]" << std::endl;     
    exit( 0 );
    }   

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<unsigned int, ImageDimension> LabelImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  LabelImageType::Pointer labelImage = LabelImageType::New();
  labelImage->SetOrigin( reader->GetOutput()->GetOrigin() );
  labelImage->SetSpacing( reader->GetOutput()->GetSpacing() );
  labelImage->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  labelImage->Allocate();
  labelImage->FillBuffer( 0 );

  typedef itk::LiveWireImageFunction<ImageType> FunctionType;
  FunctionType::Pointer function = FunctionType::New();

  if ( argc > 2*ImageDimension + 3 )
    {
    function->SetGradientMagnitudeWeight( atof( argv[2*ImageDimension + 3] ) );
    } 
  if ( argc > 2*ImageDimension + 4 )
    {
    function->SetZeroCrossingWeight( atof( argv[2*ImageDimension + 4] ) );
    } 
  if ( argc > 2*ImageDimension + 5 )
    {
    function->SetGradientDirectionWeight( atof( argv[2*ImageDimension + 5] ) );
    } 
  if ( argc > 2*ImageDimension + 6 )
    {
    function->SetFindStepEdges( atoi( argv[2*ImageDimension + 6] ) );
    } 
 
  if ( function->GetFindStepEdges() == false )
    {
    typedef itk::BinaryThinning3DImageFilter<ImageType, FunctionType::RealImageType> ThinnerType;
    ThinnerType::Pointer thinner = ThinnerType::New();
    thinner->SetInput( reader->GetOutput() );
    thinner->Update();
    function->SetZeroCrossingImage( thinner->GetOutput() );

    typedef itk::SignedMaurerDistanceMapImageFilter<ImageType, ImageType> FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput( reader->GetOutput() );
    filter->SetSquaredDistance( false );
    filter->SetUseImageSpacing( false );
    filter->SetInsideIsPositive( true );
    filter->Update();

    function->SetInputImage( filter->GetOutput() );  
    }  
  else
   {
   function->SetInputImage( reader->GetOutput() );
   }  

  ImageType::IndexType anchor;
  ImageType::IndexType free;

  for ( unsigned int d = 0; d < ImageDimension; d++ )
    {
    anchor[d] = atoi( argv[d+3] );
    free[d] = atoi( argv[d+3+ImageDimension] );
    } 

  function->SetAnchorSeed( anchor );

  FunctionType::OutputType::Pointer path = function->EvaluateAtIndex( free );
  
  typedef itk::PathIterator<LabelImageType, FunctionType::OutputType> IteratorType;
  typedef FunctionType::OutputType::ContinuousIndexType VertexType;
   
  IteratorType It( labelImage, path );
  It.GoToBegin();
  while ( !It.IsAtEnd() )
    {
    labelImage->SetPixel( It.GetIndex(), 1 );
    ++It;
    }

  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( labelImage );
  writer->Update();

  {
  typedef itk::ImageFileWriter<FunctionType::RealImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( "pathcost.nii.gz" );
  writer->SetInput( function->GetPathDirectionCostImage() );
  writer->Update();
  }

}

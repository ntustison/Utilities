#include "itkBoykovGraphTraits.h"
#include "itkGraph.h"
#include "itkBoykovImageToGraphFunctor.h"
#include "itkBoykovAlphaExpansionMRFImageFilter.h"

#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "vnl/vnl_math.h"
#include <vector>
#include <string>

#include "global.h"

#define pi 3.1415926535897932384626433832795

int main( int argc, char * argv[])
{

  if ( argc < 6 )
    {
    std::cout << argv[0] << " inputImage outputPrefix numberOfLabels likelihoodImage1 likelihoodImage2 ... likelihoodImageN" << std::endl;
    exit( 0 );
    } 
  

  typedef itk::Image<int,  ImageDimension> ImageType;
  typedef ImageType::IndexType IndexType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  typedef itk::Image<double, ImageDimension> LikelihoodImageType;
  typedef itk::ImageRegionIteratorWithIndex<LikelihoodImageType> IteratorType;
  typedef itk::CastImageFilter<ImageType, LikelihoodImageType> CasterType;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput( reader->GetOutput() );
  caster->Update();

  typedef itk::BoykovGraphTraits<double, ImageDimension>                  GraphTraitsType;  
  typedef itk::Graph<GraphTraitsType>                                     GraphType;  
  typedef itk::BoykovImageToGraphFunctor<LikelihoodImageType, GraphType>  FunctorType;  
  FunctorType::Pointer BoykovFunctor = FunctorType::New();
  BoykovFunctor->SetRadius( 1 );
  BoykovFunctor->SetLambda( 10.0 );
  BoykovFunctor->SetSigma( 10.0 );
  BoykovFunctor->SetExcludeBackground( false );
  BoykovFunctor->SetBackgroundValue( 0 );  

  typedef itk::BoykovAlphaExpansionMRFImageFilter
    <LikelihoodImageType, GraphTraitsType> FilterType;
  FilterType::Pointer filter = FilterType::New(); 
  filter->SetRandomizeInitialLabeling(false);
  filter->SetNumberOfClasses( atoi( argv[3] ) );
  filter->SetInput( caster->GetOutput() );
  filter->SetImageToGraphFunctor( BoykovFunctor );  

  for ( unsigned int i = 0; i < filter->GetNumberOfClasses(); i++ )
    {
    typedef itk::ImageFileReader<LikelihoodImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[4+i] );
    reader->Update();

    filter->SetLikelihoodImage( i+1, reader->GetOutput() );        
    } 
  filter->Update();
  
  std::string filename = std::string( argv[2] ) + std::string( "Labels.nii" );

  typedef itk::Image<int,  ImageDimension> LabelImageType;  
  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( filename.c_str() );  
  writer->SetInput(filter->GetOutput());
  writer->Write();
  
  std::cout << "Output files: out.hdr out.img" << std::endl;
  
  for ( unsigned int i = 0; i < filter->GetNumberOfClasses(); i++ )
    {  
    itk::OStringStream buf;
    buf << i;

    std::string filename = std::string( argv[2] ) + std::string( "PosteriorProbability_" )
      + buf.str() + std::string( ".nii" );

    typedef itk::ImageFileWriter<LikelihoodImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( filename.c_str() );  
    writer->SetInput( filter->GetPosteriorProbabilityImage( i ) );
    writer->Update();
    }
  
  return 0;  
};

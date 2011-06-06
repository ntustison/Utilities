#include "itkBoykovGraphTraits.h"
#include "itkBoykovImageToGraphFunctor.h"
#include "itkBoykovAlphaExpansionMRFImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkGraph.h"

#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkNumericTraits.h"
#include "vnl/vnl_math.h"
#include <vector>
#include <string>

template <unsigned int ImageDimension>
int BoykovGraphCutFilter( int argc, char * argv[])
{
  typedef float RealType;
 
  typedef itk::Image<RealType,  ImageDimension> ImageType;
  typedef typename ImageType::IndexType IndexType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::Image<int,  ImageDimension> MaskImageType;
  typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
  typename MaskReaderType::Pointer maskreader = MaskReaderType::New();
  maskreader->SetFileName( argv[3] );
  maskreader->Update();

  itk::ImageRegionIterator<ImageType> It( reader->GetOutput(),
    reader->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<MaskImageType> ItM( maskreader->GetOutput(),
    maskreader->GetOutput()->GetLargestPossibleRegion() );
  
  It.GoToBegin();
  ItM.GoToBegin();
  while( !ItM.IsAtEnd() )
    {
    if( ItM.Get() == 0 )
      {
      It.Set( itk::NumericTraits<RealType>::max() );
      }
    ++It;
    ++ItM;  
    }  

  typedef itk::Image<double, ImageDimension> LikelihoodImageType;
  typedef itk::ImageRegionIteratorWithIndex<LikelihoodImageType> IteratorType;
  typedef itk::CastImageFilter<ImageType, LikelihoodImageType> CasterType;
  typename CasterType::Pointer caster = CasterType::New();
  caster->SetInput( reader->GetOutput() );
  caster->Update();

  typedef itk::BoykovGraphTraits
    <double, ImageDimension> GraphTraitsType;  
  typedef itk::Graph<GraphTraitsType> GraphType;  
  typedef itk::BoykovImageToGraphFunctor
    <LikelihoodImageType, GraphType> FunctorType;  
  typename FunctorType::Pointer BoykovFunctor = FunctorType::New();
  BoykovFunctor->SetRadius( 1 );
  BoykovFunctor->SetSigma( 10.0 );
  BoykovFunctor->SetLambda( 10.0 );
  if( argc > 6+atoi( argv[5] ) && atoi( argv[6+atoi( argv[5] )] ) > 1 )
    {
    BoykovFunctor->SetRadius( atoi( argv[6+atoi( argv[5] )] ) ); 
    }
  BoykovFunctor->ActivateAllNeighbors();
  if( argc > 7+atoi( argv[5] ) )
    {
    BoykovFunctor->SetSigma( atof( argv[7+atoi( argv[5] )] ) );
    }
  if( argc > 8+atoi( argv[5] ) )
    {
    BoykovFunctor->SetLambda( atof( argv[8+atoi( argv[5] )] ) );
    }
  BoykovFunctor->SetExcludeBackground( true );
  BoykovFunctor->SetBackgroundValue( itk::NumericTraits<RealType>::max() ); 

  typedef itk::BoykovAlphaExpansionMRFImageFilter
    <LikelihoodImageType, GraphTraitsType> FilterType;
  typename FilterType::Pointer filter = FilterType::New(); 
  filter->SetRandomizeInitialLabeling( false );
  filter->SetNumberOfClasses( atoi( argv[5] ) );
  filter->SetInput( caster->GetOutput() );
  filter->SetImageToGraphFunctor( BoykovFunctor );  

  for ( unsigned int i = 0; i < filter->GetNumberOfClasses(); i++ )
    {
    typedef itk::ImageFileReader<LikelihoodImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[6+i] );
    reader->Update();

    filter->SetLikelihoodImage( i+1, reader->GetOutput() );        
    }
  filter->DebugOn();  
  filter->Update();
  
  std::string filename = std::string( argv[4] ) 
    + std::string( "Labels.nii.gz" );

  typedef itk::Image<int,  ImageDimension> LabelImageType;  
  typedef itk::ImageFileWriter<LabelImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( filename.c_str() );  
  writer->SetInput( filter->GetOutput() );
  writer->Write();
  
  for ( unsigned int i = 1; i <= filter->GetNumberOfClasses(); i++ )
    {  
    itk::OStringStream buf;
    buf << i;

    std::string filename = std::string( argv[4] ) 
      + std::string( "PosteriorProbability_" )
      + buf.str() + std::string( ".nii.gz" );
  
    typedef itk::ImageFileWriter<LikelihoodImageType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( filename.c_str() );  
    writer->SetInput( filter->GetPosteriorProbabilityImage( i ) );
    writer->Update();
    }
  
  return EXIT_SUCCESS;  
};

int main( int argc, char *argv[] )
{
  if ( argc < 6 + atoi( argv[5] ) )
    {
    std::cout << argv[0] << " ImageDimension inputImage maskImage outputPrefix "
      << "numberOfLabels(>=2) likelihoodImage1 ... likelihoodImageN "
      << "[radius] [sigma] [lambda] " 
      << std::endl;
    exit( 0 );
    } 

  switch( atoi( argv[1] ) ) 
   {
   case 2:
     BoykovGraphCutFilter<2>( argc, argv );
     break;
   case 3:
     BoykovGraphCutFilter<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBranchDecompositionFilter.h"
#include "itkCastImageFilter.h"
#include "itkJunctionDetectionFilter.h"


template<class TFilter>
class CommandProgressUpdate : public itk::Command 
{
public:
  typedef  CommandProgressUpdate    Self;
  typedef  itk::Command             Superclass;
  typedef  itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );
  
  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    const TFilter * filter = dynamic_cast< const TFilter * >(object);
    if( !(itk::ProgressEvent().CheckEvent( &event )) )  return;
    std::cout << (int)(100 * filter->GetProgress()) << "% completed\r" << std::flush;
  }

protected:
  CommandProgressUpdate(){};
};


template <unsigned int ImageDimension>
int JunctionAndBranch( int argc, char *argv[] )
{
  typedef unsigned short PixelType;

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  float inner = ( argc > 4 ) ? 2.0 : atof( argv[4] );
  float outer = ( argc > 5 ) ? 3.0 : atof( argv[5] );


  typedef itk::JunctionDetectionFilter<ImageType> DetectorType;
  typedef itk::BranchDecompositionFilter<ImageType> DecomposerType;




  typedef CommandProgressUpdate<DetectorType> CommandDetectorType;
  typename CommandDetectorType::Pointer observer1 = CommandDetectorType::New();
  typename DetectorType::Pointer detector = DetectorType::New();
  detector->SetInnerRadius( inner);
  detector->SetOuterRadius( outer );
  detector->SetMinNumberOfPixel( 8 );
  detector->SetInput( reader->GetOutput() );
  detector->AddObserver( itk::ProgressEvent(), observer1 );
  detector->Update();

  const typename DetectorType::JCLabelMapType& jcLabelmap = detector->GetJCLabelMap();

  typename DecomposerType::JCLabelMapType jcLabelMapInput;

  for( typename DetectorType::JCLabelMapType::const_iterator 
    jclmIt = jcLabelmap.begin(); jclmIt != jcLabelmap.end(); ++jclmIt )
    {
    int jcLabel = (*jclmIt).first;
    typename DetectorType::JCLabelPairType jclPair = (*jclmIt).second;
    jcLabelMapInput[jcLabel] = jclPair;
    }

  std::cout << "Done with junction detection" << std::endl;

  {
  std::string filename = std::string( argv[3] ) 
    + std::string( "Junctions.nii.gz" ); 
   
  typedef itk::ImageFileWriter<typename DetectorType::OutputImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( detector->GetOutput() );
  writer->SetFileName( filename.c_str() );                                          
  writer->Update();
  }

  typedef itk::CastImageFilter<typename DetectorType::OutputImageType, 
    ImageType> CasterType;
  typename CasterType::Pointer caster = CasterType::New();
  caster->SetInput( detector->GetOutput() );
  caster->Update();

  typename DecomposerType::Pointer decomposer = DecomposerType::New();
  typedef CommandProgressUpdate<DecomposerType> CommandDecomposerType;
  typename CommandDecomposerType::Pointer observer2 = CommandDecomposerType::New();
  decomposer->SetInnerRadius( inner );
  decomposer->SetOuterRadius( outer );
  decomposer->SetMinNumberOfPixel( 8 );
  decomposer->SetJCLabelMap( &jcLabelMapInput );
  decomposer->SetInput( caster->GetOutput() );
  decomposer->AddObserver( itk::ProgressEvent(), observer2 );
  decomposer->Update();

  std::cout << "Done with branch decomposition" << std::endl;
  
  {
  std::string filename = std::string( argv[3] ) 
    + std::string( "Branches.nii.gz" ); 
   
  typedef itk::ImageFileWriter<typename DecomposerType::OutputImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( decomposer->GetOutput() );
  writer->SetFileName( filename.c_str() );                                          
  writer->Update();
  }
    
  return 0;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << "Usage: " << argv[0] << " imageDimension inputImage "
      << "outputPrefix [innerRadius] [outerRadius]" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) ) 
   {
   case 2:
     JunctionAndBranch<2>( argc, argv );
     break;
   case 3:
     JunctionAndBranch<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

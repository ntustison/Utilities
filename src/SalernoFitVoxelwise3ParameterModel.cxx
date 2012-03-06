#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkAmoebaOptimizer.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include "vnl/vnl_vector.h"

class amoebaCostFunction : public itk::SingleValuedCostFunction
{

//
// We are solving the three parameter model (A, B, T1^*) at each
// (x,y) voxel:
//   S(x,y,t_n) = A(x,y) - B(x,y) \times \exp( -t_n / T1^*(x,y) )
//

public:
  typedef amoebaCostFunction                    Self;
  typedef itk::SingleValuedCostFunction         Superclass;
  typedef itk::SmartPointer<Self>               Pointer;
  typedef itk::SmartPointer<const Self>         ConstPointer;

  itkNewMacro( Self );
  itkTypeMacro( amoebaCostFunction, SingleValuedCostFunction );

  typedef Superclass::ParametersType              ParametersType;
  typedef Superclass::DerivativeType              DerivativeType;
  typedef Superclass::MeasureType                 MeasureType;

  typedef vnl_vector<double>                      VectorType;

  enum { SpaceDimension = 3 };

  amoebaCostFunction()
    {
    }

  void SetConstants( VectorType inversionTimes, VectorType intensities )
    {
    if( intensities.size() != inversionTimes.size() )
      {
      itkExceptionMacro( "constant vectors aren't the same length." );
      }

    this->m_InversionTimes = inversionTimes;
    this->m_Intensities = intensities;
    }

  double GetValue( const ParametersType & parameters ) const
    {
    double A = parameters[0];
    double B = parameters[1];
    double T1 = parameters[2];

    double val = 0.0;
    for( unsigned int n = 0; n < this->m_Intensities.size(); n++ )
      {
      val += vnl_math_abs( this->m_Intensities[n] - ( A - B * vcl_exp( -this->m_InversionTimes[n] / T1 ) ) );
      }

    return val;
    }

  void GetDerivative( const ParametersType & parameters,
                            DerivativeType & derivative ) const
    {
    itkExceptionMacro( "Not implemented." );
    }

  unsigned int GetNumberOfParameters(void) const
    {
    return SpaceDimension;
    }

private:

   VectorType               m_InversionTimes;
   VectorType               m_Intensities;
};

int main( unsigned int argc, char *argv[] )
{
  if ( argc < 6 )
    {
    std::cout
      << argv[0] << " outputImagePrefix inputImage1 inversionTime1 "
      << "inputImage2 inversionTime2 ... inputImageN inversionTimeN "
      << std::endl;
    exit( 1 );
    }

  typedef float PixelType;
  typedef itk::Image<PixelType, 2> ImageType;

  typedef itk::AmoebaOptimizer OptimizerType;
  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetMaximumNumberOfIterations( 1000 );
  optimizer->SetParametersConvergenceTolerance( 0.01 );
  optimizer->SetFunctionConvergenceTolerance( 0.001 );

  std::vector<ImageType::Pointer> inputImages;
  std::vector<float> itimes;

  for( unsigned int n = 2; n < argc; n+=2 )
    {
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( argv[n] );
    reader->Update();

    inputImages.push_back( reader->GetOutput() );
    itimes.push_back( atof( argv[n+1] ) );
    }

//  if( itimes.size() != inputImages.size() )
//    {
//    itkExceptionMacro( "Number of inversion times does not match number of images." );
//    }

  ImageType::Pointer A = ImageType::New();
  A->CopyInformation( inputImages[0] );
  A->SetRegions( inputImages[0]->GetLargestPossibleRegion() );
  A->Allocate();
  A->FillBuffer( 0.0 );

  ImageType::Pointer B = ImageType::New();
  B->CopyInformation( inputImages[0] );
  B->SetRegions( inputImages[0]->GetLargestPossibleRegion() );
  B->Allocate();
  B->FillBuffer( 0.0 );

  ImageType::Pointer T1 = ImageType::New();
  T1->CopyInformation( inputImages[0] );
  T1->SetRegions( inputImages[0]->GetLargestPossibleRegion() );
  T1->Allocate();
  T1->FillBuffer( 0.0 );

  amoebaCostFunction::VectorType inversionTimes;
  inversionTimes.set_size( itimes.size() );

  amoebaCostFunction::VectorType intensities;
  intensities.set_size( itimes.size() );

  float maxInversionTime = itimes[0];
  unsigned int max_n = 0;

  for( unsigned int n = 1; n < itimes.size(); n++ )
    {
    inversionTimes[n] = itimes[n];
    if( itimes[n] > maxInversionTime )
      {
      max_n = n;
      maxInversionTime = itimes[n];
      }
    }

  typedef amoebaCostFunction::VectorType VectorType;
  amoebaCostFunction::Pointer costFunction = amoebaCostFunction::New();

  optimizer->SetCostFunction( costFunction.GetPointer() );

  itk::ImageRegionConstIteratorWithIndex<ImageType> It( inputImages[0], inputImages[0]->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    ImageType::IndexType index = It.GetIndex();

    for( unsigned int n = 0; n < inputImages.size(); n++ )
      {
      intensities[n] = inputImages[n]->GetPixel( index );
      }
    costFunction->SetConstants( inversionTimes, intensities );

    OptimizerType::ParametersType initialPosition( 3 );
    initialPosition[0] = intensities[max_n];
    initialPosition[1] = 2 * initialPosition[0];
    initialPosition[2] = inversionTimes[max_n] / vnl_math::ln2;

    optimizer->SetInitialPosition( initialPosition );

    try
      {
      optimizer->StartOptimization();

//      std::cout << optimizer->GetStopConditionDescription() << std::endl;
      }
    catch( itk::ExceptionObject & e )
      {
      std::cerr << "Exception thrown ! " << std::endl;
      std::cerr << "An error ocurred during Optimization" << std::endl;
      std::cerr << "Location    = " << e.GetLocation()    << std::endl;
      std::cerr << "Description = " << e.GetDescription() << std::endl;
      std::cerr <<"[TEST 1 FAILURE]\n";
      return EXIT_FAILURE;
      }

    OptimizerType::ParametersType currentPosition = optimizer->GetCurrentPosition();

    A->SetPixel( index, currentPosition[0] );
    B->SetPixel( index, currentPosition[1] );
    T1->SetPixel( index, currentPosition[2] );
    }

  std::string filenameA = std::string( argv[1] ) + std::string( "A.nii.gz" );
  std::string filenameB = std::string( argv[1] ) + std::string( "B.nii.gz" );
  std::string filenameT1 = std::string( argv[1] ) + std::string( "T1.nii.gz" );

  {
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( filenameA.c_str() );
  writer->SetInput( A );
  writer->Update();
  }

  {
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( filenameB.c_str() );
  writer->SetInput( B );
  writer->Update();
  }

  {
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( filenameT1.c_str() );
  writer->SetInput( T1 );
  writer->Update();
  }

  return 0;
}


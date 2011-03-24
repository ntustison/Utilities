#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLabeledPointSetFileReader.h"

#include <string>
#include <vector>

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
							std::istringstream iss( element );
							iss >> value;
							values.push_back( value );
							}
					}
			return values;
			}

template <unsigned int ImageDimension>
int BSpline( int argc, char *argv[] )
{
  typedef float RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;

  typedef itk::ImageFileReader<RealImageType> ImageReaderType;
  typename ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName( argv[5] );
  reader->Update();

  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;

  typedef itk::PointSet<long, ImageDimension> PointSetType;

  typedef itk::LabeledPointSetFileReader<PointSetType> ReaderType;
  typename ReaderType::Pointer fixedPoints = ReaderType::New();
  fixedPoints->SetFileName( argv[2] );
  fixedPoints->Update();

  typename ReaderType::Pointer movingPoints = ReaderType::New();
  movingPoints->SetFileName( argv[3] );
  movingPoints->Update();

  if( fixedPoints->GetOutput()->GetNumberOfPoints() !=
    movingPoints->GetOutput()->GetNumberOfPoints() )
    {
    std::cerr << "The number of fixed points and moving points must be the same." << std::endl;
    return EXIT_FAILURE;
    }
  if( fixedPoints->GetNumberOfLabels() != movingPoints->GetNumberOfLabels() )
    {
    std::cerr << "The number of fixed and moving labels must be the same." << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::PointSet<VectorType, ImageDimension> DeformationFieldPointSetType;
  typename DeformationFieldPointSetType::Pointer fieldPoints =
    DeformationFieldPointSetType::New();
  fieldPoints->Initialize();
  unsigned long count = 0;

  typename PointSetType::PointsContainerConstIterator fIt =
    fixedPoints->GetOutput()->GetPoints()->Begin();
  typename PointSetType::PointsContainerConstIterator mIt =
    movingPoints->GetOutput()->GetPoints()->Begin();

  while( fIt != fixedPoints->GetOutput()->GetPoints()->End() )
    {
    typename PointSetType::PointType fpoint = fIt.Value();
    typename PointSetType::PointType mpoint = mIt.Value();

    typename DeformationFieldType::PointType point;
    VectorType vector;
    typename DeformationFieldPointSetType::PointType fieldPoint;
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      fieldPoint[i] = mpoint[i];
      vector[i] = fpoint[i] - mpoint[i];
      }
    fieldPoints->SetPoint( count, fieldPoint );
    fieldPoints->SetPointData( count, vector );
    count++;

    ++fIt;
    ++mIt;
    }

  typedef itk::BSplineScatteredDataPointSetToImageFilter
    <DeformationFieldPointSetType, DeformationFieldType> BSplineFilterType;
  typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();

  unsigned int numberOfLevels = atoi( argv[7] );
  unsigned int splineOrder = atoi( argv[6] );

  std::vector<unsigned int> cp
    = ConvertVector<unsigned int>( std::string( argv[8] ) );
  typename BSplineFilterType::ArrayType ncps;

  if ( cp.size() == 1 )
    {
    ncps.Fill( cp[0] );
    }
  else if ( cp.size() == ImageDimension )
    {
    for ( unsigned int d = 0; d < ImageDimension; d++ )
      {
      ncps[d] = cp[d];
      }
    }
  else
    {
    std::cerr << "Invalid ncps format." << std::endl;
    }


  //bspliner->DebugOn();
  bspliner->SetOrigin( reader->GetOutput()->GetOrigin() );
  bspliner->SetSpacing( reader->GetOutput()->GetSpacing() );
  bspliner->SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  bspliner->SetGenerateOutputImage( true );
  bspliner->SetNumberOfLevels( numberOfLevels );
  bspliner->SetSplineOrder( splineOrder );
  bspliner->SetNumberOfControlPoints( ncps );
  bspliner->SetInput( fieldPoints );

  typedef itk::ImageFileWriter<DeformationFieldType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[4] );
  writer->SetInput( bspliner->GetOutput() );
  writer->Update();
  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 8 )
    {
    std::cout << argv[0] << " dimension fixedPointSet movingPointSet outputField "
      << "referenceImage order numberOfevels ncps[0]xncps[1]x..." << std::endl;
    exit( 0 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     BSpline<2>( argc, argv );
     break;
   case 3:
     BSpline<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


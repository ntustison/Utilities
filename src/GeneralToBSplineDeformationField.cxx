#include "itkGeneralToBSplineDeformationFieldFilter.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"

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
int GeneralToBSplineDeformationField( int argc, char *argv[] )
{

  typedef float RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;

  typedef itk::ImageFileReader<DeformationFieldType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typedef itk::GeneralToBSplineDeformationFieldFilter<DeformationFieldType,
    DeformationFieldType> FilterType;
  typename FilterType::ArrayType ncps;

  unsigned int numberOfLevels = atoi( argv[5] );
  unsigned int splineOrder = atoi( argv[4] );

  std::vector<unsigned int> cp
    = ConvertVector<unsigned int>( std::string( argv[6] ) );

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

  typedef itk::PointSet<VectorType, ImageDimension> DeformationFieldPointSetType;
  typename DeformationFieldPointSetType::Pointer fieldPoints =
    DeformationFieldPointSetType::New();
  unsigned long count = 0;

  itk::ImageRegionIteratorWithIndex<DeformationFieldType>
    It( reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion() );
  for ( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    VectorType vector = It.Get();
    if( vector.GetSquaredNorm() == 0 )
      {
      continue;
      }
    typename DeformationFieldType::PointType point;
    reader->GetOutput()->TransformIndexToPhysicalPoint( It.GetIndex(), point );

    if( argc > 7 && atoi( argv[7] ) )
      {
      point += vector;
      itk::ContinuousIndex<double, ImageDimension> cidx;
      reader->GetOutput()->TransformPhysicalPointToContinuousIndex( point, cidx );

      if ( !reader->GetOutput()->GetLargestPossibleRegion().IsInside( cidx ) )
        {
        continue;
        }
      }

    typename DeformationFieldPointSetType::PointType fieldPoint;
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      fieldPoint[i] = point[i];
      }
    fieldPoints->SetPoint( count, fieldPoint );
    if( argc > 7 && atoi( argv[7] ) )
      {
      fieldPoints->SetPointData( count, -vector );
      }
    else
      {
      fieldPoints->SetPointData( count, vector );
      }
    count++;
    }

  typedef itk::BSplineScatteredDataPointSetToImageFilter
    <DeformationFieldPointSetType, DeformationFieldType> BSplineFilterType;
  typename BSplineFilterType::Pointer bspliner = BSplineFilterType::New();

  //bspliner->DebugOn();
  bspliner->SetOrigin( reader->GetOutput()->GetOrigin() );
  bspliner->SetSpacing( reader->GetOutput()->GetSpacing() );
  bspliner->SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  bspliner->SetDirection( reader->GetOutput()->GetDirection() );
  bspliner->SetGenerateOutputImage( true );
  bspliner->SetNumberOfLevels( numberOfLevels );
  bspliner->SetSplineOrder( splineOrder );
  bspliner->SetNumberOfControlPoints( ncps );
  bspliner->SetInput( fieldPoints );
  bspliner->Update();

  typedef itk::ImageFileWriter<DeformationFieldType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( bspliner->GetOutput() );
  writer->Update();

  typename WriterType::Pointer writer2 = WriterType::New();
  writer2->SetFileName( "controlPoints.nii.gz" );
  writer2->SetInput( bspliner->GetPhiLattice() );
  writer2->Update();

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 7 )
    {
    std::cout << argv[0] << " imageDimension inputField outputField order "
              << "numberOfevels ncps[0]xncps[1]x... [doInverse]" << std::endl;
    exit( 0 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     GeneralToBSplineDeformationField<2>( argc, argv );
     break;
   case 3:
     GeneralToBSplineDeformationField<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


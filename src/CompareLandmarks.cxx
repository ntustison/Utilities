#include "itkLabeledPointSetFileReader.h"
#include "itkPointSet.h"

template <unsigned int ImageDimension>
int CompareLandmarks( int argc, char *argv[] )
{
  typedef itk::PointSet<unsigned long, ImageDimension> PointSetType;
  typedef typename PointSetType::PointType PointType;
  typedef itk::LabeledPointSetFileReader<PointSetType> ReaderType;

  typename ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[2] );
  reader1->Update();

  typename ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[3] );
  reader2->Update();

  typename PointSetType::PointsContainerIterator It1
    = reader1->GetOutput()->GetPoints()->Begin();
  typename PointSetType::PointDataContainerIterator ItL1
    = reader1->GetOutput()->GetPointData()->Begin();


  float averageDistance = 0;

  while( It1 != reader1->GetOutput()->GetPoints()->End() )
    {
    unsigned long label = ItL1.Value(); 
    if( label == 1 )
      {
      ++It1;
      ++ItL1;
      continue; 
      } 

    typename PointSetType::PointsContainerIterator It2
      = reader2->GetOutput()->GetPoints()->Begin();
    typename PointSetType::PointDataContainerIterator ItL2
      = reader2->GetOutput()->GetPointData()->Begin();

    while( It2 != reader2->GetOutput()->GetPoints()->End() )
      {
      if( ItL2.Value() == label )
        { 
        float distance = It2.Value().EuclideanDistanceTo( It1.Value() );      
        std::cout << label << ", " << distance << std::endl;

        averageDistance += distance;
        break;
        }
      ++It2;
      ++ItL2;  
      }
    ++It1;
    ++ItL1;
    }
  float N = static_cast<float>( 
    reader1->GetOutput()->GetNumberOfPoints()-1 );
    
  averageDistance /= N;

  float variance = 0.0;

  It1 = reader1->GetOutput()->GetPoints()->Begin();
  ItL1 = reader1->GetOutput()->GetPointData()->Begin();
  while( It1 != reader1->GetOutput()->GetPoints()->End() )
    {
    unsigned long label = ItL1.Value(); 
    if( label == 1 )
      {
      ++It1;
      ++ItL1;
      continue; 
      } 

    typename PointSetType::PointsContainerIterator It2
      = reader2->GetOutput()->GetPoints()->Begin();
    typename PointSetType::PointDataContainerIterator ItL2
      = reader2->GetOutput()->GetPointData()->Begin();

    while( It2 != reader2->GetOutput()->GetPoints()->End() )
      {
      if( ItL2.Value() == label )
        { 
        float distance = It2.Value().EuclideanDistanceTo( It1.Value() );      
        variance += vnl_math_sqr( distance - averageDistance ); 
        break;
        }
      ++It2;
      ++ItL2;  
      }
    ++It1;
    ++ItL1;
    }

  std::cout << "-------------------------------" << std::endl;
  std::cout << averageDistance << " +/- " << vcl_sqrt( variance / ( N - 1 ) ) 
    << std::endl;

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " ImageDimension pointSet1 pointSet2" << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     CompareLandmarks<2>( argc, argv );
     break;
   case 3:
     CompareLandmarks<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


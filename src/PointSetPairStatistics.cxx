#include "itkLabeledPointSetFileReader.h"
#include "itkLabeledPointSetFileWriter.h"

#include <fstream.h>

#include <iomanip.h>


template <unsigned int ImageDimension>
int PointSetPairStatistics( unsigned int argc, char *argv[] )

{
  typedef double RealType;

  typedef itk::PointSet<long, ImageDimension> PointSetType;

  typedef itk::LabeledPointSetFileReader<PointSetType> ReaderType;
  typename ReaderType::Pointer reader1 = ReaderType::New();
  reader1->SetFileName( argv[2] );
  reader1->Update();

  typename ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[3] );
  reader2->Update();

  if( argc > 4 && argv[4] )
    {
    std::cout << "Point sets are paired. " << std::endl;
    if( reader1->GetOutput()->GetNumberOfPoints() !=
      reader2->GetOutput()->GetNumberOfPoints() )
      {
      std::cerr << "Point sets do not have the same number of points." << std::endl;
      return EXIT_FAILURE;
      }

    typename PointSetType::PointsContainerConstIterator ItP1 =
      reader1->GetOutput()->GetPoints()->Begin();
    typename PointSetType::PointsContainerConstIterator ItP2 =
      reader2->GetOutput()->GetPoints()->Begin();

   float minDistance = itk::NumericTraits<float>::max();
   float maxDistance = itk::NumericTraits<float>::NonpositiveMin();
   float N = 0.0;
   float A = 0.0;
   float Q = 0.0;

   while( ItP1 != reader1->GetOutput()->GetPoints()->End() )
     {
     float distance = ( ItP1.Value() ).EuclideanDistanceTo( ItP2.Value() );
     if( distance < minDistance )
       {
       minDistance = distance;
       }
     if( distance > maxDistance )
       {
       maxDistance = distance;
       }
     N += 1.0;

     Q = Q + ( N - 1.0 ) * vnl_math_sqr( distance - A ) / N;
     A = A + ( distance - A ) / N;

     ++ItP1;
     ++ItP2;
     }

    float stdDistance = vcl_sqrt( Q / ( N - 1.0 ) );
    float meanDistance = A;

    std::cout << "      min : " << minDistance << std::endl;
    std::cout << "      max : " << maxDistance << std::endl;
    std::cout << "      mean: " << meanDistance << std::endl;
    std::cout << "      std : " << stdDistance << std::endl;

    }
  else
    {
     //
  /**
   * Calculate maximum of the minimum pairwise distances between the two point-sets.
   */

    std::cout << "Min-Max pairwise distances." << std::endl;
    std::cout << "  Point-set 1" << std::endl;

    for ( unsigned int i = 0; i < reader1->GetNumberOfLabels(); i++ )
      {
      std::vector<float> stdMinDistances;

      typename PointSetType::PointsContainerConstIterator ItP1 =
        reader1->GetOutput()->GetPoints()->Begin();
      typename PointSetType::PointDataContainerIterator ItD1 =
        reader1->GetOutput()->GetPointData()->Begin();
      while( ItP1 != reader1->GetOutput()->GetPoints()->End() )
        {
        if( ItD1.Value() == reader1->GetLabelSet()->operator[](i) )
          {
          float minDistance = itk::NumericTraits<float>::max();

          typename PointSetType::PointsContainerConstIterator ItP2 =
            reader2->GetOutput()->GetPoints()->Begin();
          typename PointSetType::PointDataContainerIterator ItD2 =
            reader2->GetOutput()->GetPointData()->Begin();
          while( ItP2 != reader2->GetOutput()->GetPoints()->End() )
            {
            if( ItD2.Value() == ItD1.Value() )
              {
              float distance = ( ItP1.Value() ).EuclideanDistanceTo( ItP2.Value() );
              if( distance < minDistance )
                {
                minDistance = distance;
                }
              }
            ++ItP2;
            ++ItD2;
            }
          stdMinDistances.push_back( minDistance );
          }
        ++ItP1;
        ++ItD1;
        }

      std::cout << "    Stats for label " << reader1->GetLabelSet()->operator[](i) << std::endl;

      if( stdMinDistances.size() > 0 )
        {
        vnl_vector<float> vnlMinDistances;
        vnlMinDistances.set_size( stdMinDistances.size() );
        std::vector<float>::const_iterator it;
        unsigned long index = 0;
        for( it = stdMinDistances.begin(); it != stdMinDistances.end(); ++it )
          {
          vnlMinDistances.put( index++, *it );
          }
        float mean = vnlMinDistances.mean();
        float variance = 0;
        for( it = stdMinDistances.begin(); it != stdMinDistances.end(); ++it )
          {
          variance += vnl_math_sqr( *it - mean );
          }

        std::cout << "      min : " << vnlMinDistances.min_value() << std::endl;
        std::cout << "      max : " << vnlMinDistances.max_value() << std::endl;
        std::cout << "      mean: " << mean << std::endl;
        std::cout << "      std : " << vcl_sqrt( variance /
          static_cast<float>( vnlMinDistances.size() - 1 ) ) << std::endl;
        }
      }

    std::cout << "  Point-set 2" << std::endl;

    for ( unsigned int i = 0; i < reader2->GetNumberOfLabels(); i++ )
      {
      std::vector<float> stdMinDistances;

      typename PointSetType::PointsContainerConstIterator ItP2 =
        reader2->GetOutput()->GetPoints()->Begin();
      typename PointSetType::PointDataContainerIterator ItD2 =
        reader2->GetOutput()->GetPointData()->Begin();
      while( ItP2 != reader2->GetOutput()->GetPoints()->End() )
        {
        if( ItD2.Value() == reader2->GetLabelSet()->operator[](i) )
          {
          float minDistance = itk::NumericTraits<float>::max();

          typename PointSetType::PointsContainerConstIterator ItP1 =
            reader1->GetOutput()->GetPoints()->Begin();
          typename PointSetType::PointDataContainerIterator ItD1 =
            reader1->GetOutput()->GetPointData()->Begin();
          while( ItP1 != reader1->GetOutput()->GetPoints()->End() )
            {
            if( ItD1.Value() == ItD2.Value() )
              {
              float distance = ( ItP2.Value() ).EuclideanDistanceTo( ItP1.Value() );
              if( distance < minDistance )
                {
                minDistance = distance;
                }
              }
            ++ItP1;
            ++ItD1;
            }
          stdMinDistances.push_back( minDistance );
          }
        ++ItP2;
        ++ItD2;
        }

      std::cout << "    Stats for label " << reader2->GetLabelSet()->operator[](i) << std::endl;

      if( stdMinDistances.size() > 0 )
        {
        vnl_vector<float> vnlMinDistances;
        vnlMinDistances.set_size( stdMinDistances.size() );
        std::vector<float>::const_iterator it;
        unsigned long index = 0;
        for( it = stdMinDistances.begin(); it != stdMinDistances.end(); ++it )
          {
          vnlMinDistances.put( index++, *it );
          }
        float mean = vnlMinDistances.mean();
        float variance = 0;
        for( it = stdMinDistances.begin(); it != stdMinDistances.end(); ++it )
          {
          variance += vnl_math_sqr( *it - mean );
          }

        std::cout << "      min : " << vnlMinDistances.min_value() << std::endl;
        std::cout << "      max : " << vnlMinDistances.max_value() << std::endl;
        std::cout << "      mean: " << mean << std::endl;
        std::cout << "      std : " << vcl_sqrt( variance /
          static_cast<float>( vnlMinDistances.size() - 1 ) ) << std::endl;
        }
      }
    }
  return 0;

}


int main( int argc, char *argv[] )
{
  if ( argc < 4 )

    {

    std::cout << "Usage: " << argv[0] << " imageDimension pointSetFile1 pointSetFile2 [isPaired=0]"
      << std::endl;
    exit( 1 );

    }


  switch( atoi( argv[1] ) )
   {
   case 2:
     PointSetPairStatistics<2>( argc, argv );
     break;
   case 3:
     PointSetPairStatistics<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

#include "itkLabeledPointSetFileReader.h"
#include "itkLabeledPointSetFileWriter.h"

#include <fstream.h>

#include <iomanip.h>


template <unsigned int ImageDimension>
int ConvertLabeledPointSet( unsigned int argc, char *argv[] )

{


  typedef double RealType;

  typedef itk::PointSet<long, ImageDimension> PointSetType;

  typedef itk::LabeledPointSetFileReader<PointSetType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  if ( argc > 4 )
    {
    reader->SetRandomPercentage( atof( argv[4] ) );
    }
  if ( argc > 5 )
    {
    reader->SetExtractBoundaryPoints( atoi( argv[5] ) );
    }
  reader->Update();

  typedef itk::Image<short, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  typename ImageReaderType::Pointer imageReader = ImageReaderType::New();
  if( argc > 7 )
    {
    imageReader->SetFileName( argv[7] );
    imageReader->Update();
    }

  std::cout << "Number of points: "
            << reader->GetOutput()->GetNumberOfPoints() << std::endl;
  std::cout << "Number of labels: " << reader->GetNumberOfLabels() << std::endl;
  std::cout << "Labels: " << std::endl;
  for ( unsigned int i = 0; i < reader->GetNumberOfLabels(); i++ )
    {
    std::cout << "   " << reader->GetLabelSet()->operator[](i) << ": ";

    unsigned long N = 0;

    typename PointSetType::PointsContainerConstIterator ItP =
      reader->GetOutput()->GetPoints()->Begin();
    typename PointSetType::PointDataContainerIterator ItD =
      reader->GetOutput()->GetPointData()->Begin();
    while( ItP != reader->GetOutput()->GetPoints()->End() )
      {
      if( ItD.Value() == reader->GetLabelSet()->operator[](i) )
        {
        N++;
        }
      ++ItP;
      ++ItD;
      }
    std::cout << N << std::endl;
    }


  if( argc > 6 && atoi( argv[6] ) )
    {
    typename PointSetType::Pointer centers = PointSetType::New();
    centers->Initialize();

    for( unsigned int n = 0; n < reader->GetNumberOfLabels(); n++ )
      {
      int currentLabel = reader->GetLabelSet()->operator[](n);
      typename PointSetType::PointType center;
      center.Fill( 0 );
      float N = 0;
      typename PointSetType::PointsContainerConstIterator ItP =
        reader->GetOutput()->GetPoints()->Begin();
      typename PointSetType::PointDataContainerIterator ItD =
        reader->GetOutput()->GetPointData()->Begin();
      while( ItP != reader->GetOutput()->GetPoints()->End() )
        {
        if( ItD.Value() == currentLabel )
          {
          typename PointSetType::PointType point = ItP.Value();
          for( unsigned int d = 0; d < ImageDimension; d++ )
            {
            center[d] += point[d];
            }
          N+=1.0;
          }
        ++ItP;
        ++ItD;
        }
      for( unsigned int d = 0; d < ImageDimension; d++ )
        {
        center[d] /= N;
        }
      centers->SetPoint( n, center );
      centers->SetPointData( n, currentLabel );
      }
    typedef itk::LabeledPointSetFileWriter<PointSetType> WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( argv[3] );
    writer->SetInput( centers );
    if( argc > 7 )
      {
      writer->SetImageSize(
        imageReader->GetOutput()->GetLargestPossibleRegion().GetSize() );
      writer->SetImageOrigin( imageReader->GetOutput()->GetOrigin() );
      writer->SetImageSpacing( imageReader->GetOutput()->GetSpacing() );
      writer->SetImageDirection( imageReader->GetOutput()->GetDirection() );
      }
    writer->Update();
    }
  else
    {
    if( argc > 5 && atof( argv[5] ) < 0 )
      {
      ofstream str( argv[3] );
      str.precision( 7 );
      typename PointSetType::PointsContainerConstIterator ItP =
        reader->GetOutput()->GetPoints()->Begin();
      typename PointSetType::PointDataContainerIterator ItD =
        reader->GetOutput()->GetPointData()->Begin();
      while( ItP != reader->GetOutput()->GetPoints()->End() )
        {
        for( unsigned int d = 0; d < ImageDimension-1; d++ )
          {
          str << std::scientific << std::setw( 16 ) << ItP.Value()[d];
          }
        str << std::scientific << std::setw( 16 ) << ItP.Value()[ImageDimension-1] << std::endl;
        ++ItP;
        ++ItD;
        }
      str.close();
      }
    else
      {
      typedef itk::LabeledPointSetFileWriter<PointSetType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( argv[3] );
      writer->SetInput( reader->GetOutput() );
      if( argc > 7 )
        {
        writer->SetImageSize(
          imageReader->GetOutput()->GetLargestPossibleRegion().GetSize() );
        writer->SetImageOrigin( imageReader->GetOutput()->GetOrigin() );
        writer->SetImageSpacing( imageReader->GetOutput()->GetSpacing() );
        writer->SetImageDirection( imageReader->GetOutput()->GetDirection() );
        }
      writer->Update();
      }
    }

  return 0;

}


int main( int argc, char *argv[] )
{
  if ( argc < 4 )

    {

    std::cout << "Usage: " << argv[0] << " imageDimension inputFile outputFile "
      << "[percentage-[0,1]] [boundaryPointsOnly] [centerOfMassPointsOnly] "
      << "[referenceImage] " << std::endl;
    std::cout << "   If [boundaryPointsOnly] < 0 the output points are in the "
      << "following text file format: " << std::endl;
    std::cout << "      x_0 y_0 [z_0] " << std::endl;
    std::cout << "      x_1 y_1 [z_1] " << std::endl;
    std::cout << "      .   .   .   " << std::endl;
    std::cout << "      x_N y_N [z_N] " << std::endl;

    exit( 1 );

    }


  switch( atoi( argv[1] ) )
   {
   case 2:
     ConvertLabeledPointSet<2>( argc, argv );
     break;
   case 3:
     ConvertLabeledPointSet<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}

#include "itkGeneralToBSplineDeformationFieldFilter.h"

#include "itkImage.h"
#include "itkVectorImageFileReader.h"
#include "itkVectorImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "global.h"

int main( int argc, char *argv[] )        
{
  if ( argc < 6 )
    {
    std::cout << argv[0] << " inputField outputField order numberOfevels ncps" << std::endl;
    exit( 0 );
    } 

  typedef float RealType; 
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> DeformationFieldType;
  
  typedef itk::VectorImageFileReader<RealImageType, DeformationFieldType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->SetUseAvantsNamingConvention( true );
  reader->Update();

  typedef itk::GeneralToBSplineDeformationFieldFilter<DeformationFieldType, 
    DeformationFieldType> FilterType;
  FilterType::ArrayType ncps;
  ncps.Fill( atoi( argv[5] ) );
  FilterType::Pointer bspliner = FilterType::New();
  bspliner->SetInput( reader->GetOutput() );
  bspliner->SetNumberOfLevels( atoi( argv[4] ) );
  bspliner->SetSplineOrder( atoi( argv[3] ) );
  bspliner->SetNumberOfControlPoints( ncps );
  bspliner->SetUseFFDRegularization( false );
//  bspliner->SetIgnorePixelValue( zeroVector ); 
  bspliner->Update();


  typedef itk::VectorImageFileWriter<DeformationFieldType, RealImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( bspliner->GetOutput() );
  writer->SetUseAvantsNamingConvention( true );
  writer->Update();


  return EXIT_SUCCESS;
}     


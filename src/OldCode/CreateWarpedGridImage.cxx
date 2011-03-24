#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorImageFileReader.h"
#include "itkVector.h"

#include "itkWarpImageFilter.h"
#include "itkGridImageSource.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 3 )
    {
    std::cout << "Usage: " << argv[0] << " deformationField outputImage" << std::endl;
    exit( 1 );
    }

  typedef float RealType;
  typedef itk::Image<RealType, ImageDimension> RealImageType;
  typedef itk::Vector<RealType, ImageDimension> VectorType;
  typedef itk::Image<VectorType, ImageDimension> VectorImageType;

  /**
   * Read in vector field
   */
  typedef itk::VectorImageFileReader<RealImageType, VectorImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();
  
  typedef itk::GridImageSource<RealImageType> GridSourceType;
  GridSourceType::Pointer gridder = GridSourceType::New();
  gridder->SetSpacing( reader->GetOutput()->GetSpacing() );
  gridder->SetOrigin( reader->GetOutput()->GetOrigin() );
  gridder->SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );

  GridSourceType::ArrayType gridSpacing;
  GridSourceType::ArrayType sigma;
  GridSourceType::BoolArrayType which;
  which.Fill( false );
  for ( unsigned int i = 0; i < 2; i++ )
    {
    which[i] = true;
    }
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    gridSpacing[i] = reader->GetOutput()->GetLargestPossibleRegion().GetSize()[i]
      *reader->GetOutput()->GetSpacing()[i]/25.0;
    sigma[i] = gridSpacing[i]/10.0;
    }
  gridder->SetGridSpacing( gridSpacing );
  gridder->SetSigma( sigma );
  gridder->SetWhichDimensions( which );
  gridder->Update(); 

  typedef itk::WarpImageFilter<RealImageType, RealImageType, VectorImageType> WarperType;
  WarperType::Pointer warper = WarperType::New();
  warper->SetDeformationField( reader->GetOutput() ); 
  warper->SetInput( gridder->GetOutput() );
  warper->SetOutputOrigin( gridder->GetOutput()->GetOrigin() );
  warper->SetOutputSpacing( gridder->GetOutput()->GetSpacing() );
  warper->Update();
  
  std::string file = std::string( argv[2] );
  typedef itk::ImageFileWriter<RealImageType> ImageWriterType;
  ImageWriterType::Pointer gridWriter = ImageWriterType::New();
  gridWriter->SetFileName( file.c_str() );
  gridWriter->SetInput( warper->GetOutput() );
  gridWriter->Update();

  return 0;
}

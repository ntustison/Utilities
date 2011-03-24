#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkLabelStatisticsImageFilter.h"

#include "global.h"

int main( int argc, char *argv[] )
{
  if ( argc < 2 )
    {
    std::cout << "Usage: " << argv[0] << " image [mask] [label]" << std::endl;
    exit( 1 );
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<unsigned int, ImageDimension> MaskImageType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );  
  reader->Update();
  
  MaskImageType::Pointer mask = MaskImageType::New();
  unsigned int label = 1;
  if ( argc < 3 )
    {
    mask->SetOrigin( reader->GetOutput()->GetOrigin() );
    mask->SetSpacing( reader->GetOutput()->GetSpacing() );
    mask->SetRegions( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
    mask->Allocate();
    mask->FillBuffer( 1 );
    }
  else
    {
    typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
    MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader->SetFileName( argv[2] );  
    maskReader->Update();
    mask = maskReader->GetOutput();

    if ( argc == 4 )
      {
      label = atoi( argv[3] );
      }
    }

  
 
  std::cout << "Image information" << std::endl;
  std::cout << "  Size:     " << reader->GetOutput()->GetLargestPossibleRegion().GetSize() << std::endl;
  std::cout << "  Origin:   " << reader->GetOutput()->GetOrigin() << std::endl;
  std::cout << "  Spacing:  " << reader->GetOutput()->GetSpacing() << std::endl; 
  std::cout << "  Index:    " << reader->GetOutput()->GetLargestPossibleRegion().GetIndex() << std::endl;

  typedef itk::LabelStatisticsImageFilter<ImageType, MaskImageType> StatsType;
  StatsType::Pointer stats = StatsType::New();
  stats->SetInput( reader->GetOutput() );
  stats->SetLabelInput( mask );
  stats->Update();
 
  std::cout << "Statistics for label " << label << std::endl;
  std::cout << "  Maximum:  " << stats->GetMaximum( label ) << std::endl;
  std::cout << "  Minimum:  " << stats->GetMinimum( label ) << std::endl;
  std::cout << "  Mean:     " << stats->GetMean( label ) << std::endl;
  std::cout << "  Sigma:    " << stats->GetSigma( label ) << std::endl;
  std::cout << "  Variance: " << stats->GetVariance( label ) << std::endl;


  return 0;
}

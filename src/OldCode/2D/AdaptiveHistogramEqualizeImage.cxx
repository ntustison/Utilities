#include "itkAdaptiveHistogramEqualizationImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "global.h"



int main (const int argc, const char * argv[])
{
  if ( argc < 6 )
    {
      std::cout << "Usage " << argv[0] 
                << " inputImage outputImage"
                << " alpha(0 = classical histogram equalization -> 1 = unsharp mask)"
                << " beta(0 = unsharp mask -> 1 = pass through)"
                << " radius"
                << std::endl;
      return 1;
    }

  typedef itk::Image<PixelType, ImageDimension> ImageType; 
  
  itk::ImageFileReader <ImageType>::Pointer reader;
  reader= itk::ImageFileReader <ImageType>::New();
  reader->SetFileName ( argv[1] );
  reader->Modified();
  reader->Update();

  ImageType::SizeType radius;
  for ( unsigned int d = 0; d < ImageDimension; d++ )
    {
    radius[d]= atoi( argv[5] );
    }


  itk::AdaptiveHistogramEqualizationImageFilter<ImageType>::Pointer filter;
  filter = itk::AdaptiveHistogramEqualizationImageFilter<ImageType>::New();
  filter->SetInput ( reader->GetOutput() );
  filter->SetAlpha( atof( argv[2] ) );
  filter->SetBeta( atof( argv[3] ) );
  filter->SetRadius( radius );
  filter->Update();

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( argv[2] );
  writer->Update();
  
  return 0;
}

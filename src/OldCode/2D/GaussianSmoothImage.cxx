#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkGaussianOperator.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkNeighborhoodInnerProduct.h"

#include "global.h"

int main( int argc, char ** argv )
{
  if ( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0] << " inputImageFile outputImageFile sigma" << std::endl;
    return -1;
    }

  typedef itk::Image< PixelType, ImageDimension >  ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
 
  typedef itk::ConstNeighborhoodIterator<ImageType> NeighborhoodIteratorType;
  typedef itk::ImageRegionIterator< ImageType>      IteratorType;
  
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();

  ImageType::Pointer output = ImageType::New();
  output->SetRegions( reader->GetOutput()->GetRequestedRegion() );
  output->SetOrigin( reader->GetOutput()->GetOrigin() );
  output->SetSpacing( reader->GetOutput()->GetSpacing() );
  output->Allocate();
  
  itk::NeighborhoodInnerProduct<ImageType> innerProduct;
   
  typedef itk::NeighborhoodAlgorithm
    ::ImageBoundaryFacesCalculator<ImageType> FaceCalculatorType;
  
  FaceCalculatorType faceCalculator;
  FaceCalculatorType::FaceListType faceList;
  FaceCalculatorType::FaceListType::iterator fit;
  
  IteratorType out;
  NeighborhoodIteratorType it;

  itk::GaussianOperator<PixelType, ImageDimension> gaussianOperator;
  gaussianOperator.SetVariance( atof( argv[3] ) * atof( argv[3] ) );

  ImageType::Pointer input = reader->GetOutput();
  for ( unsigned int i = 0; i < ImageDimension; ++i )
    {
    gaussianOperator.SetDirection( i );
    gaussianOperator.CreateDirectional();
    
    faceList = faceCalculator( input, output->GetRequestedRegion(),
                               gaussianOperator.GetRadius() );

    for ( fit = faceList.begin(); fit != faceList.end(); ++fit )
      {
      it = NeighborhoodIteratorType( gaussianOperator.GetRadius(),
                                     input, *fit );

      out = IteratorType( output, *fit );
      
      for ( it.GoToBegin(), out.GoToBegin(); ! it.IsAtEnd(); ++it, ++out )
        {
        out.Set( innerProduct( it, gaussianOperator ) );
        }
      }
    
    // Swap the input and output buffers
    if ( i != ImageDimension - 1 )
      {
      ImageType::Pointer tmp = input;
      input = output;
      output = tmp;
      }
    }
  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput( output );
  writer->Update();
  return 0;
}

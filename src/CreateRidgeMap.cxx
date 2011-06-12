#include "itkAddImageFilter.h"
#include "itkConvolutionImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
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
int CreateRidgeMap3D( unsigned int argc, char *argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<PixelType, ImageDimension-1> SliceType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typename ImageType::Pointer output = ImageType::New();
  output->SetOrigin( reader->GetOutput()->GetOrigin() );
  output->SetSpacing( reader->GetOutput()->GetSpacing() );
  output->SetDirection( reader->GetOutput()->GetDirection() );
  output->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  output->Allocate();
  output->FillBuffer( 0 );

  typename SliceType::SizeType sobelSize;
  typename SliceType::IndexType sobelIndex;
  typename SliceType::RegionType sobelRegion;

  sobelSize.Fill( 3 );
  sobelIndex.Fill( 0 );
  sobelRegion.SetSize( sobelSize );
  sobelRegion.SetIndex( sobelIndex );

  typename SliceType::Pointer sobelX = SliceType::New();
  sobelX->SetRegions( sobelRegion );
  sobelX->Allocate();
  sobelX->FillBuffer( 0 );

  sobelIndex[0] = 0; sobelIndex[1] = 0; sobelX->SetPixel( sobelIndex, -1 );
  sobelIndex[0] = 1; sobelIndex[1] = 0; sobelX->SetPixel( sobelIndex,  0 );
  sobelIndex[0] = 2; sobelIndex[1] = 0; sobelX->SetPixel( sobelIndex,  1 );
  sobelIndex[0] = 0; sobelIndex[1] = 1; sobelX->SetPixel( sobelIndex, -2 );
  sobelIndex[0] = 1; sobelIndex[1] = 1; sobelX->SetPixel( sobelIndex,  0 );
  sobelIndex[0] = 2; sobelIndex[1] = 1; sobelX->SetPixel( sobelIndex,  2 );
  sobelIndex[0] = 0; sobelIndex[1] = 2; sobelX->SetPixel( sobelIndex, -1 );
  sobelIndex[0] = 1; sobelIndex[1] = 2; sobelX->SetPixel( sobelIndex,  0 );
  sobelIndex[0] = 2; sobelIndex[1] = 2; sobelX->SetPixel( sobelIndex,  1 );

  typename SliceType::Pointer sobelY = SliceType::New();
  sobelY->SetRegions( sobelRegion );
  sobelY->Allocate();
  sobelY->FillBuffer( 0 );

  sobelIndex[0] = 0; sobelIndex[1] = 0; sobelY->SetPixel( sobelIndex, -1 );
  sobelIndex[0] = 1; sobelIndex[1] = 0; sobelY->SetPixel( sobelIndex, -2 );
  sobelIndex[0] = 2; sobelIndex[1] = 0; sobelY->SetPixel( sobelIndex, -1 );
  sobelIndex[0] = 0; sobelIndex[1] = 1; sobelY->SetPixel( sobelIndex,  0 );
  sobelIndex[0] = 1; sobelIndex[1] = 1; sobelY->SetPixel( sobelIndex,  0 );
  sobelIndex[0] = 2; sobelIndex[1] = 1; sobelY->SetPixel( sobelIndex,  0 );
  sobelIndex[0] = 0; sobelIndex[1] = 2; sobelY->SetPixel( sobelIndex,  1 );
  sobelIndex[0] = 1; sobelIndex[1] = 2; sobelY->SetPixel( sobelIndex,  2 );
  sobelIndex[0] = 2; sobelIndex[1] = 2; sobelY->SetPixel( sobelIndex,  1 );

  typename SliceType::SizeType ridgeSize;
  typename SliceType::IndexType ridgeIndex;
  typename SliceType::RegionType ridgeRegion;

  ridgeSize[0] = 3;
  ridgeSize[1] = 7;
  ridgeIndex.Fill( 0 );
  ridgeRegion.SetSize( ridgeSize );
  ridgeRegion.SetIndex( ridgeIndex );

  typename SliceType::Pointer ridgeX = SliceType::New();
  ridgeX->SetRegions( ridgeRegion );
  ridgeX->Allocate();
  ridgeX->FillBuffer( 0 );

  ridgeIndex[0] = 0; ridgeIndex[1] = 0; ridgeX->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 1; ridgeIndex[1] = 0; ridgeX->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 2; ridgeIndex[1] = 0; ridgeX->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 3; ridgeIndex[1] = 0; ridgeX->SetPixel( ridgeIndex,  0 );
  ridgeIndex[0] = 4; ridgeIndex[1] = 0; ridgeX->SetPixel( ridgeIndex, -1 );
  ridgeIndex[0] = 5; ridgeIndex[1] = 0; ridgeX->SetPixel( ridgeIndex, -1 );
  ridgeIndex[0] = 6; ridgeIndex[1] = 0; ridgeX->SetPixel( ridgeIndex, -1 );
  ridgeIndex[0] = 0; ridgeIndex[1] = 1; ridgeX->SetPixel( ridgeIndex,  2 );
  ridgeIndex[0] = 1; ridgeIndex[1] = 1; ridgeX->SetPixel( ridgeIndex,  2 );
  ridgeIndex[0] = 2; ridgeIndex[1] = 1; ridgeX->SetPixel( ridgeIndex,  2 );
  ridgeIndex[0] = 3; ridgeIndex[1] = 1; ridgeX->SetPixel( ridgeIndex,  0 );
  ridgeIndex[0] = 4; ridgeIndex[1] = 1; ridgeX->SetPixel( ridgeIndex, -2 );
  ridgeIndex[0] = 5; ridgeIndex[1] = 1; ridgeX->SetPixel( ridgeIndex, -2 );
  ridgeIndex[0] = 6; ridgeIndex[1] = 1; ridgeX->SetPixel( ridgeIndex, -2 );
  ridgeIndex[0] = 0; ridgeIndex[1] = 2; ridgeX->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 1; ridgeIndex[1] = 2; ridgeX->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 2; ridgeIndex[1] = 2; ridgeX->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 3; ridgeIndex[1] = 2; ridgeX->SetPixel( ridgeIndex,  0 );
  ridgeIndex[0] = 4; ridgeIndex[1] = 2; ridgeX->SetPixel( ridgeIndex, -1 );
  ridgeIndex[0] = 5; ridgeIndex[1] = 2; ridgeX->SetPixel( ridgeIndex, -1 );
  ridgeIndex[0] = 6; ridgeIndex[1] = 2; ridgeX->SetPixel( ridgeIndex, -1 );

  ridgeSize[0] = 7;
  ridgeSize[1] = 3;
  ridgeRegion.SetSize( ridgeSize );

  typename SliceType::Pointer ridgeY = SliceType::New();
  ridgeY->SetRegions( ridgeRegion );
  ridgeY->Allocate();
  ridgeY->FillBuffer( 0 );

  ridgeIndex[0] = 0; ridgeIndex[1] = 0; ridgeY->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 1; ridgeIndex[1] = 0; ridgeY->SetPixel( ridgeIndex,  2 );
  ridgeIndex[0] = 2; ridgeIndex[1] = 0; ridgeY->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 0; ridgeIndex[1] = 1; ridgeY->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 1; ridgeIndex[1] = 1; ridgeY->SetPixel( ridgeIndex,  2 );
  ridgeIndex[0] = 2; ridgeIndex[1] = 1; ridgeY->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 0; ridgeIndex[1] = 2; ridgeY->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 1; ridgeIndex[1] = 2; ridgeY->SetPixel( ridgeIndex,  2 );
  ridgeIndex[0] = 2; ridgeIndex[1] = 2; ridgeY->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 0; ridgeIndex[1] = 3; ridgeY->SetPixel( ridgeIndex,  0 );
  ridgeIndex[0] = 1; ridgeIndex[1] = 3; ridgeY->SetPixel( ridgeIndex,  0 );
  ridgeIndex[0] = 2; ridgeIndex[1] = 3; ridgeY->SetPixel( ridgeIndex,  0 );
  ridgeIndex[0] = 0; ridgeIndex[1] = 4; ridgeY->SetPixel( ridgeIndex, -1 );
  ridgeIndex[0] = 1; ridgeIndex[1] = 4; ridgeY->SetPixel( ridgeIndex, -2 );
  ridgeIndex[0] = 2; ridgeIndex[1] = 4; ridgeY->SetPixel( ridgeIndex, -1 );
  ridgeIndex[0] = 0; ridgeIndex[1] = 5; ridgeY->SetPixel( ridgeIndex, -1 );
  ridgeIndex[0] = 1; ridgeIndex[1] = 5; ridgeY->SetPixel( ridgeIndex, -2 );
  ridgeIndex[0] = 2; ridgeIndex[1] = 5; ridgeY->SetPixel( ridgeIndex, -1 );
  ridgeIndex[0] = 0; ridgeIndex[1] = 6; ridgeY->SetPixel( ridgeIndex, -1 );
  ridgeIndex[0] = 1; ridgeIndex[1] = 6; ridgeY->SetPixel( ridgeIndex, -2 );
  ridgeIndex[0] = 2; ridgeIndex[1] = 6; ridgeY->SetPixel( ridgeIndex, -1 );

  unsigned int direction = ( argc > 5 ? atoi( argv[5] ) : 1 );

  typename ImageType::RegionType region;
  typename ImageType::SizeType size
    = output->GetLargestPossibleRegion().GetSize();
  typename ImageType::IndexType index
    = output->GetLargestPossibleRegion().GetIndex();
  size[2] = 0;
  region.SetSize( size );

  for ( int s = output->GetLargestPossibleRegion().GetSize()[2] - 1;
          s >= 0; s-- )
    {
    index[2] = s + output->GetLargestPossibleRegion().GetIndex()[2];
    region.SetIndex( index );

    typedef itk::ExtractImageFilter<ImageType, SliceType> ExtracterType;
    typename ExtracterType::Pointer extracter = ExtracterType::New();
    extracter->SetInput( reader->GetOutput() );
    extracter->SetExtractionRegion( region );
    extracter->SetDirectionCollapseToIdentity();
    extracter->Update();

    typedef itk::DiscreteGaussianImageFilter<SliceType, SliceType> GaussianType;
    typename GaussianType::Pointer gaussian = GaussianType::New();
    gaussian->SetInput( extracter->GetOutput() );
    gaussian->SetUseImageSpacing( true );
    gaussian->SetVariance( ( argc > 4 ? vcl_sqrt( atof( argv[4] ) ) : 1.0 ) );
    gaussian->Update();

    /**
     * Create the signed gradient map
     */

    typedef itk::ConvolutionImageFilter<SliceType, SliceType> ConvolverType;

    typename ConvolverType::Pointer convolverSX = ConvolverType::New();
    convolverSX->SetInput( gaussian->GetOutput() );
    convolverSX->SetImageKernelInput( sobelX );
    convolverSX->Update();

    typename ConvolverType::Pointer convolverSY = ConvolverType::New();
    convolverSY->SetInput( gaussian->GetOutput() );
    convolverSY->SetImageKernelInput( sobelY );
    convolverSY->Update();

    typename SliceType::Pointer gradientMap = SliceType::New();
    gradientMap->SetOrigin( gaussian->GetOutput()->GetOrigin() );
    gradientMap->SetSpacing( gaussian->GetOutput()->GetSpacing() );
    gradientMap->SetDirection( gaussian->GetOutput()->GetDirection() );
    gradientMap->SetRegions( gaussian->GetOutput()->GetLargestPossibleRegion() );
    gradientMap->Allocate();
    gradientMap->FillBuffer( 0 );

    itk::ImageRegionIterator<SliceType> ItG( gradientMap,
      gradientMap->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<SliceType> ItSX( convolverSX->GetOutput(),
      convolverSX->GetOutput()->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<SliceType> ItSY( convolverSY->GetOutput(),
      convolverSY->GetOutput()->GetLargestPossibleRegion() );

    ItG.GoToBegin();
    ItSX.GoToBegin();
    ItSY.GoToBegin();
    while( !ItG.IsAtEnd() )
      {
      float sobelSum = vnl_math_abs( ItSX.Get() )
        + vnl_math_abs( ItSY.Get() );
      float sign = 1.0;
      if( direction == 0 )
        {
        sign = ( ItSX.Get() > 0 ? 1.0 : -1.0 );
        }
      else
        {
        sign = ( ItSY.Get() > 0 ? 1.0 : -1.0 );
        }
      ItG.Set( sign * sobelSum );

      ++ItG;
      ++ItSX;
      ++ItSY;
      }

    /**
     * Create the signed ridge map
     */

    typename ConvolverType::Pointer convolverRX = ConvolverType::New();
    convolverRX->SetInput( gradientMap );
    convolverRX->SetImageKernelInput( ridgeX );
    convolverRX->Update();

    typename ConvolverType::Pointer convolverRY = ConvolverType::New();
    convolverRY->SetInput( gradientMap );
    convolverRY->SetImageKernelInput( ridgeY );
    convolverRY->Update();

    itk::ImageRegionIteratorWithIndex<SliceType> ItRX( convolverRX->GetOutput(),
      convolverRX->GetOutput()->GetLargestPossibleRegion() );
    itk::ImageRegionIteratorWithIndex<SliceType> ItRY( convolverRY->GetOutput(),
      convolverRY->GetOutput()->GetLargestPossibleRegion() );

    ItRX.GoToBegin();
    ItRY.GoToBegin();
    while( !ItRY.IsAtEnd() )
      {
      float ridgeSum = vnl_math_abs( ItRX.Get() )
        + vnl_math_abs( ItRY.Get() );
      float sign = 1.0;
      if( direction == 0 )
        {
        sign = ( ItRX.Get() > 0 ? 1.0 : -1.0 );
        }
      else
        {
        sign = ( ItRY.Get() > 0 ? 1.0 : -1.0 );
        }
      typename ImageType::IndexType outputIndex;
      outputIndex[0] = ItRX.GetIndex()[0];
      outputIndex[1] = ItRX.GetIndex()[1];
      outputIndex[2] = index[2];

      output->SetPixel( outputIndex, sign * ridgeSum );

      ++ItRX;
      ++ItRY;
      }
    }




  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( output );
  writer->Update();

  return 0;
}

template <unsigned int ImageDimension>
int CreateRidgeMap2D( unsigned int argc, char *argv[] )
{
  typedef float PixelType;
  typedef itk::Image<PixelType, ImageDimension> ImageType;
  typedef itk::Image<PixelType, ImageDimension> SliceType;

  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();

  typename ImageType::Pointer output = ImageType::New();
  output->SetOrigin( reader->GetOutput()->GetOrigin() );
  output->SetSpacing( reader->GetOutput()->GetSpacing() );
  output->SetDirection( reader->GetOutput()->GetDirection() );
  output->SetRegions( reader->GetOutput()->GetLargestPossibleRegion() );
  output->Allocate();
  output->FillBuffer( 0 );

  typename SliceType::SizeType sobelSize;
  typename SliceType::IndexType sobelIndex;
  typename SliceType::RegionType sobelRegion;

  sobelSize.Fill( 3 );
  sobelIndex.Fill( 0 );
  sobelRegion.SetSize( sobelSize );
  sobelRegion.SetIndex( sobelIndex );

  typename SliceType::Pointer sobelX = SliceType::New();
  sobelX->SetRegions( sobelRegion );
  sobelX->Allocate();
  sobelX->FillBuffer( 0 );

  sobelIndex[0] = 0; sobelIndex[1] = 0; sobelX->SetPixel( sobelIndex, -1 );
  sobelIndex[0] = 1; sobelIndex[1] = 0; sobelX->SetPixel( sobelIndex,  0 );
  sobelIndex[0] = 2; sobelIndex[1] = 0; sobelX->SetPixel( sobelIndex,  1 );
  sobelIndex[0] = 0; sobelIndex[1] = 1; sobelX->SetPixel( sobelIndex, -2 );
  sobelIndex[0] = 1; sobelIndex[1] = 1; sobelX->SetPixel( sobelIndex,  0 );
  sobelIndex[0] = 2; sobelIndex[1] = 1; sobelX->SetPixel( sobelIndex,  2 );
  sobelIndex[0] = 0; sobelIndex[1] = 2; sobelX->SetPixel( sobelIndex, -1 );
  sobelIndex[0] = 1; sobelIndex[1] = 2; sobelX->SetPixel( sobelIndex,  0 );
  sobelIndex[0] = 2; sobelIndex[1] = 2; sobelX->SetPixel( sobelIndex,  1 );

  typename SliceType::Pointer sobelY = SliceType::New();
  sobelY->SetRegions( sobelRegion );
  sobelY->Allocate();
  sobelY->FillBuffer( 0 );

  sobelIndex[0] = 0; sobelIndex[1] = 0; sobelY->SetPixel( sobelIndex, -1 );
  sobelIndex[0] = 1; sobelIndex[1] = 0; sobelY->SetPixel( sobelIndex, -2 );
  sobelIndex[0] = 2; sobelIndex[1] = 0; sobelY->SetPixel( sobelIndex, -1 );
  sobelIndex[0] = 0; sobelIndex[1] = 1; sobelY->SetPixel( sobelIndex,  0 );
  sobelIndex[0] = 1; sobelIndex[1] = 1; sobelY->SetPixel( sobelIndex,  0 );
  sobelIndex[0] = 2; sobelIndex[1] = 1; sobelY->SetPixel( sobelIndex,  0 );
  sobelIndex[0] = 0; sobelIndex[1] = 2; sobelY->SetPixel( sobelIndex,  1 );
  sobelIndex[0] = 1; sobelIndex[1] = 2; sobelY->SetPixel( sobelIndex,  2 );
  sobelIndex[0] = 2; sobelIndex[1] = 2; sobelY->SetPixel( sobelIndex,  1 );

  typename SliceType::SizeType ridgeSize;
  typename SliceType::IndexType ridgeIndex;
  typename SliceType::RegionType ridgeRegion;

  ridgeSize[0] = 3;
  ridgeSize[1] = 7;
  ridgeIndex.Fill( 0 );
  ridgeRegion.SetSize( ridgeSize );
  ridgeRegion.SetIndex( ridgeIndex );

  typename SliceType::Pointer ridgeX = SliceType::New();
  ridgeX->SetRegions( ridgeRegion );
  ridgeX->Allocate();
  ridgeX->FillBuffer( 0 );

  ridgeIndex[0] = 0; ridgeIndex[1] = 0; ridgeX->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 1; ridgeIndex[1] = 0; ridgeX->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 2; ridgeIndex[1] = 0; ridgeX->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 3; ridgeIndex[1] = 0; ridgeX->SetPixel( ridgeIndex,  0 );
  ridgeIndex[0] = 4; ridgeIndex[1] = 0; ridgeX->SetPixel( ridgeIndex, -1 );
  ridgeIndex[0] = 5; ridgeIndex[1] = 0; ridgeX->SetPixel( ridgeIndex, -1 );
  ridgeIndex[0] = 6; ridgeIndex[1] = 0; ridgeX->SetPixel( ridgeIndex, -1 );
  ridgeIndex[0] = 0; ridgeIndex[1] = 1; ridgeX->SetPixel( ridgeIndex,  2 );
  ridgeIndex[0] = 1; ridgeIndex[1] = 1; ridgeX->SetPixel( ridgeIndex,  2 );
  ridgeIndex[0] = 2; ridgeIndex[1] = 1; ridgeX->SetPixel( ridgeIndex,  2 );
  ridgeIndex[0] = 3; ridgeIndex[1] = 1; ridgeX->SetPixel( ridgeIndex,  0 );
  ridgeIndex[0] = 4; ridgeIndex[1] = 1; ridgeX->SetPixel( ridgeIndex, -2 );
  ridgeIndex[0] = 5; ridgeIndex[1] = 1; ridgeX->SetPixel( ridgeIndex, -2 );
  ridgeIndex[0] = 6; ridgeIndex[1] = 1; ridgeX->SetPixel( ridgeIndex, -2 );
  ridgeIndex[0] = 0; ridgeIndex[1] = 2; ridgeX->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 1; ridgeIndex[1] = 2; ridgeX->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 2; ridgeIndex[1] = 2; ridgeX->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 3; ridgeIndex[1] = 2; ridgeX->SetPixel( ridgeIndex,  0 );
  ridgeIndex[0] = 4; ridgeIndex[1] = 2; ridgeX->SetPixel( ridgeIndex, -1 );
  ridgeIndex[0] = 5; ridgeIndex[1] = 2; ridgeX->SetPixel( ridgeIndex, -1 );
  ridgeIndex[0] = 6; ridgeIndex[1] = 2; ridgeX->SetPixel( ridgeIndex, -1 );

  ridgeSize[0] = 7;
  ridgeSize[1] = 3;
  ridgeIndex.Fill( 0 );
  ridgeRegion.SetSize( ridgeSize );
  ridgeRegion.SetIndex( ridgeIndex );

  typename SliceType::Pointer ridgeY = SliceType::New();
  ridgeY->SetRegions( ridgeRegion );
  ridgeY->Allocate();
  ridgeY->FillBuffer( 0 );

  ridgeIndex[0] = 0; ridgeIndex[1] = 0; ridgeY->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 1; ridgeIndex[1] = 0; ridgeY->SetPixel( ridgeIndex,  2 );
  ridgeIndex[0] = 2; ridgeIndex[1] = 0; ridgeY->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 0; ridgeIndex[1] = 1; ridgeY->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 1; ridgeIndex[1] = 1; ridgeY->SetPixel( ridgeIndex,  2 );
  ridgeIndex[0] = 2; ridgeIndex[1] = 1; ridgeY->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 0; ridgeIndex[1] = 2; ridgeY->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 1; ridgeIndex[1] = 2; ridgeY->SetPixel( ridgeIndex,  2 );
  ridgeIndex[0] = 2; ridgeIndex[1] = 2; ridgeY->SetPixel( ridgeIndex,  1 );
  ridgeIndex[0] = 0; ridgeIndex[1] = 3; ridgeY->SetPixel( ridgeIndex,  0 );
  ridgeIndex[0] = 1; ridgeIndex[1] = 3; ridgeY->SetPixel( ridgeIndex,  0 );
  ridgeIndex[0] = 2; ridgeIndex[1] = 3; ridgeY->SetPixel( ridgeIndex,  0 );
  ridgeIndex[0] = 0; ridgeIndex[1] = 4; ridgeY->SetPixel( ridgeIndex, -1 );
  ridgeIndex[0] = 1; ridgeIndex[1] = 4; ridgeY->SetPixel( ridgeIndex, -2 );
  ridgeIndex[0] = 2; ridgeIndex[1] = 4; ridgeY->SetPixel( ridgeIndex, -1 );
  ridgeIndex[0] = 0; ridgeIndex[1] = 5; ridgeY->SetPixel( ridgeIndex, -1 );
  ridgeIndex[0] = 1; ridgeIndex[1] = 5; ridgeY->SetPixel( ridgeIndex, -2 );
  ridgeIndex[0] = 2; ridgeIndex[1] = 5; ridgeY->SetPixel( ridgeIndex, -1 );
  ridgeIndex[0] = 0; ridgeIndex[1] = 6; ridgeY->SetPixel( ridgeIndex, -1 );
  ridgeIndex[0] = 1; ridgeIndex[1] = 6; ridgeY->SetPixel( ridgeIndex, -2 );
  ridgeIndex[0] = 2; ridgeIndex[1] = 6; ridgeY->SetPixel( ridgeIndex, -1 );

  unsigned short direction = ( argc > 5 ? atoi( argv[5] ) : 1 );

  typedef itk::DiscreteGaussianImageFilter<SliceType, SliceType> GaussianType;
  typename GaussianType::Pointer gaussian = GaussianType::New();
  gaussian->SetInput( reader->GetOutput() );
  gaussian->SetUseImageSpacing( true );
  gaussian->SetVariance( ( argc > 4 ? vnl_math_sqr( atof( argv[4] ) ) : 1.0 ) );
  gaussian->Update();

  /**
    * Create the signed gradient map
    */

  typedef itk::ConvolutionImageFilter<SliceType, SliceType> ConvolverType;

  typename ConvolverType::Pointer convolverSX = ConvolverType::New();
  convolverSX->SetInput( gaussian->GetOutput() );
  convolverSX->SetImageKernelInput( sobelX );
  convolverSX->Update();

  typename ConvolverType::Pointer convolverSY = ConvolverType::New();
  convolverSY->SetInput( gaussian->GetOutput() );
  convolverSY->SetImageKernelInput( sobelY );
  convolverSY->Update();

  typename SliceType::Pointer gradientMap = SliceType::New();
  gradientMap->SetOrigin( gaussian->GetOutput()->GetOrigin() );
  gradientMap->SetSpacing( gaussian->GetOutput()->GetSpacing() );
  gradientMap->SetDirection( gaussian->GetOutput()->GetDirection() );
  gradientMap->SetRegions( gaussian->GetOutput()->GetLargestPossibleRegion() );
  gradientMap->Allocate();
  gradientMap->FillBuffer( 0 );

  itk::ImageRegionIterator<SliceType> ItG( gradientMap,
    gradientMap->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<SliceType> ItSX( convolverSX->GetOutput(),
    convolverSX->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<SliceType> ItSY( convolverSY->GetOutput(),
    convolverSY->GetOutput()->GetLargestPossibleRegion() );

  ItG.GoToBegin();
  ItSX.GoToBegin();
  ItSY.GoToBegin();
  while( !ItG.IsAtEnd() )
    {
    float sobelSum = vnl_math_abs( ItSX.Get() )
      + vnl_math_abs( ItSY.Get() );
    float sign = 1.0;
    if( direction == 0 )
      {
      sign = ( ItSX.Get() > 0 ? 1.0 : -1.0 );
      }
    else
      {
      sign = ( ItSY.Get() > 0 ? 1.0 : -1.0 );
      }
    ItG.Set( sign * sobelSum );

    ++ItG;
    ++ItSX;
    ++ItSY;
    }

  /**
    * Create the signed ridge map
    */

  typename ConvolverType::Pointer convolverRX = ConvolverType::New();
  convolverRX->SetInput( gradientMap );
  convolverRX->SetImageKernelInput( ridgeX );
  convolverRX->Update();

  typename ConvolverType::Pointer convolverRY = ConvolverType::New();
  convolverRY->SetInput( gradientMap );
  convolverRY->SetImageKernelInput( ridgeY );
  convolverRY->Update();

  itk::ImageRegionIteratorWithIndex<SliceType> ItRX( convolverRX->GetOutput(),
    convolverRX->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<SliceType> ItRY( convolverRY->GetOutput(),
    convolverRY->GetOutput()->GetLargestPossibleRegion() );

  ItRX.GoToBegin();
  ItRY.GoToBegin();
  while( !ItRY.IsAtEnd() )
    {
    float ridgeSum = vnl_math_abs( ItRX.Get() )
      + vnl_math_abs( ItRY.Get() );
    float sign = 1.0;
    if( direction == 0 )
      {
      sign = ( ItRX.Get() > 0 ? 1.0 : -1.0 );
      }
    else
      {
      sign = ( ItRY.Get() > 0 ? 1.0 : -1.0 );
      }
    typename ImageType::IndexType outputIndex;
    outputIndex[0] = ItRX.GetIndex()[0];
    outputIndex[1] = ItRX.GetIndex()[1];

    output->SetPixel( outputIndex, sign * ridgeSum );
    output->SetPixel( outputIndex, gradientMap->GetPixel( outputIndex ) );

    ++ItRX;
    ++ItRY;
    }




  typedef itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( output );
  writer->Update();

  return 0;
}


int main( int argc, char *argv[] )
{
  if ( argc < 4 )
    {
    std::cout << argv[0] << " imageDimension inputImage outputImage [sigma] [direction]"
      << std::endl;
    exit( 1 );
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     CreateRidgeMap2D<2>( argc, argv );
     break;
   case 3:
     CreateRidgeMap3D<3>( argc, argv );
     break;
   default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
   }
}


#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include "vnl/vnl_math.h"
#include "vnl/algo/vnl_matrix_inverse.h"

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
      std::istringstream iss2( element );
      iss2 >> value;
      values.push_back( value );
      }
    }
  return values;
}

typedef float RealType;

std::vector<RealType> FitRegressionLine(
  std::vector<RealType> X, std::vector<RealType> Y )
{
  if( X.size() != Y.size() )
    {
    std::cerr << "Vectors are not the same size" << std::endl;
    exit( 1 );
    }

  RealType N = X.size();

  std::vector<RealType>::const_iterator itX;
  std::vector<RealType>::const_iterator itY;

  RealType sumX = 0.0;
  RealType sumY = 0.0;
  RealType sumX2 = 0.0;
  RealType sumXY = 0.0;

  for( itX = X.begin(), itY = Y.begin(); itX != X.end(); ++itX, ++itY )
    {
    sumX  += (*itX);
    sumY  += (*itY);
    sumXY += (*itX) * (*itY);
    sumX2 += (*itX) * (*itX);
    }

  std::vector<RealType> line( 2 );
  line[0] = ( N * sumXY - sumX*sumY ) / ( N *sumX2 - sumX*sumX );
  if( sumX2 == 0 )
    {
    line[1] = sumY / N;
    }
  else
    {
    line[1] = ( sumY - line[0] * sumX ) / N;
    }

  return line;
}

int main( int argc, char *argv[] )
{
  if ( argc < 13 )
    {
    std::cerr << "Usage: " << argv[0] << " longTermImage4D longTermTE1Image4D "
      << "shortTermImage4D maskImage3D description deltaTEValues TEvalues Ycoef "
      << "R2Image3D R2PrimeImage3D LambdaImage3D OEFImage3D"
      << std::endl;
    exit( 1 );
    }

  const unsigned int InputImageDimension = 4;
  const unsigned int OutputImageDimension = 3;

  typedef itk::Image<RealType, InputImageDimension> InputImageType;
  typedef itk::Image<unsigned int, OutputImageDimension> MaskImageType;
  typedef itk::Image<RealType, OutputImageDimension> OutputImageType;

  // Read in the input images

  typedef itk::ImageFileReader<InputImageType> InputImageReaderType;

  InputImageReaderType::Pointer longTermReader = InputImageReaderType::New();
  longTermReader->SetFileName( argv[1] );
  longTermReader->Update();

  unsigned int longTermN =
    longTermReader->GetOutput()->GetLargestPossibleRegion().GetSize()[3];

  InputImageReaderType::Pointer longTermTE1Reader = InputImageReaderType::New();
  longTermTE1Reader->SetFileName( argv[2] );
  longTermTE1Reader->Update();

  unsigned int longTermTE1N =
    longTermTE1Reader->GetOutput()->GetLargestPossibleRegion().GetSize()[3];

  InputImageReaderType::Pointer shortTermReader = InputImageReaderType::New();
  shortTermReader->SetFileName( argv[3] );
  shortTermReader->Update();

  unsigned int shortTermN =
    shortTermReader->GetOutput()->GetLargestPossibleRegion().GetSize()[3];

  // Read in the input mask images

  typedef itk::ImageFileReader<MaskImageType> MaskImageReaderType;
  MaskImageReaderType::Pointer maskReader = MaskImageReaderType::New();
  maskReader->SetFileName( argv[4] );
  maskReader->Update();

  // Handle the strings of values that are passed on the command line
  std::vector<unsigned int> whichTerms = ConvertVector<unsigned int>( std::string( argv[5] ) );
  std::vector<RealType> dDTE = ConvertVector<RealType>( std::string( argv[6] ) );
  std::vector<RealType> TEs = ConvertVector<RealType>( std::string( argv[7] ) );
  RealType Ycoef = atof( argv[8] );

  unsigned int nTEs = 0.5 * whichTerms.size();

  std::vector<RealType> TE_S;
  std::vector<RealType> dTEsq_S;
  std::vector<RealType> TE_LTE1;
  std::vector<RealType> dTE_LTE1;

  unsigned int count = 0;
  vnl_matrix<RealType> A( longTermN, 3, 1.0 );

  std::vector<unsigned int>::const_iterator itw;
  for( itw = whichTerms.begin(); itw != whichTerms.end(); ++itw )
    {
    unsigned int i = itw-whichTerms.begin();
    if( *itw == 1 )
      {
      dTEsq_S.push_back( dDTE[i] * dDTE[i] );
      TE_S.push_back( ( i < nTEs ) ? TEs[0] : TEs[1] );
      }
    if( *itw == 2 || *itw == 3 )
      {
      A(count, 1) = -dDTE[i];
      A(count, 2) = ( i < nTEs ) ? -TEs[0] : -TEs[1];
      count++;
      }
    if( *itw == 3 )
      {
      dTE_LTE1.push_back( dDTE[i] );
      TE_LTE1.push_back( ( i < nTEs ) ? TEs[0] : TEs[1] );
      }
    }

  // Allocate memory for the output images

  OutputImageType::Pointer oefImage = OutputImageType::New();
  oefImage->CopyInformation( maskReader->GetOutput() );
  oefImage->SetRegions( maskReader->GetOutput()->GetLargestPossibleRegion() );
  oefImage->Allocate();
  oefImage->FillBuffer( 0.0 );

  OutputImageType::Pointer r2Image = OutputImageType::New();
  r2Image->CopyInformation( maskReader->GetOutput() );
  r2Image->SetRegions( maskReader->GetOutput()->GetLargestPossibleRegion() );
  r2Image->Allocate();
  r2Image->FillBuffer( 0.0 );

  OutputImageType::Pointer r2pImage = OutputImageType::New();
  r2pImage->CopyInformation( maskReader->GetOutput() );
  r2pImage->SetRegions( maskReader->GetOutput()->GetLargestPossibleRegion() );
  r2pImage->Allocate();
  r2pImage->FillBuffer( 0.0 );

  OutputImageType::Pointer lambdaImage = OutputImageType::New();
  lambdaImage->CopyInformation( maskReader->GetOutput() );
  lambdaImage->SetRegions( maskReader->GetOutput()->GetLargestPossibleRegion() );
  lambdaImage->Allocate();
  lambdaImage->FillBuffer( 0.0 );

  // Iterate over the mask image

  itk::ImageRegionConstIteratorWithIndex<MaskImageType> It(
    maskReader->GetOutput(), maskReader->GetOutput()->GetLargestPossibleRegion() );
  for( It.GoToBegin(); !It.IsAtEnd(); ++It )
    {
    if( It.Get() == 0 )
      {
      continue;
      }
    MaskImageType::IndexType maskIndex = It.GetIndex();
    InputImageType::IndexType index;
    for( unsigned int d = 0; d < OutputImageDimension; d++ )
      {
      index[d] = maskIndex[d];
      }

    vnl_vector<RealType> L( longTermN );
    for( unsigned int i = 0; i < longTermN; i++ )
      {
      index[InputImageDimension-1] = i;
      L(i) = longTermReader->GetOutput()->GetPixel( index );
      }

    vnl_vector<RealType> B = vnl_svd<RealType>( A ).solve( L );

    RealType r2  = B[2];

    std::vector<RealType> LTE1( longTermTE1N );
    for( unsigned int i = 0; i < longTermTE1N; i++ )
      {
      index[InputImageDimension-1] = i;
      RealType tmp = longTermTE1Reader->GetOutput()->GetPixel( index );
      LTE1[i] = tmp + r2 * TE_LTE1[i];
      }

    std::vector<RealType> S( shortTermN );
    for( unsigned int i = 0; i < shortTermN; i++ )
      {
      index[InputImageDimension-1] = i;
      RealType tmp = shortTermReader->GetOutput()->GetPixel( index );
      S[i] = tmp + r2 * TE_S[i];
      }

    std::vector<RealType> lineS = FitRegressionLine( dTEsq_S, S );
    std::vector<RealType> lineLTE1 = FitRegressionLine( dTE_LTE1, LTE1 );

    RealType lambda = lineLTE1[1] - lineS[1];
    RealType r2p = -lineLTE1[0];
    RealType oef = r2p / lambda / Ycoef;

    // I have no idea why the constants are what they are in the lines that
    // follow.

    if( r2p >= 0.0 && r2p <= 15 && oef >= 0.05 && oef <= 1.0 && lambda >= 0.005 )
      {
      oefImage->SetPixel( maskIndex, oef*1000 );
      lambdaImage->SetPixel( maskIndex, lambda*10000 );
      r2Image->SetPixel( maskIndex, r2*100 );
      r2pImage->SetPixel( maskIndex, r2p*100 );
      }
    }

  typedef itk::ImageFileWriter<OutputImageType> WriterType;

  WriterType::Pointer r2Writer = WriterType::New();
  r2Writer->SetInput( r2Image );
  r2Writer->SetFileName( argv[9] );
  r2Writer->Update();

  WriterType::Pointer r2pWriter = WriterType::New();
  r2pWriter->SetInput( r2pImage );
  r2pWriter->SetFileName( argv[10] );
  r2pWriter->Update();

  WriterType::Pointer lambdaWriter = WriterType::New();
  lambdaWriter->SetInput( lambdaImage );
  lambdaWriter->SetFileName( argv[11] );
  lambdaWriter->Update();

  WriterType::Pointer oefWriter = WriterType::New();
  oefWriter->SetInput( oefImage );
  oefWriter->SetFileName( argv[12] );
  oefWriter->Update();
}

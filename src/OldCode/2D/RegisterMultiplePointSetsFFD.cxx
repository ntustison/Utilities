/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: RegisterImagesFFD.cxx,v $
  Language:  C++
  Date:      $Date: $
  Version:   $Revision:  $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itkArray.h"
#include "itkFFDMultiplePointSetsRegistrationFilter.h"
#include "itkImageFileReader.h"
#include "itkLabeledPointSetFileReader.h"
#include "itkLabeledPointSetFileWriter.h"
#include "itkVectorImageFileWriter.h"
#include "itkVariableSizeMatrix.h"

#include <getopt.h>
#include <string>
#include <fstream.h>
#include <iostream>
#include <iomanip>

#include "global.h"

std::vector<unsigned int> parseUIntVector( const std::string & str)
{
   std::vector<unsigned int> vect;

   std::string::size_type crosspos = str.find('x',0);

   if (crosspos == std::string::npos)
   {
      // only one uint
      vect.push_back( static_cast<unsigned int>( atoi(str.c_str()) ));
      return vect;
   }

   // first uint
   vect.push_back( static_cast<unsigned int>(
                      atoi( (str.substr(0,crosspos)).c_str()  ) ));

   while ( true )
   {
      std::string::size_type crossposfrom = crosspos;
      crosspos =  str.find('x',crossposfrom+1);

      if (crosspos == std::string::npos)
      {
         vect.push_back( static_cast<unsigned int>(
                            atoi( (str.substr(crossposfrom+1,str.length()-crossposfrom-1)).c_str()  ) ));
         return vect;
      }

      vect.push_back( static_cast<unsigned int>(
                         atoi( (str.substr(crossposfrom+1,crosspos)).c_str()  ) ));
   }
}

std::vector<float> parseFloatVector( const std::string & str)
{
   std::vector<float> vect;

   std::string::size_type crosspos = str.find('x',0);

   if (crosspos == std::string::npos)
   {
      // only one uint
      vect.push_back( static_cast<unsigned int>( atof(str.c_str()) ));
      return vect;
   }

   // first float
   vect.push_back( atof( (str.substr(0,crosspos)).c_str() ));

   while ( true )
   {
      std::string::size_type crossposfrom = crosspos;
      crosspos =  str.find('x',crossposfrom+1);

      if (crosspos == std::string::npos)
      {
         vect.push_back( atof( (str.substr(crossposfrom+1,str.length()-crossposfrom-1)).c_str() ));
         return vect;
      }

      vect.push_back( atof( (str.substr(crossposfrom+1,crosspos)).c_str()  ));
   }
}

template <unsigned int Dimension>
int RegisterPointSets( unsigned int argc, char *argv[] )
  {

  typedef double                                                           RealType;

  // Set up the registration filter
  std::cout << "Set up the registration filter (Dimension = " << Dimension << ")." << std::endl;

  typedef itk::PointSet<long, Dimension> PointSetType;

  typedef itk::FFDMultiplePointSetsRegistrationFilter<PointSetType> RegistrationFilterType;
  typename RegistrationFilterType::Pointer registrationFilter
    = RegistrationFilterType::New();

  typename RegistrationFilterType::ArrayType array;
  std::vector<unsigned int> meshResolution = parseUIntVector( argv[5] );
  for ( unsigned int i = 0; i < Dimension; i++ )
    {
    array[i] = meshResolution[i];
    }
  registrationFilter->SetInitialMeshResolution( array );

  for ( unsigned int i = 0; i < Dimension; i++ )
    {
    array[i] = 1;
    }
  registrationFilter->SetDirectionality( array );

  typedef itk::Image<float, Dimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  reader->Update();
  typename RegistrationFilterType::SizeType size
    = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
  typename RegistrationFilterType::OriginType origin
    = reader->GetOutput()->GetOrigin();
  typename RegistrationFilterType::SpacingType spacing
    = reader->GetOutput()->GetSpacing();

  registrationFilter->SetSize( size );
  registrationFilter->SetOrigin( origin );
  registrationFilter->SetSpacing( spacing );


  std::vector<unsigned int> iters = parseUIntVector( argv[4] );
  typename RegistrationFilterType::ResizableArrayType numberOfIterations;
  numberOfIterations.SetSize( iters.size() );
  for ( unsigned int i = 0; i < iters.size(); i++ )
    {
    numberOfIterations[i] = iters[i];
    }
  registrationFilter->SetAlpha( atof( argv[6] ) );
  registrationFilter->SetAnnealingRate( atof( argv[3] ) );
  registrationFilter->SetFilePrefix( argv[2] );
  registrationFilter->SetUseInputAsSamples( atoi( argv[8] ) );
  registrationFilter->SetUseAnisotropicCovariances( atoi( argv[7] ) );
  registrationFilter->SetProlificacy( true );
  registrationFilter->SetWhichGradient( 0 );
  registrationFilter->SetNumberOfLevels( iters.size() );
  registrationFilter->SetMaximumNumberOfIterations( numberOfIterations );
  registrationFilter->SetSplineOrder( 3 );
  registrationFilter->SetEmployTerm2( atoi( argv[9] ) );

  registrationFilter->SetRegularizationSigma( atof( argv[10] ) );
  registrationFilter->SetKernelSigma( atof( argv[11] ) );
  registrationFilter->SetNumberOfSamples( atoi( argv[12] ) );
  registrationFilter->SetEvaluationKNeighborhood( atoi( argv[13] ) );
  registrationFilter->SetCovarianceKNeighborhood( atoi( argv[14] ) );
  registrationFilter->SetEmploySteepestDescent( false );
  registrationFilter->SetLineSearchMaximumIterations( 3 );

  // Read point sets
  std::cout << "Reading point sets." << std::endl;
  

  for ( unsigned int i = 15; i < argc; i++ )
    {
    typedef itk::LabeledPointSetFileReader<PointSetType> ReaderType;
    typename ReaderType::Pointer pointsReader = ReaderType::New();
    pointsReader->SetFileName( argv[i] );
    pointsReader->Update();
    
    registrationFilter->SetInput( i-15, pointsReader->GetOutput() );
    }

  std::cout << "Registering " << registrationFilter->GetNumberOfInputs() << " point sets." << std::endl;

  registrationFilter->Update();
  
  // Write the outputs

  std::cout << "Writing the output." << std::endl;


  for ( unsigned int i = 0; i < registrationFilter->GetNumberOfInputs(); i++ )
    {
    std::string file;
  
    typedef itk::VectorImageFileWriter
      <typename RegistrationFilterType::ControlPointLatticeType,
       typename RegistrationFilterType::RealImageType> ControlPointLatticeWriterType;
  
    itk::OStringStream buf;
    buf << "_" << i;
  
    std::string filename = std::string( argv[2] )
      + std::string( "ControlPointLattice" ) + buf.str()
      + std::string( "_3.nii" );
  
    typename ControlPointLatticeWriterType::Pointer cpwriter
      = ControlPointLatticeWriterType::New();
    cpwriter->SetFileName( filename.c_str() );
    cpwriter->SetInput( registrationFilter->GetTotalDeformationFieldControlPoints( i ) );
    cpwriter->Update();
  
    filename = std::string( argv[2] )
      + std::string( "WarpedPoints" ) + buf.str()
      + std::string( ".vtk" );
  
    typedef itk::LabeledPointSetFileWriter<PointSetType> PointSetWriterType;
    typename PointSetWriterType::Pointer pswriter = PointSetWriterType::New();
    pswriter->SetFileName( filename.c_str() );
    pswriter->SetInput( registrationFilter->GetOutput( i ) );
    pswriter->Update();
    }

  return EXIT_SUCCESS;

}

int main( unsigned int argc, char *argv[] )
{
  if ( argc == 1 )
    {
    std::cerr << "Usage: " << argv[0] << " inputDomainImage outputPrefix "
              << "annealingRate iterations bsplineResolution "
              << "alpha useAnisotropicCovariances useInputAsSamples "
              << "useRegularizationTerm pointSetSigma kernelSigma numberOfSamples "
              << "evaluationKNeighborhood covarianceKNeighborhood pointSet[0] ... "
              << "pointSet[n] " << std::endl;
    exit( 0 );
    }

  RegisterPointSets<ImageDimension>( argc, argv );

  return EXIT_SUCCESS;
}


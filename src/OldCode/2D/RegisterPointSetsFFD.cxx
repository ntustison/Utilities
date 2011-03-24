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
#include "itkFFDPointSetRegistrationFilter.h"
#include "itkImageFileReader.h"
#include "itkLabeledPointSetFileReader.h"
#include "itkLabeledPointSetFileWriter.h"
#include "itkTimeProbe.h"
#include "itkVectorImageFileWriter.h"
#include "itkVariableSizeMatrix.h"

#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"

#include <getopt.h>
#include <string>
#include <fstream.h>
#include <iostream>
#include <iomanip>

#include "global.h"

static const char *optString = "f:m:p:x:a:r:c:d:k:g:l:n:i:j:h:s:t:v:R:B:G:D:S:L:I:K:P:Z:O:o:at:ev::h?";

static const struct option longOpts[] = {
   { "fixed-point-set",            required_argument, NULL, 'f' },
   { "moving-point-set",           required_argument, NULL, 'm' },
   { "alpha",                      optional_argument, NULL, 'a' },
   { "expansion-factor",           optional_argument, NULL, 'x' },
   { "annealing-rate",             optional_argument, NULL, 'r' },
   { "use-input-as-samples",       optional_argument, NULL, 'h' },
   { "use-anisotropic-covariances",optional_argument, NULL, 'v' },
   { "prolificacy",                optional_argument, NULL, 'p' },
   { "employ-term2",               optional_argument, NULL, 't' },
   { "regularization-sigma",       optional_argument, NULL, 'c' },
   { "kernel-sigma",               optional_argument, NULL, 'd' },
   { "evaluation k-neighborhood",  optional_argument, NULL, 'k' },
   { "covariance k-neighborhood",  optional_argument, NULL, 'K' },

   { "num-levels",                 required_argument, NULL, 'n' },
   { "num-iterations",             required_argument, NULL, 'i' },
   { "num-samples",                required_argument, NULL, 's' },

   { "bspline-mesh-resolution",    required_argument, NULL, 'R' },
   { "bspline-order",              required_argument, NULL, 'B' },
   { "rueckert-gradient",          optional_argument, NULL, 'G' },
   { "directionality",             optional_argument, NULL, 'D' },
   { "steepest-descent",           optional_argument, NULL, 'S' },
   { "use-line-search",            optional_argument, NULL, 'L' },
   { "max-number-line-search-iterations", optional_argument, NULL, 'l' },

   { "domain-image",               optional_argument, NULL, 'I' },
   { "domain-spacing",             optional_argument, NULL, 'P' },
   { "domain-size",                optional_argument, NULL, 'Z' },
   { "domain-origin",              optional_argument, NULL, 'O' },

   { "outputPrefix",               required_argument, NULL, 'o' },

   { NULL,                         no_argument, NULL, 0 }
};

struct arguments
  {
   std::string  fixedPointSet;                  /* -f option */
   std::string  movingPointSet;                 /* -m option */
   float        alpha;                          /* -a option */
   float        expansionFactor;                /* -x option */
   float        annealingRate;                  /* -r option */
   bool         useInputAsSamples;              /* -h option */
   bool         useAnisotropicCovariances;      /* -v option */
   bool         prolificacy;                    /* -p option */
   bool         employTerm2;                    /* -t option */
   float        regularizationSigma;            /* -c option */
   float        kernelSigma;                     /* -d option */
   unsigned int evaluationkNeighborhood;        /* -k option */
   unsigned int covariancekNeighborhood;        /* -K option */

   unsigned int numLevels;                      /* -n option */
   std::vector<unsigned int> numIterations;     /* -i option */
   unsigned int numFixedSamples;                /* -s option */
   unsigned int numMovingSamples;               /* -j option */

   std::string  bsplineMeshResolution;          /* -R option */
   unsigned int bsplineOrder;                   /* -B option */
   bool         useRueckertGradient;            /* -G option */

   bool employSteepestDescent;                  /* -S option */
   bool employLineSearch;                       /* -L option */
   unsigned int maxNumberOfLineSearchIterations; /* -l option */
   std::string directionality;                  /* -D option */
   std::string domainSpacing;                   /* -P option */
   std::string domainOrigin;                    /* -O option */
   std::string domainSize;                      /* -Z option */
   std::string domainImage;                     /* -I option */
   std::string  outputPrefix;                   /* -o option */

   arguments () :
     fixedPointSet( "" ),
     movingPointSet( "" ),
     alpha( 1.5f ),
     expansionFactor( 0.0 ),
     annealingRate( 0.93 ),
     useInputAsSamples( false ),
     useAnisotropicCovariances( false ),
     prolificacy( true ),
     employTerm2( true ),
     regularizationSigma( 1.0f ),
     kernelSigma( 1.0f ),
     evaluationkNeighborhood( 50u ),
     covariancekNeighborhood( 4u ),

     numLevels( 3u ),
     numFixedSamples( 1000 ),
     numMovingSamples( 1000 ),

     bsplineMeshResolution( "4x4x4" ),
     bsplineOrder( 3 ),
     useRueckertGradient( false ),
     employSteepestDescent( false ),
     employLineSearch( true ),
     maxNumberOfLineSearchIterations( 10 ),
     directionality( "1x1x1" ),

     domainSpacing( "1.0x1.0x1.0" ),
     domainOrigin( "0.0x0.0x0.0" ),
     domainSize( "100x100x100" ),
     domainImage( "" ),

     outputPrefix( "registrationPointSetFFD" )
       {
       numIterations = std::vector<unsigned int>( numLevels, 10u );
       }

   friend std::ostream& operator<< (std::ostream& o, const arguments& args)
     {
     std::ostringstream osstr;
     for ( unsigned int i = 0; i < args.numIterations.size(); ++i )
       {
       osstr << args.numIterations[i] << " ";
       }
     std::string iterstr = "[ " + osstr.str() + "]";

     return o
        << "Arguments structure:" << std::endl
        << "  Fixed point set: " << args.fixedPointSet << std::endl
        << "  Moving point set: " << args.movingPointSet << std::endl
        << "  Alpha: " << args.alpha << std::endl
        << "  Expansion factor: " << args.expansionFactor << std::endl
        << "  Annealing rate: " << args.annealingRate << std::endl
        << "  Use input as samples: " << args.useInputAsSamples << std::endl
        << "  Use anisotropic covariances: " << args.useAnisotropicCovariances << std::endl
        << "  Prolificacy: " << args.prolificacy << std::endl
        << "  Employ term2: " << args.employTerm2 << std::endl
        << "  Regularization sigma: " << args.regularizationSigma << std::endl
        << "  Initialization sigma: " << args.kernelSigma << std::endl
        << "  evaluation k Neighborhood: " << args.evaluationkNeighborhood << std::endl
        << "  covariance k Neighborhood: " << args.covariancekNeighborhood << std::endl

        << "  Output prefix: " << args.outputPrefix << std::endl

        << "  Domain image: " << args.domainImage << std::endl
        << "  Domain spacing: " << args.domainSpacing << std::endl
        << "  Domain size: " << args.domainSize << std::endl
        << "  Domain origin: " << args.domainOrigin << std::endl

        << "  Number of multiresolution levels: " << args.numLevels << std::endl
        << "  Number of  iterations: " <<iterstr << std::endl


        << "  Bspline mesh resolution: " << args.bsplineMeshResolution << std::endl
        << "  Bspline order: " << args.bsplineOrder << std::endl
        << "  Use Reuckert Gradient: " << args.useRueckertGradient << std::endl
        << "  Directionality: " << args.directionality << std::endl
        << "  Employ Steepest Descent: " << args.employSteepestDescent << std::endl
        << "  Employ Line Search: " << args.employLineSearch << std::endl
        << "  Maximum number of line search iterations: " << args.maxNumberOfLineSearchIterations << std::endl;
     }
};

/* Display program usage, and exit.
 */
void display_usage( const std::string progname )
{
   struct arguments defargs = arguments();

   std::ostringstream osstr;
   for ( unsigned int i = 0; i < defargs.numIterations.size(); ++i )
     {
     osstr << defargs.numIterations[i] << " ";
     }
   std::string iterstr = "[ " + osstr.str() + "]";

   std::cout << std::endl;
   std::cout << progname << " - register N point sets using RegisterPointSetsFFD algorithm" << std::endl;
   std::cout << "Usage: "<<progname<<" [OPTION...]" << std::endl
        << "  -f/ Fixed point set: " << defargs.fixedPointSet << std::endl
        << "  -m/ Moving point set: " << defargs.movingPointSet << std::endl
        << "  -a/ Alpha: " << defargs.alpha << std::endl
        << "  -x/ Expansion factor: " << defargs.expansionFactor << std::endl
        << "  -r/ Annealing rate: " << defargs.annealingRate << std::endl
        << "  -h/ Use input as samples: " << defargs.useInputAsSamples << std::endl
        << "  -v/ Use anisotropic covariances: " << defargs.useAnisotropicCovariances << std::endl
        << "  -p/ Prolificacy: " << defargs.prolificacy << std::endl
        << "  -t/ Employ term 2: " << defargs.employTerm2 << std::endl
        << "  -c/ Regularization sigma: " << defargs.regularizationSigma << std::endl
        << "  -d/ Initialization sigma: " << defargs.kernelSigma << std::endl
        << "  -k/ evaluation k Neighborhood: " << defargs.evaluationkNeighborhood << std::endl
        << "  -K/ covariance k Neighborhood: " << defargs.covariancekNeighborhood << std::endl

        << "  -o/ Output prefix: " << defargs.outputPrefix << std::endl

        << "  -I/ Domain image: " << defargs.domainImage << std::endl
        << "  -P/ Domain spacing: " << defargs.domainSpacing << std::endl
        << "  -Z/ Domain size: " << defargs.domainSize << std::endl
        << "  -O/ Domain origin: " << defargs.domainOrigin << std::endl

        << "  -n/ Number of multiresolution levels: " << defargs.numLevels << std::endl
        << "  -i/ Number of  iterations: " <<iterstr << std::endl
        << "  -s/ Number of fixed samples: " <<defargs.numFixedSamples << std::endl
        << "  -j/ Number of moving samples: " <<defargs.numMovingSamples << std::endl

        << "  -R/ Bspline mesh resolution: " << defargs.bsplineMeshResolution << std::endl
        << "  -B/ Bspline order: " << defargs.bsplineOrder << std::endl
        << "  -G/ Use Reuckert Gradient: " << defargs.useRueckertGradient << std::endl
        << "  -D/ Directionality: " << defargs.directionality << std::endl
        << "  -S/ Employ Steepest Descent: " << defargs.employSteepestDescent << std::endl
        << "  -L/ Employ Line Search: " << defargs.employLineSearch << std::endl
        << "  -l/ Max. number of line search iterations: " << defargs.maxNumberOfLineSearchIterations << std::endl;

   std::cout << std::endl;

   exit( EXIT_FAILURE );
};

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



void parseOpts (int argc, char **argv, struct arguments & args)
{

  //itk::Systemfilesystem::path progpath(argv[0]);
   const std::string progname( "FFD Point Set Registration/Normalization" );

   // Default values.
   args = arguments();

   std::vector<unsigned int> defiter = args.numIterations;
   args.numIterations.clear();

   if ( argc == 1 )
     {
     display_usage( progname );
     }

   int opt = 0; /* it's actually going to hold a char */
   int longIndex = 0;

   while ( (opt = getopt_long(argc, argv, optString, longOpts, &longIndex)) != -1 )
     {
      switch( opt )
      {

      case 'f':
        if (! optarg) display_usage(progname);
        args.fixedPointSet = optarg;
        break;

      case 'm':
        if (! optarg) display_usage(progname);
        args.movingPointSet = optarg;
        break;

      case 'a':
         if (! optarg) display_usage(progname);
         args.alpha = atof( optarg );
         break;

      case 'x':
         if (! optarg) display_usage(progname);
         args.expansionFactor = atof( optarg );
         break;

      case 'r':
         if (! optarg) display_usage(progname);
         args.annealingRate = atof( optarg );
         break;

      case 'c':
         if (! optarg) display_usage(progname);
         args.regularizationSigma = atof( optarg );
         break;

      case 'd':
         if (! optarg) display_usage(progname);
         args.kernelSigma = atof( optarg );
         break;

      case 'k':
         if (! optarg) display_usage(progname);
         args.evaluationkNeighborhood = static_cast<unsigned int>( atoi( optarg ) );
         break;

      case 'K':
         if (! optarg) display_usage(progname);
         args.covariancekNeighborhood = static_cast<unsigned int>( atoi( optarg ) );
         break;

      case 'h':
         if (! optarg) display_usage(progname);
         args.useInputAsSamples = static_cast<bool>( atoi( optarg ) );
         break;

      case 'v':
         if (! optarg) display_usage(progname);
         args.useAnisotropicCovariances = static_cast<bool>( atoi( optarg ) );
         break;

      case 'p':
         if (! optarg) display_usage(progname);
         args.prolificacy = static_cast<bool>( atoi( optarg ) );
         break;

      case 't':
         if (! optarg) display_usage(progname);
         args.employTerm2 = static_cast<bool>( atoi( optarg ) );
         break;

      case 'o':
         if (! optarg) display_usage(progname);
         args.outputPrefix = optarg;
         break;

      case 'n':
         if (! optarg) display_usage(progname);
         args.numLevels = static_cast<unsigned int>( atoi(optarg) );
         break;

      case 's':
         if (! optarg) display_usage(progname);
         args.numFixedSamples = static_cast<unsigned int>( atoi(optarg) );
         break;

      case 'j':
         if (! optarg) display_usage(progname);
         args.numMovingSamples = static_cast<unsigned int>( atoi(optarg) );
         break;

      case 'i':
         if (! optarg) display_usage(progname);
         args.numIterations = parseUIntVector(std::string(optarg));
         break;

      case 'R':
         if (! optarg) display_usage( progname );
         args.bsplineMeshResolution = optarg;
         break;

      case 'B':
         if (! optarg) display_usage( progname );
         args.bsplineOrder = atoi( optarg );
         break;

      case 'G':
         if (! optarg) display_usage( progname );
         args.useRueckertGradient = static_cast<bool>( atoi( optarg ) );
         break;
      case 'D':
         if (! optarg) display_usage( progname );
         args.directionality = optarg;
         break;
      case 'S':
         if (! optarg) display_usage( progname );
         args.employSteepestDescent = static_cast<bool>( atoi( optarg ) );
         break;
      case 'L':
         if (! optarg) display_usage( progname );
         args.employLineSearch = static_cast<bool>( atoi( optarg ) );
         break;

      case 'l':
         if (! optarg) display_usage( progname );
         args.maxNumberOfLineSearchIterations = static_cast<unsigned int>( atoi( optarg ) );
         break;

      case 'P':
         if (! optarg) display_usage( progname );
         args.domainSpacing = optarg;
         break;
      case 'O':
         if (! optarg) display_usage( progname );
         args.domainOrigin = optarg;
         break;
      case 'Z':
         if (! optarg) display_usage( progname );
         args.domainSize = optarg;
         break;
      case 'I':
         if (! optarg) display_usage( progname );
         args.domainImage = optarg;
         break;


      case '?':   /* fall-through is intentional */
      default:
         display_usage(progname);
         break;
      }
   }

   if ( args.numLevels == 0)
   {
      std::cout<<"The number of levels should be at least one." << std::endl;
      display_usage(progname);
   }
   if ( args.numIterations.empty() )
   {
      // set a default number of iterations per level
      args.numIterations = std::vector<unsigned int>(args.numLevels, defiter[0]);
   }
   else if ( args.numLevels != args.numIterations.size() )
   {
      std::cout<<"The number of levels and the number of iterations do not match." << std::endl;
      display_usage(progname);
   }
}


template <unsigned int Dimension>
int RegisterPointSets( struct arguments & args )
  {

  std::cout << "Starting RegisterPointSetsFFD with the following arguments: " << std::endl;
  std::cout << args << std::endl << std::endl;

  typedef double                                                           RealType;

  // Set up the registration filter
  std::cout << "Set up the registration filter (Dimension = " << Dimension << ")." << std::endl;

  typedef itk::PointSet<long, Dimension> PointSetType;

  typedef itk::FFDPointSetRegistrationFilter<PointSetType> RegistrationFilterType;
  typename RegistrationFilterType::Pointer registrationFilter
    = RegistrationFilterType::New();

  typename RegistrationFilterType::ArrayType array;
  std::vector<unsigned int> meshResolution = parseUIntVector( args.bsplineMeshResolution );
  for ( unsigned int i = 0; i < Dimension; i++ )
    {
    array[i] = meshResolution[i];
    }
  registrationFilter->SetInitialMeshResolution( array );

  std::vector<unsigned int> directionality = parseUIntVector( args.directionality );
  for ( unsigned int i = 0; i < Dimension; i++ )
    {
    array[i] = directionality[i];
    }
  registrationFilter->SetDirectionality( array );


  if ( !args.domainImage.empty() )
    {
    typedef itk::Image<float, Dimension> ImageType;
    typedef itk::ImageFileReader<ImageType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( args.domainImage );
    reader->Update();
    typename RegistrationFilterType::SizeType size
      = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
    typename RegistrationFilterType::OriginType origin
      = reader->GetOutput()->GetOrigin();
    typename RegistrationFilterType::SpacingType spacing
      = reader->GetOutput()->GetSpacing();

    RealType factor = args.expansionFactor;
    for ( unsigned int d = 0; d < Dimension; d++ )
      {
      origin[d] -= factor * spacing[d];
      size[d] += static_cast<unsigned int>( 2*factor );
      }

    registrationFilter->SetSize( size );
    registrationFilter->SetOrigin( origin );
    registrationFilter->SetSpacing( spacing );
    }
  else if ( args.expansionFactor == 0 )
    {
    std::vector<unsigned int> domainSize = parseUIntVector( args.domainSize );
    typename RegistrationFilterType::SizeType size;
    for ( unsigned int i = 0; i < Dimension; i++ )
      {
      size[i] = domainSize[i];
      }
    registrationFilter->SetSize( size );
    std::vector<float> domainOrigin = parseFloatVector( args.domainOrigin );
    typename RegistrationFilterType::OriginType origin;
    for ( unsigned int i = 0; i < Dimension; i++ )
      {
      origin[i] = domainOrigin[i];
      }
    registrationFilter->SetOrigin( origin );
    std::vector<float> domainSpacing = parseFloatVector( args.domainSpacing );
    typename RegistrationFilterType::SpacingType spacing;
    for ( unsigned int i = 0; i < Dimension; i++ )
      {
      spacing[i] = domainSpacing[i];
      }
    registrationFilter->SetSpacing( spacing );
    }
  else
    {
    registrationFilter->SetExpansionFactor( args.expansionFactor );
    }

  typename RegistrationFilterType::ResizableArrayType numberOfIterations;
  numberOfIterations.SetSize( args.numLevels );
  for ( unsigned int i = 0; i < args.numLevels; i++ )
    {
    numberOfIterations[i] = args.numIterations[i];
    }
  registrationFilter->SetAlpha( args.alpha );
  registrationFilter->SetAnnealingRate( args.annealingRate );
  registrationFilter->SetFilePrefix( args.outputPrefix );
  registrationFilter->SetUseInputAsSamples( args.useInputAsSamples );
  registrationFilter->SetUseAnisotropicCovariances( args.useAnisotropicCovariances );
  registrationFilter->SetProlificacy( args.prolificacy );
  registrationFilter->SetWhichGradient( args.useRueckertGradient );
  registrationFilter->SetNumberOfLevels( args.numLevels );
  registrationFilter->SetMaximumNumberOfIterations( numberOfIterations );
  registrationFilter->SetNumberOfFixedSamples( args.numFixedSamples );
  registrationFilter->SetNumberOfMovingSamples( args.numMovingSamples );
  registrationFilter->SetSplineOrder( args.bsplineOrder );
  registrationFilter->SetRegularizationSigma( args.regularizationSigma );
  registrationFilter->SetKernelSigma( args.kernelSigma );
  registrationFilter->SetEvaluationKNeighborhood( args.evaluationkNeighborhood );
  registrationFilter->SetCovarianceKNeighborhood( args.covariancekNeighborhood );
  registrationFilter->SetEmploySteepestDescent( args.employSteepestDescent );
  registrationFilter->SetEmployTerm2( args.employTerm2 );
  registrationFilter->SetDirectionality( array );
  if ( args.employLineSearch )
    {
    registrationFilter->SetLineSearchMaximumIterations( args.maxNumberOfLineSearchIterations );
    }
  else
    {
    registrationFilter->SetLineSearchMaximumIterations( 0 );
    }

  // Read point sets
  std::cout << "Reading point sets." << std::endl;
  
  typedef itk::LabeledPointSetFileReader<PointSetType> ReaderType;
  
  typename ReaderType::Pointer fixedPointsReader = ReaderType::New();
  fixedPointsReader->SetFileName( args.fixedPointSet.c_str() );
  fixedPointsReader->Update();

  typename ReaderType::Pointer movingPointsReader = ReaderType::New();
  movingPointsReader->SetFileName( args.movingPointSet.c_str() );
  movingPointsReader->Update();
  
  registrationFilter->SetInput( 0, fixedPointsReader->GetOutput() );
  registrationFilter->SetInput( 1, movingPointsReader->GetOutput() );

  std::cout << "Number of fixed points: " << fixedPointsReader->GetOutput()->GetNumberOfPoints() << std::endl;
  std::cout << "    Number of fixed labels: "
            << fixedPointsReader->GetLabelSet()->size() << std::endl;
  std::cout << "    Distinct fixed labels: ";
  for ( unsigned int n = 0;
        n < fixedPointsReader->GetLabelSet()->size(); n++ )
    {
    std::cout << fixedPointsReader->GetLabelSet()->operator[]( n ) << " ";
    }
  std::cout << endl;  
  std::cout << "Number of moving points: " << movingPointsReader->GetOutput()->GetNumberOfPoints() << std::endl;
  std::cout << "    Number of moving labels: "
            << movingPointsReader->GetLabelSet()->size() << std::endl;
  std::cout << "    Distinct moving labels: ";
  for ( unsigned int n = 0;
        n < movingPointsReader->GetLabelSet()->size(); n++ )
    {
    std::cout << movingPointsReader->GetLabelSet()->operator[]( n ) << " ";
    }
  std::cout << endl;  


  std::cout << "Registering " << registrationFilter->GetNumberOfInputs() << " point sets." << std::endl;

  itk::TimeProbe timer;
  timer.Start();
  registrationFilter->Update();
  timer.Stop();
  std::cout << "Registration filter run time = " << timer.GetMeanTime() << std::endl;

  // Write the outputs

  std::cout << "Writing the output." << std::endl;

  std::string file;

  typedef itk::VectorImageFileWriter
    <typename RegistrationFilterType::ControlPointLatticeType,
     typename RegistrationFilterType::RealImageType> ControlPointLatticeWriterType;

  itk::OStringStream buf;
  buf << "_" << args.bsplineOrder;

  std::string filename = std::string( args.outputPrefix )
    + std::string( "ControlPointLattice" ) + buf.str()
    + std::string( ".nii" );

  typename ControlPointLatticeWriterType::Pointer cpwriter
    = ControlPointLatticeWriterType::New();
  cpwriter->SetFileName( filename.c_str() );
  cpwriter->SetInput( registrationFilter->GetTotalDeformationFieldControlPoints() );
  cpwriter->Update();

  filename = std::string( args.outputPrefix )
    + std::string( "WarpedPoints" ) + buf.str()
    + std::string( ".vtk" );

  typedef itk::LabeledPointSetFileWriter<PointSetType> PointSetWriterType;
  typename PointSetWriterType::Pointer pswriter = PointSetWriterType::New();
  pswriter->SetFileName( filename.c_str() );
  pswriter->SetInput( registrationFilter->GetOutput() );
  pswriter->Update();

  return EXIT_SUCCESS;

}

int main( int argc, char *argv[] )
{
  if ( argc == 1 )
    {
    display_usage( "RegisterPointSetsFFD" );
    exit( 0 );
    }

  struct arguments args;
  parseOpts( argc, argv, args );

  RegisterPointSets<ImageDimension>( args );

  return EXIT_SUCCESS;
}


/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: RegisterImagesFFD.cxx,v $
  Language:  C++
  Date:      $Date: 2008/11/20 19:03:58 $
  Version:   $Revision: 1.6 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itkArray.h"
#include "itkCheckerBoardImageFilter.h"
#include "itkVectorFieldGradientImageFunction.h"
#include "itkVectorImageFileWriter.h"
#include "itkVectorImageFileReader.h"
#include "itkFFDRegistrationFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkTimeProbe.h"
#include "itkVectorImageFileReader.h"
#include "itkWarpImageFilter.h"

// Image similarity metrics
#include "itkAvantsMutualInformationRegistrationFunction.h"
#include "itkMeanSquareRegistrationFunction.h"
#include "itkAvantsPDEDeformableRegistrationFunction.h"
#include "itkProbabilisticRegistrationFunction.h"
#include "itkRobustDemonsRegistrationFunction.h"
#include "itkRobustOpticalFlow.h"
#include "itkSectionMutualInformationRegistrationFunction.h"

#include <getopt.h>
#include <string>
#include <iostream>
#include <iomanip>

static const char *optString = "a:f:m:j:k:w:c:d:o:O:r:b:i:s:g:x:p:y:z:l:R:B:G:J:X:N:I:C:L:D:S:P:Q:F:A:at:ev::h?";

static const struct option longOpts[] = {
   { "fixed-image",                required_argument, NULL, 'f' },
   { "moving-image",               required_argument, NULL, 'm' },
   { "fixed-image2",               optional_argument, NULL, 'j' },   /* not used */
   { "moving-image2",              optional_argument, NULL, 'k' },   /* not used */
   { "weight-image",               optional_argument, NULL, 'w' },   /* not used */
   { "mask-image",                 optional_argument, NULL, 'x' },
   { "fixed-landmarks",            optional_argument, NULL, 'y' },
   { "moving-landmarks",           optional_argument, NULL, 'z' },
   { "similarity-term",            required_argument, NULL, 'c' },
   { "similarity-term2",           optional_argument, NULL, 'd' },   /* not used */
   { "output-image",               required_argument, NULL, 'o' },
   { "output-field",               optional_argument, NULL, 'O' },
   { "true-field",                 required_argument, NULL, 'r' },   /* not used */
   { "num-levels",                 required_argument, NULL, 'n' },
   { "num-iterations",             required_argument, NULL, 'i' },
   { "number-to-interpolate",      required_argument, NULL, 'b' },   /* not used */
   { "def-field-sigma",            required_argument, NULL, 's' },   /* not used */
   { "up-field-sigma",             required_argument, NULL, 'g' },   /* not used */
   { "affine-transform",           optional_argument, NULL, 'a' },
   { "use-tensors",                optional_argument, NULL, 't' },   /* not used */
   { "use-histogram-matching",     no_argument,       NULL, 'e' },
   { "prolificacy",                optional_argument, NULL, 'p' },
   { "help",                       no_argument,       NULL, 'h' },

   { "input-control-point-lattice",optional_argument, NULL, 'A' },
   { "bspline-mesh-resolution",    required_argument, NULL, 'R' },
   { "bspline-order",              required_argument, NULL, 'B' },
   { "rueckert-gradient",          optional_argument, NULL, 'G' },
   { "jacobian-regularization",    optional_argument, NULL, 'J' },
   { "max-jacobian",               optional_argument, NULL, 'X' },
   { "min-jacobian",               optional_argument, NULL, 'N' },
   { "interior-penalty-weight",    optional_argument, NULL, 'I' },
   { "control-point-file",         optional_argument, NULL, 'C' },
   { "directionality",             optional_argument, NULL, 'D' },
   { "directionality-field",       optional_argument, NULL, 'F' },
   { "steepest-descent",           optional_argument, NULL, 'S' },
   { "use-line-search",            optional_argument, NULL, 'L' },
   { "max-number-line-search-iterations", optional_argument, NULL, 'l' },
   { "output-performance-file",    optional_argument, NULL, 'P' },


   { NULL,                         no_argument, NULL, 0 }
};

struct arguments
  {

   // Brian's arguments
   std::string  fixedImageFile;                 /* -f option */
   std::string  movingImageFile;                /* -m option */
   std::string  fixedImageFile2;                /* -j option */
   std::string  movingImageFile2;               /* -k option */
   std::string  weightImageFile;                /* -w option */
   std::string  outputImageFile;                /* -o option */
   std::string  outputFieldFile;                /* -O option */
   std::string  trueFieldFile;                  /* -r option */  /* not used */
   std::vector<unsigned int> numIterations;     /* -i option */
   float sigmaDef;                              /* -s option */  /* not used */
   float sigmaUp;                               /* -g option */  /* not used */
   float maxStepLength;                         /* -l option */  /* not used */
   bool affineTransform;                        /* -a option */
   bool useTensors;                             /* -t option */  /* not used */
   bool useHistogramMatching;                   /* -e option */
   unsigned int verbosity;                      /* -v option */  /* not used */
   unsigned int  whichMetric[2];                /* -c/-d option */
   unsigned int  interpolateNimages;            /* -b option */  /* not used */
   std::string fixedLandmarksFile;              /* -y option */
   std::string movingLandmarksFile;             /* -z option */

   // Nick's arguments
   std::string bsplineMeshResolution;           /* -R option */
   unsigned int bsplineOrder;                   /* -B option */
   bool useRueckertGradient;                    /* -G option */

   bool useJacobianRegularization;              /* -J option */
   float maxJacobian;                           /* -X option */
   float minJacobian;                           /* -N option */
   float interiorPenaltyWeight;                 /* -I option */
   std::string inputControlPointLatticeFile;    /* -Q option */
   std::string directionality;                  /* -D option */
   
   bool employSteepestDescent;                  /* -S option */
   bool employLineSearch;                       /* -L option */
   unsigned int maxNumberOfLineSearchIterations; /* -l option */
   std::string outputPerformanceFile;           /* -P option */
   bool prolificacy;                             /* -p option */

   arguments () :
     fixedImageFile( "" ),
     movingImageFile( "" ),
     fixedImageFile2( "" ),
     movingImageFile2( "" ),
     weightImageFile( "" ),
     outputImageFile( "registrationFFD" ),
     outputFieldFile( "" ),
     trueFieldFile( "" ),
     sigmaDef( 3.0f ),
     sigmaUp( 0.0f ),
     maxStepLength( 0.5f ),
     affineTransform( false ),
     useTensors( 0u ),
     useHistogramMatching( false ),
     verbosity( 0u ),
     interpolateNimages( 0 ),
     fixedLandmarksFile( "" ),
     movingLandmarksFile( "" ),

     bsplineMeshResolution( "4x4x4" ),
     bsplineOrder( 3 ),
     useRueckertGradient( false ),
     useJacobianRegularization( false ),
     maxJacobian( 10.0f ),
     minJacobian( 0.1f ),
     interiorPenaltyWeight( 1.0f ),
     inputControlPointLatticeFile( "" ),
     directionality( "1x1x1" ),

     employSteepestDescent( false ),
     employLineSearch( true ),
     maxNumberOfLineSearchIterations( 10 ),
     outputPerformanceFile( "" ),
     prolificacy( false )
       {
       numIterations = std::vector<unsigned int>( 3, 10u );
       whichMetric[0] = 4;
       whichMetric[1] = 99;
       }


   friend std::ostream& operator<< (std::ostream& o, const arguments& args)
     {
     std::ostringstream osstr;
     for ( unsigned int i = 0; i < args.numIterations.size(); ++i )
       {
       osstr << args.numIterations[i] << " ";
       }
     std::string iterstr = "[ " + osstr.str() + "]";

     if ( args.affineTransform )
       {
       return o
          << "Arguments structure:" << std::endl
          << "  Fixed image file: " << args.fixedImageFile << std::endl
          << "  Moving image file: " << args.movingImageFile << std::endl
          << "  2nd fixed image file: " << args.fixedImageFile2 << std::endl
          << "  2nd moving image file: " << args.movingImageFile2 << std::endl
          << "  Fixed landmark file: " << args.fixedLandmarksFile << std::endl
          << "  Moving landmark file: " << args.movingLandmarksFile << std::endl
          << "  Weight image file: " << args.weightImageFile << std::endl
          << "  Output image file prefix: " << args.outputImageFile << std::endl
          << "  Input control point lattice file: " << args.inputControlPointLatticeFile << std::endl

          << std::endl

          << "  Use histogram matching: " << args.useHistogramMatching << std::endl
          << "  Number of  iterations: " <<iterstr << std::endl
          << "  Cost/Similarity term: " << args.whichMetric[0] << std::endl
          << "  Cost/Similarity term2: " << args.whichMetric[1] << std::endl
          << "  Use Reuckert Gradient: " << args.useRueckertGradient << std::endl
          << "  Directionality: " << args.directionality << std::endl
          << "  Employ Steepest Descent: " << args.employSteepestDescent << std::endl
          << "  Output performance file: " << args.outputPerformanceFile << std::endl
          << "  Prolificacy: " << args.prolificacy << std::endl

          << std::endl

          << "Affine transformation: " << std::endl
          << "  B-spline order: 1" << std::endl
          << "  B-spline mesh resolution: 2^(n-D) patches" << std::endl
          << "  No line search " << std::endl;
       }
     else
       {
       return o
          << "Arguments structure:" << std::endl
          << "  Fixed image file: " << args.fixedImageFile << std::endl
          << "  Moving image file: " << args.movingImageFile << std::endl
          << "  2nd fixed image file: " << args.fixedImageFile2 << std::endl
          << "  2nd moving image file: " << args.movingImageFile2 << std::endl
          << "  Fixed landmark file: " << args.fixedLandmarksFile << std::endl
          << "  Moving landmark file: " << args.movingLandmarksFile << std::endl
          << "  Weight image file: " << args.weightImageFile << std::endl
          << "  Output image file prefix: " << args.outputImageFile << std::endl
          << "  Input control point lattice file: " << args.inputControlPointLatticeFile << std::endl

          << std::endl

          << "  Use histogram matching: " << args.useHistogramMatching << std::endl
          << "  Number of  iterations: " <<iterstr << std::endl
          << "  Cost/Similarity term: " << args.whichMetric[0] << std::endl
          << "  Cost/Similarity term2: " << args.whichMetric[1] << std::endl
          << "  Use Reuckert Gradient: " << args.useRueckertGradient << std::endl
          << "  Directionality: " << args.directionality << std::endl
          << "  Employ Steepest Descent: " << args.employSteepestDescent << std::endl
          << "  Output performance file: " << args.outputPerformanceFile << std::endl
          << "  Prolificacy: " << args.prolificacy << std::endl

          << std::endl

          << "General B-spline non-rigid registration" << std::endl
          << "  B-spline order: " << args.bsplineOrder << std::endl
          << "  B-spline mesh resolution: " << args.bsplineMeshResolution << std::endl
          << "  Use Jacobian Regularization: " << args.useJacobianRegularization << std::endl
          << "  Max Jacobian: " << args.maxJacobian << std::endl
          << "  Min Jacobian: " << args.minJacobian << std::endl
          << "  Interior penalty weight: " << args.interiorPenaltyWeight << std::endl
          << "  Employ Line Search: " << args.employLineSearch << std::endl
          << "  Maximum number of line search iterations: " << args.maxNumberOfLineSearchIterations << std::endl;
       }
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
   std::cout << progname << " - register 2 images using MVSyN algorithm" << std::endl;
   std::cout << "Usage: "<<progname<<" [OPTION...]" << std::endl;

   std::cout << " IMPORTANT stuff "  << std::endl;
   std::cout << "  -f/--fixed-image = STRING                     Fixed image filename (MANDATORY)"            << std::endl;
   std::cout << "  -m/--moving-image = STRING                    Moving image filename (MANDATORY)"           << std::endl;
   std::cout << "  -o/--output-image = STRING                    Output image filename, default: "            << defargs.outputImageFile << std::endl;
   std::cout << "  -i/--num-iterations = UINTx...xUINT           Number of iterations - default: "            << iterstr << std::endl;
   std::cout << "  -R/--bspline-mesh-resolution = UINTx...xUINT  B-spline mesh resolution, default: "         << defargs.bsplineMeshResolution << std::endl;
   std::cout << "  -B/--bspline-order = UINT                     B-spline order, default: "                   << defargs.bsplineOrder << std::endl;
   std::cout << "  -A/--input-control-point = STRING             Input control point file, default: "         << defargs.inputControlPointLatticeFile << std::endl;
   std::cout << "  -c/--similarity-term = UINT                   Cost/similarity term, default:  "            << defargs.whichMetric[0]  << std::endl;
   std::cout << "  -d/--similarity-term = UINT                   Cost/similarity term 2, default:  "          << defargs.whichMetric[1] << std::endl;
   std::cout << "  Cost/similarity Options: " << std::endl;
   std::cout << "     0: optical flow" << std::endl;
   std::cout << "     1: n.a." << std::endl;
   std::cout << "     2: Mutual Information (not good)" << std::endl;
   std::cout << "     3: MutualInformation (good - 2nd most recommended) " << std::endl;
   std::cout << "     4: Cross Correlation of radius 5 (recommended) "     << std::endl;
   std::cout << "     5: Cross Correlation of radius 2 (similar to 4) "    << std::endl;
   std::cout << "  -j/--fixed-image2 = STRING           2nd Fixed image filename (OPTIONAL)"                  << std::endl;
   std::cout << "  -k/--moving-image2 = STRING          2nd Moving image filename (OPTIONAL)"                 << std::endl;
//   std::cout << "  -t/--use-tensors = UINT              boolean determining if 2nd images are tensors - OPT " << std::endl;
   std::cout << "  -w/--weight-image = STRING           weight image filename (OPTIONAL), default 0.5 "       << std::endl;

   std::cout << " MEDIUM importance stuff " << std::endl;
   std::cout << "  -x/--mask-image = STRING                      Mask image filename "                       << std::endl;
   std::cout << "  -y/--fixed-landmark = STRING                  Fixed landmark file "                       << std::endl;
   std::cout << "  -z/--moving-landmark = STRING                 Moving landmark file "                      << std::endl;
   std::cout << "  -G/rueckert-gradient = BOOL                   Use Rueckert Gradient, default: false "     << std::endl;
   std::cout << "  -J/jacobian-regularization = BOOL             Use jacobian reg., default: false "         << std::endl;
   std::cout << "  -X/max-jacobian = FLOAT                       Max jacobian, default: 10 "                 << std::endl;
   std::cout << "  -N/min-jacobian = FLOAT                       Max jacobian, default: 0.1 "                << std::endl;
   std::cout << "  -I/interior-penalty-weight = FLOAT            Jacobian reg. weight, default: 1.0 "        << std::endl;
   std::cout << "  -S/employ-steepest-descent = BOOL             Use steepest descent, default: false "      << std::endl;
   std::cout << "  -L/employ-line-search = BOOL                  Use line search, default: true "      << std::endl;
   std::cout << "  -l/Max. number of line search iterations = INT " << std::endl;

   std::cout << " LESS important stuff " << std::endl;
   std::cout << "  -O/--output-field = STRING                    Output field filename, default: deformationField.nii"  << std::endl;
   std::cout << "  -e/--use-histogram-matching = BOOL            Use histogram matching, default: false"                << std::endl;
   std::cout << "  -D/directionality = STRING                    Directionality, default: 1x1x1 "                       << std::endl;
   std::cout << "  -P/performanceFile = STRING                   Output performance filename "                          << std::endl;
   std::cout << "  -p/prolificacy = BOOL                         Print intermediate stage images "                      << std::endl;
   std::cout << "  -h/--help                                     Display this message and exit"                         << std::endl;


//   std::cout << " NOTES ON LANDMARKS : if the number of landmarks in each file is the same, then " << std::endl;
//   std::cout << "  we assume the landmark correspondence is one to one.  otherwise, use grouped ICP " << std::endl;
//   std::cout << " This implementation is NOT doing exact landmark matching, but it is close to that, in general " << std::endl;

//   std::cout << "  -r/--true-field = STRING     True field filename - default: not used" << std::endl;
   std::cout << "  -a/--affineTransform       restrict transformation to affine" << std::endl;
///   std::cout << "  -t/--gradient-type=UINT    Type of gradient used for computing the demons force"<<std::endl
//            <<"                             (0 is symmetrized, 1 is fixed image, 2 is moving image) - default: "<<defargs.gradientType << std::endl;
//   std::cout<<"  -v/--verbose(=UINT)        Verbosity - default: "<<defargs.verbosity<<"; without argurment: 1" << std::endl;


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

void parseOpts (int argc, char **argv, struct arguments & args)
{
  //itk::Systemfilesystem::path progpath(argv[0]);
   const std::string progname( "FFD Registration/Normalization" );

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
        args.fixedImageFile = optarg;
        break;

      case 'm':
        if (! optarg) display_usage(progname);
        args.movingImageFile = optarg;
        break;

      case 'j':
         if (! optarg) display_usage(progname);
         args.fixedImageFile2 = optarg;
         break;

      case 'k':
         if (! optarg) display_usage(progname);
         args.movingImageFile2 = optarg;
         break;

      case 'y':
         if (! optarg) display_usage(progname);
         args.fixedLandmarksFile = optarg;
         break;

      case 'z':
         if (! optarg) display_usage(progname);
         args.movingLandmarksFile = optarg;
         break;

      case 'w':
         if (! optarg) display_usage(progname);
         args.weightImageFile = optarg;
         break;

      case 'c':
         if (! optarg) display_usage(progname);
         args.whichMetric[0] = atoi(optarg);
         break;

      case 'b':
         if (! optarg) display_usage(progname);
         args.interpolateNimages = atoi(optarg);
         break;

      case 'd':
         if (! optarg) display_usage(progname);
         args.whichMetric[1] = atoi(optarg);
         break;

      case 'o':
         if (! optarg) display_usage(progname);
         args.outputImageFile = optarg;
         break;

      case 'O':
         if (! optarg) display_usage( progname ); //args.outputFieldFile = "CHANGETHISSTRING";
         args.outputFieldFile = optarg;
         break;

      case 'r':
         if (! optarg) display_usage(progname);
         else args.trueFieldFile = optarg;
         break;

      case 'i':
         if (! optarg) display_usage(progname);
         args.numIterations = parseUIntVector(std::string(optarg));
         break;

      case 's':
         if (! optarg) display_usage(progname);
         args.sigmaDef = atof(optarg);
         if ( args.sigmaDef<0.5 and args.sigmaDef>0.0 )
         {
            std::cout<<"Sigma is too small (min=0.5). We set it to 0.0 (no smoothing)."
                     <<std::endl << std::endl;
            args.sigmaDef = 0.0;
         }
         break;

      case 'g':
         if (! optarg) display_usage(progname);
         args.sigmaUp = atof(optarg);
         if ( args.sigmaUp<0.5 and args.sigmaUp>0.0 )
         {
             std::cout<<"Sigma is too small (min=0.5). We set it to 0.0 (no smoothing)."
                 <<std::endl << std::endl;
             args.sigmaUp = 0.0;
         }
         break;

      case 't':
         if (! optarg) display_usage(progname);
         args.useTensors = static_cast<unsigned int>( atoi(optarg) );
         break;

      case 'e':
         args.useHistogramMatching = true;
         break;

      case 'v':
         if (! optarg) args.verbosity++;
         else args.verbosity = static_cast<unsigned int>( atoi(optarg) );
         break;

      case 'p':
         if (! optarg) display_usage( progname );
         else args.prolificacy = static_cast<bool>( atoi(optarg) );
         break;

      case 'a':
         if (! optarg) display_usage( progname );
         else args.affineTransform = static_cast<bool>( atoi(optarg) );
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
      case 'J':
         if (! optarg) display_usage( progname );
         args.useJacobianRegularization = static_cast<bool>( atoi( optarg ) );
         break;
      case 'X':
         if (! optarg) display_usage( progname );
         args.maxJacobian = atof( optarg );
         break;
      case 'N':
         if (! optarg) display_usage( progname );
         args.minJacobian = atof( optarg );
         break;
      case 'I':
         if (! optarg) display_usage( progname );
         args.interiorPenaltyWeight = atof( optarg );
         break;
      case 'D':
         if (! optarg) display_usage( progname );
         args.directionality = optarg;
         break;
      case 'Q':
         if (! optarg) display_usage( progname );
         args.inputControlPointLatticeFile = optarg;
         break;
      case 'P':
         if (! optarg) display_usage( progname );
         args.outputPerformanceFile = optarg;
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


      case 'h': /* fall-through is intentional */
      case '?':   /* fall-through is intentional */
      default:
         display_usage(progname);
         break;
  }
   }

   if ( args.outputFieldFile=="CHANGETHISSTRING" )
   {

      unsigned int pos = args.outputImageFile.find(".");
      if ( pos < args.outputFieldFile.size() )
      {
         args.outputFieldFile = args.outputImageFile;
         args.outputFieldFile.replace(pos, args.outputFieldFile.size(), "-field.mha");
      }
      else
      {
         args.outputFieldFile = args.outputImageFile + "-field.mha";
      }

   }

   if ( args.numIterations.empty() )
   {
      // set a default number of iterations per level
      args.numIterations = std::vector<unsigned int>(3, defiter[0]);
   }
}

template <unsigned int ImageDimension>
int RegisterImages( struct arguments & args )
  {

  std::cout << "Starting RegisterImagesFFD with the following arguments: " << std::endl;
  std::cout << args << std::endl << std::endl;

  typedef float                                                           RealType;

  typedef itk::Image<RealType, ImageDimension>                            ImageType;
  typedef itk::Image<RealType, ImageDimension>                            RealImageType;

  // Read the images
  if ( args.fixedImageFile.empty() || args.movingImageFile.empty() )
    {
    std::cerr << "Error:  One or more image files not specified. " << std::endl;
    exit( 1 );
    }

  // Set up the registration filter
  std::cout << "Set up the registration filter (ImageDimension = " << ImageDimension << ")." << std::endl;

  typedef itk::FFDRegistrationFilter<RealImageType,
    RealImageType, RealImageType>                                         RegistrationFilterType;
  typename RegistrationFilterType::Pointer registrationFilter
    = RegistrationFilterType::New();

  typedef typename RegistrationFilterType::DeformationFieldType           DeformationFieldType;

  typedef itk::ImageFileReader<RealImageType>                             ImageReaderType;

  typename ImageReaderType::Pointer fixedImage = ImageReaderType::New();
  fixedImage->SetFileName( args.fixedImageFile );
  fixedImage->Update();

  std::cout << "Fixed Image: " << std::endl;
  std::cout << "  origin : " << fixedImage->GetOutput()->GetOrigin() << std::endl;
  std::cout << "  spacing: " << fixedImage->GetOutput()->GetSpacing() << std::endl;
  std::cout << "  size   : " << fixedImage->GetOutput()->GetLargestPossibleRegion().GetSize() << std::endl;

  typename ImageReaderType::Pointer movingImage = ImageReaderType::New();
  movingImage->SetFileName( args.movingImageFile );
  movingImage->Update();

  std::cout << "Moving Image: " << std::endl;
  std::cout << "  origin : " << movingImage->GetOutput()->GetOrigin() << std::endl;
  std::cout << "  spacing: " << movingImage->GetOutput()->GetSpacing() << std::endl;
  std::cout << "  size   : " << movingImage->GetOutput()->GetLargestPossibleRegion().GetSize() << std::endl;


  // Rescale the image intensities so that they fall between 0 and 1
  //const RealType desiredMinimum =  0.0;
  //const RealType desiredMaximum =  1.0;
  //
  //typedef itk::RescaleIntensityImageFilter<RealImageType, RealImageType>      RescalerType;
  //typename RescalerType::Pointer movingrescalefilter = RescalerType::New();
  //typename RescalerType::Pointer fixedrescalefilter = RescalerType::New();
  //
  //movingrescalefilter->SetInput( movingImage->GetOutput() );
  //fixedrescalefilter->SetInput( fixedImage->GetOutput() );
  //
  //movingrescalefilter->SetOutputMinimum( desiredMinimum );
  //movingrescalefilter->SetOutputMaximum( desiredMaximum );
  //movingrescalefilter->UpdateLargestPossibleRegion();
  //fixedrescalefilter->SetOutputMinimum( desiredMinimum );
  //fixedrescalefilter->SetOutputMaximum( desiredMaximum );
  //fixedrescalefilter->UpdateLargestPossibleRegion();

  if ( args.useHistogramMatching )
    {
    std::cout << "Histogram match the images." << std::endl;

    // Histogram match the images
    typedef itk::HistogramMatchingImageFilter<RealImageType, RealImageType> HEFilterType;
    typename HEFilterType::Pointer IntensityEqualizeFilter = HEFilterType::New();

//    IntensityEqualizeFilter->SetReferenceImage( fixedrescalefilter->GetOutput() );
//    IntensityEqualizeFilter->SetInput( movingrescalefilter->GetOutput() );
    IntensityEqualizeFilter->SetReferenceImage( fixedImage->GetOutput() );
    IntensityEqualizeFilter->SetInput( movingImage->GetOutput() );
    IntensityEqualizeFilter->SetNumberOfHistogramLevels( 255 );
    IntensityEqualizeFilter->SetNumberOfMatchPoints( 12 );
    IntensityEqualizeFilter->ThresholdAtMeanIntensityOn();
    IntensityEqualizeFilter->Update();

//    registrationFilter->SetInput( 0, fixedrescalefilter->GetOutput() );
    registrationFilter->SetInput( 0, fixedImage->GetOutput() );
    registrationFilter->SetInput( 1, IntensityEqualizeFilter->GetOutput() );
    }
  else
    {
//    registrationFilter->SetInput( 0, fixedrescalefilter->GetOutput() );
//    registrationFilter->SetInput( 1, movingrescalefilter->GetOutput() );
    registrationFilter->SetInput( 0, fixedImage->GetOutput() );
    registrationFilter->SetInput( 1, movingImage->GetOutput() );
    }

  if ( !args.fixedImageFile2.empty() && !args.movingImageFile2.empty() )
    {
    typename ImageReaderType::Pointer fixedImage = ImageReaderType::New();
    fixedImage->SetFileName( args.fixedImageFile2 );
    fixedImage->Update();


    std::cout << "Fixed Image 2: " << std::endl;
    std::cout << "  origin : " << fixedImage->GetOutput()->GetOrigin() << std::endl;
    std::cout << "  spacing: " << fixedImage->GetOutput()->GetSpacing() << std::endl;
    std::cout << "  size   : " << fixedImage->GetOutput()->GetLargestPossibleRegion().GetSize() << std::endl;

    typename ImageReaderType::Pointer movingImage = ImageReaderType::New();
    movingImage->SetFileName( args.movingImageFile2 );
    movingImage->Update();

    std::cout << "Moving Image 2: " << std::endl;
    std::cout << "  origin : " << movingImage->GetOutput()->GetOrigin() << std::endl;
    std::cout << "  spacing: " << movingImage->GetOutput()->GetSpacing() << std::endl;
    std::cout << "  size   : " << movingImage->GetOutput()->GetLargestPossibleRegion().GetSize() << std::endl;

    // Rescale the image intensities so that they fall between 0 and 1
    const RealType desiredMinimum =  0.0;
    const RealType desiredMaximum =  1.0;

    typedef itk::RescaleIntensityImageFilter<RealImageType, RealImageType>      RescalerType;
    typename RescalerType::Pointer movingrescalefilter = RescalerType::New();
    typename RescalerType::Pointer fixedrescalefilter = RescalerType::New();

    movingrescalefilter->SetInput( movingImage->GetOutput() );
    fixedrescalefilter->SetInput( fixedImage->GetOutput() );

    movingrescalefilter->SetOutputMinimum( desiredMinimum );
    movingrescalefilter->SetOutputMaximum( desiredMaximum );
    movingrescalefilter->UpdateLargestPossibleRegion();
    fixedrescalefilter->SetOutputMinimum( desiredMinimum );
    fixedrescalefilter->SetOutputMaximum( desiredMaximum );
    fixedrescalefilter->UpdateLargestPossibleRegion();

    if ( args.useHistogramMatching )
      {
      std::cout << "Histogram match the images." << std::endl;

      // Histogram match the images
      typedef itk::HistogramMatchingImageFilter<RealImageType, RealImageType> HEFilterType;
      typename HEFilterType::Pointer IntensityEqualizeFilter = HEFilterType::New();

      IntensityEqualizeFilter->SetReferenceImage( fixedrescalefilter->GetOutput() );
      IntensityEqualizeFilter->SetInput( movingrescalefilter->GetOutput() );
      IntensityEqualizeFilter->SetNumberOfHistogramLevels( 255 );
      IntensityEqualizeFilter->SetNumberOfMatchPoints( 12 );
      IntensityEqualizeFilter->ThresholdAtMeanIntensityOn();
      IntensityEqualizeFilter->Update();

      registrationFilter->SetInput( 2, fixedrescalefilter->GetOutput() );
      registrationFilter->SetInput( 3, IntensityEqualizeFilter->GetOutput() );
      }
    else
      {
      registrationFilter->SetInput( 2, fixedrescalefilter->GetOutput() );
      registrationFilter->SetInput( 3, movingrescalefilter->GetOutput() );
      }
    }

  // Set up the similarity metrics
  std::cout << "Set up image similarity metric." << std::endl;

  typedef itk::AvantsPDEDeformableRegistrationFunction<RealImageType,
    RealImageType, DeformationFieldType> MetricType;
  typename MetricType::Pointer metric[2];
  typename MetricType::RadiusType radius;
  bool maximizeMetric[2];

  unsigned int numberOfMISamples = 7500;
  unsigned int numberOfHistogramBins = 64;
  if ( fixedImage->GetOutput()->GetLargestPossibleRegion().GetSize()[0] < 80 || ImageDimension == 2 )
    {
    numberOfHistogramBins = 32;
    }
  else if ( fixedImage->GetOutput()->GetLargestPossibleRegion().GetSize()[0] > 256 )
    {
    numberOfHistogramBins = 128;
    }

  for ( unsigned int i = 0; i < 2; i++ )
    {
    if ( args.whichMetric[i] > 8 )
      {
      continue;
      }
    switch( args.whichMetric[i] )
      {
      case 0:
        {
        typedef itk::MeanSquareRegistrationFunction
          <RealImageType, RealImageType, DeformationFieldType> Metric0Type;
        typename Metric0Type::Pointer metric0 = Metric0Type::New();
        metric0->SetIntensityDifferenceThreshold( 1e-3 );
        metric0->SetRobust( true );
        metric0->SetSymmetric( false );
        metric0->SetNormalizeGradient( false );
        metric0->SetGradientStep( 1e6 );
        maximizeMetric[i] = true;
        radius.Fill( 1 );
        metric[i] = metric0;
        break;
        }
      case 1:
        {
        typedef itk::ProbabilisticRegistrationFunction
          <RealImageType, RealImageType, DeformationFieldType> Metric1Type;
        typename Metric1Type::Pointer metric1 = Metric1Type::New();
        metric1->SetFullyRobust( true );
        metric1->SetNormalizeGradient( false );
        metric1->SetGradientStep( 1e6 );
        radius.Fill( 2 );
        maximizeMetric[i] = true;
        metric[i] = metric1;
        break;
        }
      case 2:
        {
        typedef itk::AvantsMutualInformationRegistrationFunction
          <RealImageType, RealImageType, DeformationFieldType> Metric2Type;
        typename Metric2Type::Pointer metric2 = Metric2Type::New();
        metric2->SetNumberOfSpatialSamples( numberOfMISamples );
        metric2->SetNumberOfHistogramBins( numberOfHistogramBins );
        metric2->SetNormalizeGradient( true );
        metric2->SetGradientStep( 1e6 );
        radius.Fill( 1 );
        maximizeMetric[i] = true;
        metric[i] = metric2;
        break;
        }
      case 3:
        {
        typedef itk::AvantsMutualInformationRegistrationFunction
          <RealImageType, RealImageType, DeformationFieldType> Metric3Type;
        typename Metric3Type::Pointer metric3 = Metric3Type::New();
        metric3->SetNumberOfSpatialSamples( numberOfMISamples );
        metric3->SetNumberOfHistogramBins( numberOfHistogramBins );
        metric3->SetNormalizeGradient( false );
        metric3->SetGradientStep( 1e5 );
        radius.Fill( 1 );
        maximizeMetric[i] = true;
        metric[i] = metric3;
        break;
        }
      case 4:
        {
        typedef itk::ProbabilisticRegistrationFunction
          <RealImageType, RealImageType, DeformationFieldType> Metric4Type;
        typename Metric4Type::Pointer metric4 = Metric4Type::New();
        metric4->SetNormalizeGradient( false );
        metric4->SetGradientStep( 1e6 );
        radius.Fill( 5 );
        maximizeMetric[i] = true;
        metric[i] = metric4;
        break;
        }
      case 5:
        {
        typedef itk::ProbabilisticRegistrationFunction
          <RealImageType, RealImageType, DeformationFieldType> Metric5Type;
        typename Metric5Type::Pointer metric5 = Metric5Type::New();
        metric5->SetNormalizeGradient( false );
        metric5->SetGradientStep( 1e6 );
        radius.Fill( 2 );
        maximizeMetric[i] = true;
        metric[i] = metric5;
        break;
        }
      case 6: default:
        {
        typedef itk::RobustOpticalFlow
          <RealImageType, RealImageType, DeformationFieldType> Metric6Type;
        typename Metric6Type::Pointer metric6 = Metric6Type::New();
        radius.Fill( 2 );
        maximizeMetric[i] = false;
        metric[i] = metric6;
        break;
        }
      case 7:
        {
        typedef itk::AvantsMutualInformationRegistrationFunction
          <RealImageType, RealImageType, DeformationFieldType> Metric7Type;
        typename Metric7Type::Pointer metric7 = Metric7Type::New();
        metric7->SetNumberOfSpatialSamples( numberOfMISamples );
        metric7->SetNumberOfHistogramBins( numberOfHistogramBins );
        metric7->SetOpticalFlow( true );
        metric7->SetGradientStep( 1e5 );
        maximizeMetric[i] = true;
        radius.Fill( 1 );
        metric[i] = metric7;
        break;
        }
      case 8:
        {
        typedef itk::SectionMutualInformationRegistrationFunction
          <RealImageType, RealImageType, DeformationFieldType> Metric8Type;
        typename Metric8Type::Pointer metric8 = Metric8Type::New();
        metric8->SetNumberOfSpatialSamples( 7000 );
        metric8->SetNumberOfHistogramBins( 26 );
        metric8->SetOpticalFlow( false );
        metric8->SetNormalizeGradient( true );
        metric8->SetZeroInZ( true );
        metric8->SetGradientStep( 1e2 );
        radius.Fill( 1 );
        maximizeMetric[i] = true;
        metric[i] = metric8;
        break;
        }
      }
    metric[i]->SetRadius( radius );
    metric[i]->SetFixedPointSet( NULL);
    metric[i]->SetMovingPointSet( NULL);
//    maximizeMetric[i] = false;
    }

  if ( !args.weightImageFile.empty() )
    {
    std::cout << "Reading weight image." << std::endl;
    typedef itk::ImageFileReader<typename RegistrationFilterType::WeightImageType> WeightImageReaderType;
    typename WeightImageReaderType::Pointer reader = WeightImageReaderType::New();
    reader->SetFileName( args.weightImageFile );
    try
      {
      reader->Update();
      registrationFilter->SetWeightImage( reader->GetOutput() );
      }
    catch(...)
      {
      }
    }



  // Read landmarks
  if ( !args.movingLandmarksFile.empty() && !args.fixedLandmarksFile.empty() )
    {
    std::cout << "Reading landmarks." << std::endl;
    }

  typename RegistrationFilterType::LandmarkType landmark;
  typename RegistrationFilterType::LandmarkContainerType::Pointer movingLandmarks
    = RegistrationFilterType::LandmarkContainerType::New();
  movingLandmarks->Initialize();
  typename RegistrationFilterType::LandmarkContainerType::Pointer fixedLandmarks
    = RegistrationFilterType::LandmarkContainerType::New();
  fixedLandmarks->Initialize();

  std::fstream strM( args.movingLandmarksFile.c_str(), std::ios::in );
  std::fstream strF( args.fixedLandmarksFile.c_str(), std::ios::in );

  if ( strM.is_open() && strF.is_open() )
    {
    unsigned int i = 0;
    while ( strM >> landmark[0] >> landmark[1] >> landmark[2] >> landmark[3] )
      {
      if ( landmark[3] != 0 )
        {
//          std::cout << "moving landmark: " << landmark << std::endl;
        movingLandmarks->InsertElement( i++, landmark );
        }
      }
    std::cout << i << " moving landmarks read. " << std::endl;
    i = 0;
    while ( strF >> landmark[0] >> landmark[1]  >> landmark[2] >> landmark[3] )
      {
      if ( landmark[3] != 0 )
        {
//          std::cout << "fixed landmark: " << landmark << std::endl;
        fixedLandmarks->InsertElement( i++, landmark );
        }
      }
    std::cout << i << " fixed landmarks read. " << std::endl;
    registrationFilter->SetFixedLandmarkContainer( fixedLandmarks );
    registrationFilter->SetMovingLandmarkContainer( movingLandmarks );
    registrationFilter->SetInitializeWithLandmarks( true );
    registrationFilter->SetLandmarkWeighting( 10000.0 );
    }
  else
    {
    std::cout << "No landmarks read.  One or both of the files \"" << args.movingLandmarksFile
              << "\" or \"" << args.fixedLandmarksFile << "\" does not exist." << std::endl;
    }

  // Read initial control point lattice
  if ( !args.inputControlPointLatticeFile.empty() )
    {
    typedef itk::VectorImageFileReader<RealImageType, DeformationFieldType> ReaderType;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( args.inputControlPointLatticeFile );
    reader->Update();
    registrationFilter->SetInitialControlPointLattice( reader->GetOutput() );
    }


  // Set the rest of the registration filter parameters and produce the output
  typename RegistrationFilterType::ArrayType array;
//  registrationFilter->DebugOn();
  registrationFilter->SetFilePrefix( args.outputImageFile );
  registrationFilter->SetPDEDeformableMetric( metric[0], 0 );
  registrationFilter->SetMetricRadius( metric[0]->GetRadius(), 0 );
  registrationFilter->SetMaximizeMetric( maximizeMetric[0], 0 );
  if ( args.whichMetric[0] != args.whichMetric[1] && args.whichMetric[1] < 99 )
    {
    registrationFilter->SetPDEDeformableMetric( metric[1], 1 );
    registrationFilter->SetMetricRadius( metric[1]->GetRadius(), 1 );
    registrationFilter->SetMaximizeMetric( maximizeMetric[1], 1 );
    }
  registrationFilter->SetWhichGradient( args.useRueckertGradient );
  registrationFilter->SetNumberOfLevels( args.numIterations.size() );
  registrationFilter->SetProlificacy( args.prolificacy );
  registrationFilter->SetEmploySteepestDescent( args.employSteepestDescent );
  std::vector<unsigned int> directionality = parseUIntVector( args.directionality );
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    array[i] = directionality[i];
    }
  registrationFilter->SetDirectionality( array );
  for ( unsigned int i = 0; i < args.numIterations.size(); i++ )
    {
    registrationFilter->SetMaximumNumberOfIterations( args.numIterations[i], i );
    }
  typename RegistrationFilterType::ArrayType fixedImageShrinkFactors;
  typename RegistrationFilterType::ArrayType movingImageShrinkFactors;
  for ( unsigned int i = 0; i < ImageDimension; i++ )
    {
    RealType fixedSize = static_cast<RealType>(
      fixedImage->GetOutput()->GetLargestPossibleRegion().GetSize()[i] );
    RealType movingSize = static_cast<RealType>(
      movingImage->GetOutput()->GetLargestPossibleRegion().GetSize()[i] );
    RealType fixedFactor = vnl_math_min( pow( 2, args.numIterations.size()-1 ), fixedSize/32.0 );
    RealType movingFactor = vnl_math_min( pow( 2, args.numIterations.size()-1 ), movingSize/32.0 );

    fixedImageShrinkFactors[i] = vnl_math_max( 1u, static_cast<unsigned int>( fixedFactor ) );
    movingImageShrinkFactors[i] = vnl_math_max( 1u, static_cast<unsigned int>( movingFactor ) );
    }
  registrationFilter->SetFixedImageShrinkFactors( fixedImageShrinkFactors );
  registrationFilter->SetMovingImageShrinkFactors( movingImageShrinkFactors );
  std::cout << "Fixed image shrink factors: " << fixedImageShrinkFactors << std::endl;
  std::cout << "Moving image shrink factors: " << movingImageShrinkFactors << std::endl;

  if ( args.affineTransform )
    {
    registrationFilter->SetEnforceDiffeomorphism( false );
    registrationFilter->SetSplineOrder( 1 );
    registrationFilter->SetLineSearchMaximumIterations( 0 );
    registrationFilter->SetNumberOfAffineNeighborhoodSamplesPerIteration( 1 );
    registrationFilter->SetAffineNeighborhoodRadius( 2 );

    typename RegistrationFilterType::ArrayType array;
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      array[i] = fixedImage->GetOutput()->GetLargestPossibleRegion().GetSize()[i];
      }
    registrationFilter->SetMeshResolution( array );
    registrationFilter->SetDoubleMeshResolutionAtEachLevel( false );
    }
  else
    {
    registrationFilter->SetEnforceDiffeomorphism( args.useJacobianRegularization );
    registrationFilter->SetMinimumJacobian( args.minJacobian );
    registrationFilter->SetMinimumJacobian( args.maxJacobian );
    registrationFilter->SetInteriorPenaltyParameter( args.interiorPenaltyWeight );
    registrationFilter->SetSplineOrder( args.bsplineOrder );
    registrationFilter->SetDoubleMeshResolutionAtEachLevel( true );
    if ( args.employLineSearch )
      {
      registrationFilter->SetLineSearchMaximumIterations( args.maxNumberOfLineSearchIterations );
      }
    else
      {
      registrationFilter->SetLineSearchMaximumIterations( 0 );
      }
    typename RegistrationFilterType::ArrayType array;
    std::vector<unsigned int> meshResolution = parseUIntVector( args.bsplineMeshResolution );
    for ( unsigned int i = 0; i < ImageDimension; i++ )
      {
      array[i] = meshResolution[i];
      }
    registrationFilter->SetMeshResolution( array );
    }

  itk::TimeProbe timer;
  timer.Start();
  registrationFilter->Update();
  timer.Stop();
  std::cout << "Registration filter run time = " << timer.GetMeanTime() << std::endl;

  // Write the outputs

  std::cout << "Writing the output." << std::endl;

  //std::cout << "Total landmark error = " << registrationFilter->CalculateLandmarkError() << std::endl;

  std::string file;

/*
  if ( !args.outputImageFile.empty() )
    {
    typedef itk::WarpImageFilter<RealImageType,
                                 RealImageType,
                                 DeformationFieldType> WarperType;
    typename WarperType::Pointer warper = WarperType::New();

    warper->SetInput( movingImage->GetOutput() );
    warper->SetDeformationField( registrationFilter->GetDeformationField() );
    warper->SetInterpolator( registrationFilter->GetImageInterpolator() );
    warper->SetOutputSpacing( fixedImage->GetOutput()->GetSpacing() );
    warper->SetOutputOrigin( fixedImage->GetOutput()->GetOrigin() );
    warper->SetEdgePaddingValue( 1.0 );
    warper->Update();

    file = std::string( args.outputImageFile );
    typedef itk::ImageFileWriter<RealImageType> ImageWriterType;
    typename ImageWriterType::Pointer writer = ImageWriterType::New();
    writer->SetFileName( file.c_str() );
    writer->SetInput( warper->GetOutput() );
    writer->Update();
    }

  typedef itk::VectorImageFileWriter<DeformationFieldType, RealImageType> DeformationFieldWriterType;


  if ( !args.outputFieldFile.empty() )
    {
    file = std::string( args.outputFieldFile );
    typename DeformationFieldWriterType::Pointer dfwriter = DeformationFieldWriterType::New();
    dfwriter->SetFileName( file.c_str() );
    dfwriter->SetInput( registrationFilter->GetDeformationField() );
    dfwriter->Update();
    }

//  if ( !args.outputPerformanceFile.empty() )
    {
    std::string::size_type Pos = args.outputFieldFile.rfind( "." );
    std::string extension( args.outputFieldFile, Pos, args.outputFieldFile.length()-1 );
    std::string filename( args.outputFieldFile, 0, Pos );
    filename += std::string( "PerformanceFile.txt" );

//    ofstream str( args.outputPerformanceFile.c_str() );
    ofstream str( filename.c_str() );
    str << std::setw( 20 ) << "Iteration Number" << std::setw( 10 ) << "Level" << std::setw( 15 ) << "Energy Value" << std::setw( 25 ) << "Gradient Comp. Time" << std::endl;
    str.setf( std::ios::showpoint );
    for ( unsigned int i = 0; i < registrationFilter->GetLevelNumbers().size(); i++ )
      {
      str << std::setw( 20 ) << i
          << std::setw( 10 ) << registrationFilter->GetLevelNumbers()[i]
          << std::setw( 15 ) << registrationFilter->GetEnergyValues()[i]
          << std::setw( 25 ) << std::setprecision( 4 )
          << registrationFilter->GetGradientComputationTimes()[i]
          << std::endl;
      }
    }
*/

  {
  std::string filename = std::string( args.outputImageFile )
    + std::string( ".nii.gz" );


  typedef itk::BSplineControlPointImageFilter
    <DeformationFieldType, DeformationFieldType> BSplineControlPointsFilterType;
  typename BSplineControlPointsFilterType::Pointer bspliner 
    = BSplineControlPointsFilterType::New();
  typename BSplineControlPointsFilterType::ArrayType close;

  close.Fill( false );
  bspliner->SetSplineOrder( registrationFilter->GetSplineOrder() );
  bspliner->SetCloseDimension( close );
  bspliner->SetInput( registrationFilter->GetTotalDeformationFieldControlPoints() );
  bspliner->SetOrigin( fixedImage->GetOutput()->GetOrigin() );
  bspliner->SetSize( fixedImage->GetOutput()->GetLargestPossibleRegion().GetSize() );
  bspliner->SetSpacing( fixedImage->GetOutput()->GetSpacing() );
  bspliner->Update();

  typedef itk::WarpImageFilter<RealImageType, RealImageType,
                          DeformationFieldType> WarperType;
  typename WarperType::Pointer warper = WarperType::New();

  warper->SetInput( movingImage->GetOutput() );
  warper->SetDeformationField( bspliner->GetOutput() );
//  warper->SetInterpolator( this->m_ImageInterpolator );
  warper->SetOutputSpacing( fixedImage->GetOutput()->GetSpacing() );
  warper->SetOutputOrigin( fixedImage->GetOutput()->GetOrigin() );
  warper->SetEdgePaddingValue( 0.0 );
  warper->Update();

  typedef itk::ImageFileWriter<RealImageType> ImageWriterType;
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInput( warper->GetOutput() );
  writer->Update();
  }

  typedef itk::VectorImageFileWriter<DeformationFieldType, RealImageType> DeformationFieldWriterType;

  {
  itk::OStringStream buf;
  buf << args.bsplineOrder;

  std::string filename = std::string( args.outputImageFile )
    + std::string( "CPLattice_" ) + buf.str()
    + std::string( ".nii.gz" );

  typename DeformationFieldWriterType::Pointer cpwriter = DeformationFieldWriterType::New();
  cpwriter->SetFileName( filename.c_str() );
  cpwriter->SetInput( registrationFilter->GetTotalDeformationFieldControlPoints() );
  cpwriter->Update();
  }

  {
  std::string filename = std::string( args.outputImageFile )
    + std::string( "Performance.txt" );

  ofstream str( filename.c_str() );
  str << std::setw( 20 ) << "Iteration Number" << std::setw( 10 ) << "Level" << std::setw( 15 ) << std::setprecision( 8 ) << "Energy Value" << std::setw( 25 ) << "Gradient Comp. Time" << std::endl;
  str.setf( std::ios::showpoint );
  for ( unsigned int i = 0; i < registrationFilter->GetLevelNumbers().size(); i++ )
    {
    str << std::setw( 20 ) << i
        << std::setw( 10 ) << registrationFilter->GetLevelNumbers()[i]
        << std::setw( 15 ) << std::setprecision( 8 )
        << registrationFilter->GetEnergyValues()[i]
        << std::setw( 25 ) << std::setprecision( 4 )
        << registrationFilter->GetGradientComputationTimes()[i]
        << " (" << registrationFilter->GetNumberOfGradientPoints()[i] << ")"
        << std::endl;
    }
  }
  return EXIT_SUCCESS;  
}

int main( int argc, char *argv[] )
{
  if ( argc == 1 )
    {
    display_usage( "RegisterImagesFFD" );
    exit( 0 );
    }

  struct arguments args;
  parseOpts( argc, argv, args );

  // Get the image dimension
  itk::ImageIOBase::Pointer imageIO =
     itk::ImageIOFactory::CreateImageIO(
        args.fixedImageFile.c_str(), itk::ImageIOFactory::ReadMode);
  imageIO->SetFileName(args.fixedImageFile.c_str());
  imageIO->ReadImageInformation();

  switch ( imageIO->GetNumberOfDimensions() )
    {
    case 2:
      RegisterImages<2>( args );
      break;
    case 3:
      RegisterImages<3>( args );
      break;
    default:
      std::cerr << "Unsupported dimension" << std::endl;
      exit( EXIT_FAILURE );
    }

  return EXIT_SUCCESS;
}

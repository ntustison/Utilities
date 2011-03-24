#include "imagewriting.h"
#include "itkimagereader.h"
#include "itkErodedMaskPointSetToImageFilter.h"
#include "itkESMMeanSquaresImageToImageMetric.h"
#include "itkESMOptimizer.h"
#include "itkESMDogLegOptimizer.h"
#include "itkESMRigid2DTransform.h"
#include "itkOptimizerCommandIterationUpdate.h"
#include "normxcorr2.h"
#include "parsingtools.h"
#include "timer.h"

#include <itkFRPROptimizer.h>
#include <vnl/algo/vnl_lbfgs.h>
#include <vnl/algo/vnl_lbfgsb.h>

#include <boost/filesystem/operations.hpp>
#include <boost/random.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/tuple/tuple.hpp>

#include <QFileInfo>
#include <QStringList>
#include <QTextStream>

#include <numeric>
#include <string>
#include <iostream>
#include <fstream>

#include "mzmain_noqapp.h"
#include "qxtcommandoptions.h"

namespace mz {

struct arguments
{
   QString fileIn;
   int first;
   int last;
   int skip;
   double angleScale;

   arguments()
      :fileIn("")
      ,first(0)
      ,last(-1)
      ,skip(0)
      ,angleScale(1.0/0.5)
   {
   }

   friend QTextStream& operator<< (QTextStream& o, const arguments& args)
   {
      std::ostringstream osstr;

      return o
         <<"Arguments structure:"<<endl
         <<"  Input file: "<<args.fileIn<<endl
         <<"  First frame: "<<args.first<<endl
         <<"  Last frame: "<<args.last<<endl
         <<"  Skip frames: "<<args.skip<<endl
         <<"  Angle scale: "<<args.angleScale<<endl;
  }
};


/* Display program usage, and exit.
 */
void display_usage( const std::string progname )
{
   struct arguments defargs = arguments();

   std::cout<<std::endl;
   std::cout<<progname<<" - test"<<std::endl;
   std::cout<<"Usage: "<<progname<<" [OPTION...]"<<std::endl;

   std::cout<<"  -i movie.mkt       input mkt file - MANDATORY"<<std::endl;
   std::cout<<"  -n _first_x_last_  first and last frames - default: [0 numFrames-1]"<<std::endl;
   std::cout<<"  -s skip            skip n frames between each input frame - default: "<<defargs.skip<<std::endl;
   std::cout<<"  -a scale           Scaling factor for angles - default: "<<defargs.angleScale<<std::endl;

   std::cout<<std::endl;
   std::cout<<"   Copyright (c) 2008 Mauna Kea Technologies."<<std::endl;
   std::cout<<"   Code: Tom Vercauteren."<<std::endl;
   std::cout<<"   Report bugs to <tom@maunakeatech.com>."<<std::endl;
   std::cout<<std::endl;

   exit( EXIT_FAILURE );
}


void parseOpts (QStringList mzQArgs, QxtCommandOptions & options, struct arguments & args)
{
   boost::filesystem::path progpath(mzQArgs[0].toStdString());
   const std::string progname(progpath.leaf());

   // Default values.
   args = arguments();

   if (mzQArgs.size() == 1)
   {
      display_usage(progname);
   }

  // Explanatory message
   const QString helpMessage = QString::fromStdString(progname) + QString(":   Reconstruction test");
   options.addSection(helpMessage);

   // Available options
   options.addSection(" Options:");

   options.add( "help", "show this help text",
                QxtCommandOptions::NoValue );
   options.alias( "help", "h" );

   options.add( "input", "input file",
                QxtCommandOptions::ValueRequired );
   options.alias( "input", "i" );

   options.add( "frameindex", "fixed and moving frame index",
                QxtCommandOptions::ValueRequired );
   options.alias( "frameindex", "n" );

   options.add( "skippedframes", "number of skipped frame",
                QxtCommandOptions::ValueRequired );
   options.alias( "skippedframes", "s" );

   options.add( "scalingfactor", "scaling factor",
                QxtCommandOptions::ValueRequired );
   options.alias( "scalingfactor", "a" );

   // Parse the command line
   options.parse( mzQArgs );

  // Check whether "help" was required or if any unrecognized
   // options were processed
   if( options.count("help") || options.showUnrecognizedWarning() )
   {
      options.showUsage();
   }

   // Useful streams
   QTextStream qcout( stdout, QIODevice::WriteOnly );
   QTextStream qcerr( stderr, QIODevice::WriteOnly );

   // Retrieve the options
   if ( options.count("input") )
   {
      args.fileIn = options.value("input").toString();
      if (args.fileIn.isEmpty())
      {
            printf ("Error, you need to provide a -i file.mkt parameter.\n");
            exit( EXIT_FAILURE );
      }
      //qcout << "getStr " << getStr << endl;
   }

   QString frame_index("");
   if ( options.count("frameindex") )
   {
      frame_index = options.value("frameindex").toString();
      if (frame_index.isEmpty()) {
         printf ("Error, you need to provide a fixed and moving frame indexes.\n");
         exit (EXIT_FAILURE);
      }

      try {
         boost::tuples::tie(args.first, args.last) = parse_int_pair( frame_index.toStdString() );
      } catch (std::exception & e) {
         std::cout<<"Error: "<<e.what()<<std::endl;
         exit (EXIT_FAILURE);
      }
   }

   QString skipstring("");
   if ( options.count("skippedframes") )
   {
      skipstring = options.value("skippedframes").toString();
      if (skipstring.isEmpty()) {
         printf ("Error, you need to provide a number of skipped frames.\n");
         exit (EXIT_FAILURE);
      }

      try {
         args.skip = atoi_check( skipstring.toStdString().c_str() );
      } catch (std::exception & e) {
         std::cout<<"Error: "<<e.what()<<std::endl;
         exit( EXIT_FAILURE );
      }
   }

   QString scaling("");
   if ( options.count("scalingfactor") )
   {
      scaling = options.value("scaling").toString();
      if (scaling.isEmpty()) {
            qcerr<<endl;
            qcerr<<"Error, you need to provide a scaling factor."<<endl;
            display_usage(progname);
         }
         args.angleScale = atof_check(scaling.toStdString().c_str());
   }
}


template <class ImageType, class MaskImageType>
void maskAndNormalize( const double & mean,
                       const MaskImageType * const mask,
                       const ImageType * const image,
                       ImageType * const normimage )
{
   typedef typename ImageType::PixelType      PixelType;
   typedef typename MaskImageType::PixelType  MaskPixelType;

   const unsigned int numPix = image->GetLargestPossibleRegion().GetNumberOfPixels();

   const PixelType * imgptr = image->GetBufferPointer();
   const PixelType * const endptr = imgptr + numPix;

   const MaskPixelType * maskptr = mask->GetBufferPointer();

   PixelType * nimgptr = normimage->GetBufferPointer();

   while ( imgptr!=endptr )
   {
      if (*maskptr++)
      {
         *nimgptr++ = (*imgptr++) - mean;
      }
      else
      {
         *nimgptr++ = 0;
         ++imgptr;
      }
   }
}



template <class ImageType, class IndexVectorType, class MaskType>
double registerESM( const ImageType * const fixim,
                    const ImageType * const movim,
                    const IndexVectorType * const fibidx,
                    const MaskType * const mask,
                    const double & tx, const double & ty,
                    const arguments & args )
{
   // Set up metric
   const itk::ImageRegion<2> region = fixim->GetLargestPossibleRegion();

   const double lx = region.GetSize(0)*fixim->GetSpacing()[0];
   const double ly = region.GetSize(1)*fixim->GetSpacing()[1];

   const double diagonallength = sqrt( lx*lx + ly*ly );

   typename itk::ESMRigid2DTransform<>::Pointer transform =
      itk::ESMRigid2DTransform<>::New();

   typename itk::ESMRigid2DTransform<>::PointType center;
   center[0] = fixim->GetOrigin()[0] + lx/2.0;
   center[1] = fixim->GetOrigin()[1] + ly/2.0;

   transform->SetCenter( center );

   typename itk::ESMMeanSquaresImageToImageMetric<ImageType,ImageType>::Pointer metric =
      itk::ESMMeanSquaresImageToImageMetric<ImageType,ImageType>::New();
   metric->SetMovingImage( movim );
   metric->SetFixedImage( fixim );
   metric->SetESMTransform( transform );
   metric->SetMovingImageMask( mask );
   metric->SetFixedImageRegion( fixim->GetLargestPossibleRegion() );
   metric->SetFixedImageIndexes( fibidx->CastToSTLConstContainer() );
   metric->Initialize();

   // Initial parameters
   itk::Array<double> initRigidParam(3);
   initRigidParam(0) = 0.0;
   initRigidParam(1) = tx;
   initRigidParam(2) = ty;

   // Set up scales
   itk::Optimizer::ScalesType scales(3);
   scales[0] = args.angleScale;
   scales[1] = 1.0/diagonallength;
   scales[2] = 1.0/diagonallength;

   bool usescale = false;

   // Create observer
   itk::OptimizerCommandIterationUpdate::Pointer observer =
      itk::OptimizerCommandIterationUpdate::New();

   // Set up optimizer
   //typename itk::ESMOptimizer<ImageType,ImageType>::Pointer optimizer
   //   = itk::ESMOptimizer<ImageType,ImageType>::New();
   typename itk::ESMDogLegOptimizer<ImageType,ImageType>::Pointer optimizer
      = itk::ESMDogLegOptimizer<ImageType,ImageType>::New();

   if (usescale) optimizer->SetScales( scales );
   optimizer->AddObserver( itk::IterationEvent(), observer );
   optimizer->SetESMCostFunction( metric );
   optimizer->SetInitialPosition( initRigidParam );

   optimizer->SetMaximumNumberOfIterations( 10 );

   optimizer->StartOptimization();

   std::cout<<std::endl
            <<"ESM end. end_code="<<optimizer->GetStopConditionDescription()
            <<"  finalCost="<<optimizer->GetOptimalValue()<<std::endl
            <<"  finalParams="<<optimizer->GetOptimalParameters()
            <<std::endl<<std::endl;

   return optimizer->GetOptimalValue();
}


template <class ImageType, class IndexVectorType, class MaskType>
double registerFRPR( const ImageType * const fixim,
                     const ImageType * const movim,
                     const IndexVectorType * const fibidx,
                     const MaskType * const mask,
                     const double & tx, const double & ty,
                     const arguments & args )
{
   // Set up metric
   const itk::ImageRegion<2> region = fixim->GetLargestPossibleRegion();

   const double lx = region.GetSize(0)*fixim->GetSpacing()[0];
   const double ly = region.GetSize(1)*fixim->GetSpacing()[1];

   const double diagonallength = sqrt( lx*lx + ly*ly );

   typename itk::ESMRigid2DTransform<>::Pointer transform =
      itk::ESMRigid2DTransform<>::New();

   typename itk::ESMRigid2DTransform<>::PointType center;
   center[0] = fixim->GetOrigin()[0] + lx/2.0;
   center[1] = fixim->GetOrigin()[1] + ly/2.0;

   transform->SetCenter( center );

   typename itk::ESMMeanSquaresImageToImageMetric<ImageType,ImageType>::Pointer metric =
      itk::ESMMeanSquaresImageToImageMetric<ImageType,ImageType>::New();
   metric->SetMovingImage( movim );
   metric->SetFixedImage( fixim );
   metric->SetESMTransform( transform );
   metric->SetMovingImageMask( mask );
   metric->SetFixedImageRegion( fixim->GetLargestPossibleRegion() );
   metric->SetFixedImageIndexes( fibidx->CastToSTLConstContainer() );
   metric->Initialize();

   // Initial parameters
   itk::Array<double> initRigidParam(3);
   initRigidParam(0) = 0.0;
   initRigidParam(1) = tx;
   initRigidParam(2) = ty;

   // Set up scales
   itk::Optimizer::ScalesType scales(3);
   scales[0] = args.angleScale;
   scales[1] = 1.0/diagonallength;
   scales[2] = 1.0/diagonallength;

   bool usescale = false;

   // Create observer
   itk::OptimizerCommandIterationUpdate::Pointer observer =
      itk::OptimizerCommandIterationUpdate::New();

   // Set up optimizer
   itk::FRPRPatchedOptimizer::Pointer optimizer = itk::FRPRPatchedOptimizer::New();
   //itk::FRPROptimizer::Pointer optimizer = itk::FRPROptimizer::New();

   if (usescale) optimizer->SetScales( scales );
   optimizer->AddObserver( itk::IterationEvent(), observer );
   optimizer->SetCostFunction( metric );
   optimizer->SetInitialPosition( initRigidParam );

   optimizer->SetMaximumIteration( 200 );
   optimizer->SetStepLength( 1.0e-1 );
   optimizer->SetStepTolerance( 1.0e-4 );
   optimizer->SetValueTolerance( 1.0e-6 );

   optimizer->StartOptimization();

   return optimizer->GetValue();
}



template <class ImageType, class IndexVectorType, class MaskType>
double registerLBFGS_vnl( const ImageType * const fixim,
                          const ImageType * const movim,
                          const IndexVectorType * const fibidx,
                          const MaskType * const mask,
                          const double & tx, const double & ty,
                          const arguments & args )
{
   // Set up metric
   const itk::ImageRegion<2> region = fixim->GetLargestPossibleRegion();

   const double lx = region.GetSize(0)*fixim->GetSpacing()[0];
   const double ly = region.GetSize(1)*fixim->GetSpacing()[1];

   const double diagonallength = sqrt( lx*lx + ly*ly );

   typename itk::ESMRigid2DTransform<>::Pointer transform =
      itk::ESMRigid2DTransform<>::New();

   typename itk::ESMRigid2DTransform<>::PointType center;
   center[0] = fixim->GetOrigin()[0] + lx/2.0;
   center[1] = fixim->GetOrigin()[1] + ly/2.0;

   transform->SetCenter( center );

   typename itk::ESMMeanSquaresImageToImageMetric<ImageType,ImageType>::Pointer metric =
      itk::ESMMeanSquaresImageToImageMetric<ImageType,ImageType>::New();
   metric->SetMovingImage( movim );
   metric->SetFixedImage( fixim );
   metric->SetESMTransform( transform );
   metric->SetMovingImageMask( mask );
   metric->SetFixedImageRegion( fixim->GetLargestPossibleRegion() );
   metric->SetFixedImageIndexes( fibidx->CastToSTLConstContainer() );
   metric->Initialize();

   // Initial parameters
   itk::Array<double> initRigidParam(3);
   initRigidParam(0) = 0.0;
   initRigidParam(1) = tx;
   initRigidParam(2) = ty;

   // Set up scales
   itk::Optimizer::ScalesType scales(3);
   scales[0] = args.angleScale;
   scales[1] = 1.0/diagonallength;
   scales[2] = 1.0/diagonallength;

   bool usescale = true;

   // Set up cost function adaptor
   typedef itk::SingleValuedVnlCostFunctionAdaptor CostAdaptorType;

   boost::scoped_ptr<CostAdaptorType> adaptor( new CostAdaptorType(3) );
   adaptor->SetCostFunction( metric );
   if (usescale) adaptor->SetScales(scales);

   // Set up optimizer
   boost::scoped_ptr<vnl_lbfgs> optimizer(  new vnl_lbfgs( *adaptor ) );

   vnl_vector<double> params(3);
   if (usescale)
   {
      params[0] = initRigidParam(0) * scales[0];
      params[1] = initRigidParam(1) * scales[1];
      params[2] = initRigidParam(2) * scales[2];
   }
   else
   {
      params[0] = initRigidParam(0);
      params[1] = initRigidParam(1);
      params[2] = initRigidParam(2);
   }

   optimizer->set_trace( true );
   //optimizer->set_max_function_evals( 1000 );
   //optimizer->set_f_tolerance( 1e-2 ); unused
   //optimizer->set_x_tolerance( 1e-2 ); unused
   optimizer->set_g_tolerance( 1e-1 );
   //optimizer->line_search_accuracy = 0.9;
   optimizer->default_step_length  = 0.3;

   //optimizer->set_verbose( true );
   //optimizer->set_epsilon_function( 1e-2 );
   //optimizer->set_check_derivatives(1);

   bool ok = optimizer->minimize( params );

   itk::Array<double> finalRigidParam(3);
   if (usescale)
   {
      finalRigidParam(0) = params[0] / scales[0];
      finalRigidParam(1) = params[1] / scales[1];
      finalRigidParam(2) = params[2] / scales[2];
   }
   else
   {
      finalRigidParam(0) = params[0];
      finalRigidParam(1) = params[1];
      finalRigidParam(2) = params[2];
   }

   const double finalCost = metric->GetValue(finalRigidParam);

   std::cout<<std::endl
            <<"vnl_lbfgs end. ok="<<ok<<"  end_code="<<optimizer->get_failure_code()
            <<"  end_error="<<optimizer->get_end_error()
            <<"  finalCost="<<finalCost<<std::endl
            <<"  finalParams="<<finalRigidParam
            <<std::endl<<std::endl;

   return finalCost;
}



template <class ImageType, class IndexVectorType, class MaskType>
double registerLBFGSB_vnl( const ImageType * const fixim,
                           const ImageType * const movim,
                           const IndexVectorType * const fibidx,
                           const MaskType * const mask,
                           const double & tx, const double & ty,
                           const arguments & args )
{
   // Set up metric
   const itk::ImageRegion<2> region = fixim->GetLargestPossibleRegion();

   const double lx = region.GetSize(0)*fixim->GetSpacing()[0];
   const double ly = region.GetSize(1)*fixim->GetSpacing()[1];

   const double diagonallength = sqrt( lx*lx + ly*ly );

   typename itk::ESMRigid2DTransform<>::Pointer transform =
      itk::ESMRigid2DTransform<>::New();

   typename itk::ESMRigid2DTransform<>::PointType center;
   center[0] = fixim->GetOrigin()[0] + lx/2.0;
   center[1] = fixim->GetOrigin()[1] + ly/2.0;

   transform->SetCenter( center );

   typename itk::ESMMeanSquaresImageToImageMetric<ImageType,ImageType>::Pointer metric =
      itk::ESMMeanSquaresImageToImageMetric<ImageType,ImageType>::New();
   metric->SetMovingImage( movim );
   metric->SetFixedImage( fixim );
   metric->SetESMTransform( transform );
   metric->SetMovingImageMask( mask );
   metric->SetFixedImageRegion( fixim->GetLargestPossibleRegion() );
   metric->SetFixedImageIndexes( fibidx->CastToSTLConstContainer() );
   metric->Initialize();

   // Initial parameters
   itk::Array<double> initRigidParam(3);
   initRigidParam(0) = 0.0;
   initRigidParam(1) = tx;
   initRigidParam(2) = ty;

   // Set up scales
   itk::Optimizer::ScalesType scales(3);
   scales[0] = args.angleScale;
   scales[1] = 1.0/diagonallength;
   scales[2] = 1.0/diagonallength;

   bool usescale = true;

   // Set up cost function adaptor
   typedef itk::SingleValuedVnlCostFunctionAdaptor CostAdaptorType;

   boost::scoped_ptr<CostAdaptorType> adaptor( new CostAdaptorType(3) );
   adaptor->SetCostFunction( metric );
   if (usescale) adaptor->SetScales(scales);

   // Set up optimizer
   boost::scoped_ptr<vnl_lbfgsb> optimizer(  new vnl_lbfgsb( *adaptor ) );

   vnl_vector<double> params(3);
   if (usescale)
   {
      params[0] = initRigidParam(0) * scales[0];
      params[1] = initRigidParam(1) * scales[1];
      params[2] = initRigidParam(2) * scales[2];
   }
   else
   {
      params[0] = initRigidParam(0);
      params[1] = initRigidParam(1);
      params[2] = initRigidParam(2);
   }

   optimizer->set_trace( true );
   //optimizer->set_verbose( true );
   //optimizer->set_max_function_evals( 1000 );
   //optimizer->set_f_tolerance( 1e-2 ); unused
   //optimizer->set_x_tolerance( 1e-2 ); unused
   //optimizer->set_g_tolerance( 1e-1 ); unused
   //optimizer->line_search_accuracy = 0.9; unused
   optimizer->set_projected_gradient_tolerance( 1e-1 );
   optimizer->set_cost_function_convergence_factor( 1e10 );

   bool ok = optimizer->minimize( params );

   itk::Array<double> finalRigidParam(3);
   if (usescale)
   {
      finalRigidParam(0) = params[0] / scales[0];
      finalRigidParam(1) = params[1] / scales[1];
      finalRigidParam(2) = params[2] / scales[2];
   }
   else
   {
      finalRigidParam(0) = params[0];
      finalRigidParam(1) = params[1];
      finalRigidParam(2) = params[2];
   }

   const double finalCost = metric->GetValue(finalRigidParam);

   std::cout<<std::endl
            <<"vnl_lbfgsb end. ok="<<ok<<"  end_code="<<optimizer->get_failure_code()
            <<"  end_error="<<optimizer->get_end_error()
            <<"  finalCost="<<finalCost<<std::endl
            <<"  finalParams="<<finalRigidParam
            <<std::endl<<std::endl;

   return finalCost;
}



template <class ImageType, class IndexVectorType, class MaskType>
double registerLBFGSB_vnlb( const ImageType * const fixim,
                            const ImageType * const movim,
                            const IndexVectorType * const fibidx,
                            const MaskType * const mask,
                            const double & tx, const double & ty,
                            const arguments & args )
{
   // Set up metric
   const itk::ImageRegion<2> region = fixim->GetLargestPossibleRegion();

   const double lx = region.GetSize(0)*fixim->GetSpacing()[0];
   const double ly = region.GetSize(1)*fixim->GetSpacing()[1];

   const double diagonallength = sqrt( lx*lx + ly*ly );

   typename itk::ESMRigid2DTransform<>::Pointer transform =
      itk::ESMRigid2DTransform<>::New();

   typename itk::ESMRigid2DTransform<>::PointType center;
   center[0] = fixim->GetOrigin()[0] + lx/2.0;
   center[1] = fixim->GetOrigin()[1] + ly/2.0;

   transform->SetCenter( center );

   typename itk::ESMMeanSquaresImageToImageMetric<ImageType,ImageType>::Pointer metric =
      itk::ESMMeanSquaresImageToImageMetric<ImageType,ImageType>::New();
   metric->SetMovingImage( movim );
   metric->SetFixedImage( fixim );
   metric->SetESMTransform( transform );
   metric->SetMovingImageMask( mask );
   metric->SetFixedImageRegion( fixim->GetLargestPossibleRegion() );
   metric->SetFixedImageIndexes( fibidx->CastToSTLConstContainer() );
   metric->Initialize();

   // Initial parameters
   itk::Array<double> initRigidParam(3);
   initRigidParam(0) = 0.0;
   initRigidParam(1) = tx;
   initRigidParam(2) = ty;

   // Set up scales
   itk::Optimizer::ScalesType scales(3);
   scales[0] = args.angleScale;
   scales[1] = 1.0/diagonallength;
   scales[2] = 1.0/diagonallength;

   bool usescale = true;


   // Set up bounds
   double thetabound = vnl_math::pi/10.0;
   double txbound = lx*5.0/8.0;
   double tybound = ly*5.0/8.0;

   if ( usescale )
   {
      thetabound *= scales[0];
      txbound *= scales[1];
      tybound *= scales[2];
   }

   vnl_vector<long> boundtypes(3,2);
   vnl_vector<double> upperbound(3);
   vnl_vector<double> lowerbound(3);

   upperbound[0] =  thetabound;
   lowerbound[0] = -thetabound;

   upperbound[1] =  txbound;
   lowerbound[1] = -txbound;

   upperbound[2] =  tybound;
   lowerbound[2] = -tybound;

    bool usebounds = true;


   // Set up cost function adaptor
   typedef itk::SingleValuedVnlCostFunctionAdaptor CostAdaptorType;

   boost::scoped_ptr<CostAdaptorType> adaptor( new CostAdaptorType(3) );
   adaptor->SetCostFunction( metric );
   if (usescale) adaptor->SetScales(scales);

   // Set up optimizer
   boost::scoped_ptr<vnl_lbfgsb> optimizer(  new vnl_lbfgsb( *adaptor ) );

   vnl_vector<double> params(3);
   if (usescale)
   {
      params[0] = initRigidParam(0) * scales[0];
      params[1] = initRigidParam(1) * scales[1];
      params[2] = initRigidParam(2) * scales[2];
   }
   else
   {
      params[0] = initRigidParam(0);
      params[1] = initRigidParam(1);
      params[2] = initRigidParam(2);
   }

   if ( usebounds )
   {
      optimizer->set_bound_selection( boundtypes );
      optimizer->set_upper_bound( upperbound );
      optimizer->set_lower_bound( lowerbound );
   }

   optimizer->set_trace( true );
   //optimizer->set_verbose( true );
   //optimizer->set_max_function_evals( 1000 );
   //optimizer->set_f_tolerance( 1e-2 ); unused
   //optimizer->set_x_tolerance( 1e-2 ); unused
   //optimizer->set_g_tolerance( 1e-1 ); unused
   //optimizer->line_search_accuracy = 0.9; unused
   optimizer->set_projected_gradient_tolerance( 1e-1 );
   optimizer->set_cost_function_convergence_factor( 1e10 );

   bool ok = optimizer->minimize( params );

   itk::Array<double> finalRigidParam(3);
   if (usescale)
   {
      finalRigidParam(0) = params[0] / scales[0];
      finalRigidParam(1) = params[1] / scales[1];
      finalRigidParam(2) = params[2] / scales[2];
   }
   else
   {
      finalRigidParam(0) = params[0];
      finalRigidParam(1) = params[1];
      finalRigidParam(2) = params[2];
   }

   const double finalCost = metric->GetValue(finalRigidParam);

   std::cout<<std::endl
            <<"vnl_lbfgsb b end. ok="<<ok<<"  end_code="<<optimizer->get_failure_code()
            <<"  end_error="<<optimizer->get_end_error()
            <<"  finalCost="<<finalCost<<std::endl
            <<"  finalParams="<<finalRigidParam
            <<std::endl<<std::endl;

   return finalCost;
}



template <class PixelType>
void processingFunction( arguments args )
{
   typedef mz3::ItkImageReader<PixelType>          ReaderType;
   typedef typename ReaderType::ImageType          ImageType;
   typedef typename ImageType::Pointer             ImagePointer;

   ReaderType ird(args.fileIn);

   const int nbFrames = ird.nbFrames();
   if ( args.first == -1 )  args.first=0;
   if ( args.last  == -1 )  args.last=nbFrames-1;
   if ( args.last > (nbFrames -1)) args.last=nbFrames-1;

   const int numberOfKeptFrames = (args.last - args.first)/(args.skip+1) + 1;


   // Define images
   ImagePointer fixim = ImageType::New();
   if ( !ird.getITKFrame(fixim, args.first) )
   {
      std::cout<<"Could not get the ITK frame"<<std::endl;
      exit (EXIT_FAILURE);
   }
   const itk::ImageRegion<2> region = fixim->GetLargestPossibleRegion();

   ImagePointer movim = ImageType::New();

   ImagePointer normfixim = ImageType::New();
   normfixim->SetOrigin(fixim->GetOrigin());
   normfixim->SetSpacing(fixim->GetSpacing());
   normfixim->SetRegions(region);
   normfixim->Allocate();

   ImagePointer normmovim = ImageType::New();
   normmovim->SetOrigin(fixim->GetOrigin());
   normmovim->SetSpacing(fixim->GetSpacing());
   normmovim->SetRegions(region);
   normmovim->Allocate();


   // Define vector frames
   typedef typename ReaderType::PointSetType                PointSetType;
   typedef typename PointSetType::Pointer                   PointSetPointer;
   typedef typename ReaderType::PointsContainer             PointsContainer;
   typedef typename PointsContainer::Pointer                PointsContainerPointer;
   typedef typename ReaderType::IndexVectorType             IndexVectorType;
   typedef typename IndexVectorType::Pointer                IndexVectorPointer;

   PointsContainerPointer fiberCenters = PointsContainer::New();
   ird.getITKPointsContainer(fiberCenters);

   IndexVectorPointer fiberIndexes = IndexVectorType::New();
   ird.getITKFibersIndexContainer(fiberIndexes);


   // Compute mean of images
   std::vector<double> meanvect(nbFrames);

   std::vector<unsigned int> fiberoffsets;
   ird.getFiberOffsets( fiberoffsets );

   std::vector<PixelType> pixval( fiberoffsets.size() );

   const double nbFibersf = static_cast<double>(fiberoffsets.size());

   for ( int i=args.first; i<=args.last; i+=(1+args.skip) )
   {
      ird.getPointDataFromOffsets( pixval, fiberoffsets, i);

      meanvect[i] = std::accumulate(pixval.begin(), pixval.end(), 0.0) / nbFibersf;
   }

   // Define mask
   typedef unsigned char                MaskPixelType;
   typedef itk::Image<MaskPixelType, 2> MaskImageType;
   typedef MaskImageType::Pointer       MaskImagePointer;
   typedef itk::SpatialObject<2>        MaskType;
   typedef MaskType::Pointer            MaskPointer;

   MaskImagePointer maskim;
   MaskPointer mask;

   {
   // Get a point set from the fibers (let spacing to one)
   PointSetPointer ptset = PointSetType::New();

   ptset->SetPoints( fiberCenters );

   // compute mask
   typedef itk::ErodedMaskPointSetToImageFilter
      <PointSetType, MaskImageType> MaskFilterType;
   typename MaskFilterType::Pointer maskfilter = MaskFilterType::New();
   maskfilter->SetInput(ptset);
   maskfilter->SetSpacing(fixim->GetSpacing());
   maskfilter->SetOrigin(fixim->GetOrigin());
   maskfilter->SetSize(region.GetSize());
   maskfilter->SetErosionRadius( 2.5 );
   maskfilter->UpdateLargestPossibleRegion();

   maskim = maskfilter->GetOutput();
   maskim->DisconnectPipeline();
   }

   // Normalize first image
   maskAndNormalize<ImageType,MaskImageType>( meanvect[args.first], maskim, fixim, normfixim );



   // Define correlator
   typedef NormXCorr2<PixelType, PixelType> NormXCorr2FilterType;
   NormXCorr2FilterType correlator;
   correlator.setDenominatorThreshold( 1.0e-5 );
   correlator.setCorrelationSize( NormXCorr2FilterType::Full );


   // Define random generator for perturbing normxcorr2 input
   boost::mt19937 rng;
   //rng.seed(time(NULL));
   boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >
      gen (rng, boost::normal_distribution<>(0.0, 3.0*fixim->GetSpacing()[0]));


   mz_timer timer;
   double tESM(0.0),tFRPR(0.0),tLBFGS_vnl(0.0),tLBFGSB_vnl(0.0),tLBFGSB_vnlb(0.0);
   double cESM(0.0),cFRPR(0.0),cLBFGS_vnl(0.0),cLBFGSB_vnl(0.0),cLBFGSB_vnlb(0.0);
   double mESM(0.0),mFRPR(0.0),mLBFGS_vnl(0.0),mLBFGSB_vnl(0.0),mLBFGSB_vnlb(0.0);

   ImagePointer swapptr;
   ImagePointer normswapptr;

   for ( int i=args.first+1+args.skip;
         i<=args.last; i+=(1+args.skip) )
   {
      // Read and normalize
      ird.getITKFrame(movim, i);

      maskAndNormalize<ImageType,MaskImageType>( meanvect[i], maskim, movim, normmovim );

      //writeRaw<ImageType>(
      //   normmovim, std::string("test_") + boost::lexical_cast<std::string>(i) + ".mha");


      // Register using normxcorr2
      const double nxcc = correlator.computeCorrelation(
         normmovim, region, normfixim, region );

      // Store correlator translation
      const double tx = correlator.getWorldTranslation()[0];
      const double ty = correlator.getWorldTranslation()[1];

      //const double txc = correlator.getCenterTranslation()[0];
      //const double tyc = correlator.getCenterTranslation()[1];

      std::cout<<"NormXCorr2: "<<nxcc<<" ["<<tx<<","<<ty<<"]"<<std::endl;
      //std::cout<<"            "<<nxcc<<" ["<<txc<<","<<tyc<<"]"<<std::endl;

      //const double txp = tx;
      //const double typ = ty;
      const double txp = tx + gen();
      const double typ = ty + gen();



      // Register using different optimizers
      timer.restart();
      const double aESM =
         registerESM<ImageType,IndexVectorType,MaskType>(
            fixim, movim, fiberIndexes, mask, txp, typ, args);
      tESM += timer.elapsed_seconds();
      cESM += aESM;
      mESM = std::max( mESM,  std::fabs(aESM) );

      timer.restart();
      const double aFRPR =
         registerFRPR<ImageType,IndexVectorType,MaskType>(
            fixim, movim, fiberIndexes, mask, txp, typ, args);
      tFRPR += timer.elapsed_seconds();
      cFRPR += aFRPR;
      mFRPR = std::max( mFRPR,  std::fabs(aFRPR) );

      timer.restart();
      const double aLBFGS_vnl =
         registerLBFGS_vnl<ImageType,IndexVectorType,MaskType>(
            fixim, movim, fiberIndexes, mask, txp, typ, args);
      tLBFGS_vnl += timer.elapsed_seconds();
      cLBFGS_vnl += aLBFGS_vnl;
      mLBFGS_vnl = std::max( mLBFGS_vnl,  std::fabs(aLBFGS_vnl) );

      timer.restart();
      const double aLBFGSB_vnl =
         registerLBFGSB_vnl<ImageType,IndexVectorType,MaskType>(
            fixim, movim, fiberIndexes, mask, txp, typ, args);
      tLBFGSB_vnl += timer.elapsed_seconds();
      cLBFGSB_vnl += aLBFGSB_vnl;
      mLBFGSB_vnl = std::max( mLBFGSB_vnl,  std::fabs(aLBFGSB_vnl) );

      timer.restart();
      const double aLBFGSB_vnlb =
         registerLBFGSB_vnlb<ImageType,IndexVectorType,MaskType>(
            fixim, movim, fiberIndexes, mask, txp, typ, args);
      tLBFGSB_vnlb += timer.elapsed_seconds();
      cLBFGSB_vnlb += aLBFGSB_vnlb;
      mLBFGSB_vnlb = std::max( mLBFGSB_vnlb,  std::fabs(aLBFGSB_vnlb) );


      // swap images
      swapptr = fixim;
      fixim = movim;
      movim = swapptr;

      normswapptr = normfixim;
      normfixim = normmovim;
      normmovim = normswapptr;
   }

   ird.close();

   cESM = std::fabs(cESM)/(numberOfKeptFrames-1);
   cFRPR = std::fabs(cFRPR)/(numberOfKeptFrames-1);
   cLBFGS_vnl = std::fabs(cLBFGS_vnl)/(numberOfKeptFrames-1);
   cLBFGSB_vnl = std::fabs(cLBFGSB_vnl)/(numberOfKeptFrames-1);
   cLBFGSB_vnlb = std::fabs(cLBFGSB_vnlb)/(numberOfKeptFrames-1);

   std::cout<<std::endl;
   std::cout<<"ESM:          metric="<<cESM<<"   worst="<<mESM<<"   time="<<tESM<<" seconds"<<std::endl;
   std::cout<<"FRPR:         metric="<<cFRPR<<"   worst="<<mFRPR<<"   time="<<tFRPR<<" seconds"<<std::endl;
   std::cout<<"LBFGS_vnl:    metric="<<cLBFGS_vnl<<"   worst="<<mLBFGS_vnl<<"   time="<<tLBFGS_vnl<<" seconds"<<std::endl;
   std::cout<<"LBFGSB_vnl:   metric="<<cLBFGSB_vnl<<"   worst="<<mLBFGSB_vnl<<"   time="<<tLBFGSB_vnl<<" seconds"<<std::endl;
   std::cout<<"LBFGSB_vnlb:  metric="<<cLBFGSB_vnlb<<"   worst="<<mLBFGSB_vnlb<<"   time="<<tLBFGSB_vnlb<<" seconds"<<std::endl;
   std::cout<<std::endl;
}



int mzMain( QStringList & mzQArgs )
{
   // Get the name of the application
   const QString appName = QFileInfo( mzQArgs[0] ).fileName();

   // Define the command-line parser
   QxtCommandOptions options;

   // Set flag style to double dash to also work on windows
   options.setFlagStyle(QxtCommandOptions::DoubleDash);

   struct arguments args;
   parseOpts (mzQArgs, options, args);

   // Useful streams
   QTextStream qcout( stdout, QIODevice::WriteOnly );
   QTextStream qcerr( stderr, QIODevice::WriteOnly );

   qcout<<"Starting processing with the following arguments:"<<endl;
   qcout<<args<<endl<<endl;

   mz3::MZImagePixelType pixtype = mz3::MZ_UNKNOWN;
   mz3::ImageReader ird(args.fileIn);

   if (!ird.pixType(pixtype) and pixtype == mz3::MZ_UNKNOWN)
   {
      qcerr<<endl;
      qcerr<<"Movie"<< args.fileIn << "is not correct"<<endl;
      exit(EXIT_FAILURE);
   }
   ird.close();

   switch(pixtype) {
   case mz3::MZ_SHORT:
   case mz3::MZ_USHORT:
      processingFunction<float> ( args );
      break;
   default:
      qcerr << "Unknown pixel type " << pixtype << endl;
      exit (EXIT_FAILURE);
      break;
   }
   exit(EXIT_SUCCESS);
}

} // namespace

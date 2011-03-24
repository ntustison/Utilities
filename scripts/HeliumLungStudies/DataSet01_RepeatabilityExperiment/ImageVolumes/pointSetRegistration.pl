#!/usr/bin/perl

##
# Register point-sets
##

$utilityDirectory = "/Users/nick/pkg/Utilities/bin/";

##
# Perform refinement with point-set registration
##

$PREDIRECTORY = $ARGV[0];
$POSTDIRECTORY = $ARGV[1];

$movingPointSet = $POSTDIRECTORY . "/points.vtk";
$fixedPointSet = $PREDIRECTORY . "/points.vtk";

$domainImage = $POSTDIRECTORY . "/segmentation.nii.gz";
  
$outputDirectory = $POSTDIRECTORY . "/";

$executable = "/mnt/data1/tustison/Utilities/bin/RegisterPointSetsFFD";

$prefix1 = "psr_";
@args = ( $executable, 3,
                     "-f", $fixedPointSet, 
    																	"-m", $movingPointSet,
    																	"-a", 2.0,                                # alpha [1.0, 2.0]      
    																	"-c", 30,                                 # point-set sigma
    																	"-d", 10,                                 # kernel sigma
    																	"-r", 0.95,                               # annealing rate
    																	"-v", 1,                                  # use anisotropic covariances  
    																	"-t", 0,                                  # employ term 2
    																	"-K", 10,                                 # covariance k-neighborhood 
    																	"-k", 100,                                # evaulation k-neighborhood    
    																	"-i", "30",                               # number of iterations at each level
    																	"-B", 3,                                  # spline order
    																	"-R", "2x2x1",
    																	"-l", 3,  
    																	"-L", "3x2",                              # point labels to use
    																	"-W", "1.0x1.0",                          # gradient weights (for dmffd approach)
    																	"-P", "1.0x1.0",                          # number of points used (in percentages)
     																	"-o", $outputDirectory . $prefix1,
    																	"-p", 0,
    																	"-I", $domainImage, 
    																	"-h", 1,                                  # use input as samples
    																	"-G", 0,                                  # use Rueckert gradient
                     "-s", 5000,
                     "-j", 5000,
                     "-x", 0
          																	);
system( @args ) == 0 or die "system @args failed: $?";

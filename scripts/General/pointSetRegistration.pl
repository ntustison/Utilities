#!/usr/bin/perl

#$algorithm = "/Users/nick/pkg/Utilities/bin/3D/float/RegisterPointSetsFFD";
$algorithm = "/Users/nick/pkg/Utilities/bin/2D/float/RegisterPointSetsFFD";

$fixedPointSet = "/Users/nick/pkg/Projects/PointSetMetric/Simple/left.nii.gz";
$movingPointSet = "/Users/nick/pkg/Projects/PointSetMetric/Simple/right.nii.gz";

$outputDirectory = "/Users/nick/pkg/Projects/PointSetMetric/Simple/";
$naming = "registered";

@args = ( $algorithm, "-f", $fixedPointSet, 
																						"-m", $movingPointSet,
                      "-a", 2.0,                                # alpha [1.0, 2.0]      
																						"-c", 10,                                  # point-set (regularization) sigma
                      "-r", 1.0,                                # annealing rate
                      "-h", 1,                                  # use input as samples
                      "-s", 1000,                               # number of fixed samples ( if -h = 0 )
                      "-j", 1000,                               # number of moving samples ( if -h = 0 )
                      "-k", 50,                                 # evaluation k-neighborhood 
																						"-v", 0,                                  # use anisotropic covariances
                      "-d", 100,                                # kernel sigma ( if -v = 1 )
                      "-K", 4,                                  # covariance k-neighborhood  ( if -v = 1 ) 
                      "-t", 1,                                  # use term 2 (regularization term)
																						"-n", 3, 
																						"-i", "10x5x2", 
																						"-B", 3,
																						"-R", "3x3",
                      "-l", 20,  
																						"-o", $outputDirectory . $naming,
                      "-p", 0,                                  # prolificacy
                      "-D", "1x1", 
                      "-I", $fixedPointSet,
#                      "-x", 2,                                  # expansion factor
#                      "-P", 1.0x1.0x1.0",
#                      "-Z", 100x100x100",
#                      "-O", 0.0x0.0x0.0"                       
																						);
system( @args ) == 0 or die "system @args failed: $?"

#!/usr/bin64/perl

##
# Register point-sets
##

$utilityDirectory = "/home/tustison/Utilities/bin64/";

##
# Perform refinement with point-set registration
##

$PREDIRECTORY = $ARGV[0];
$POSTDIRECTORY = $ARGV[1];
$PRESERIES = $ARGV[2];
$POSTSERIES = $ARGV[3];

print $PREDIRECTORY . "\n";
print $POSTDIRECTORY . "\n";
print $PRESERIES . "\n";
print $POSTSERIES . "\n";


$movingPointSet = $POSTDIRECTORY . "/points.vtk";
$fixedPointSet = $PREDIRECTORY . "/points.vtk";

$domainImage = $POSTDIRECTORY . "/mask_leftandright_run01.nii.gz";
  
$outputDirectory = $POSTDIRECTORY . "/";

$executable = "/home/tustison/Utilities/bin64/RegisterPointSetsFFD";

$prefix1 = "psr_";
@args = ( $executable, 3, "--point-sets", "[".$fixedPointSet.",".$movingPointSet."]", 
  "--transformation", "[3,2x2x1,".$domainImage."]", "--similarity", "[2.0,0.97,30,0,50,1,8,10.0]",
  "--optimization", "[50x25x10x10,0.5,4,3]", "--output", $outputDirectory . $prefix1, "--verbose");
system( @args ) == 0 or die "system @args failed: $?";


$executable = "/home/tustison/ANTS/bin64/WarpImageMultiTransform";
@args = ( $executable, 3, $PREDIRECTORY . "/" . $PRESERIES . ".nii.gz", 
  $POSTDIRECTORY . "/psr_Warped.nii.gz", "-R", $POSTDIRECTORY . "/" . $POSTSERIES . ".nii.gz",
  $POSTDIRECTORY  . "/" . "psr_Warp.nii.gz" ); 
system( @args ) == 0 or die "system @args failed: $?";

$executable = "/home/tustison/Utilities/bin64/ConvertScalarImageToRGB";
@args = ( $executable, 3, 
  $POSTDIRECTORY . "/psr_Warped.nii.gz", $POSTDIRECTORY . "/psr_Warped_hot.mha",
  "hot" ); 
system( @args ) == 0 or die "system @args failed: $?";

$executable = "/home/tustison/Utilities/bin64/ConvertScalarImageToRGB";
@args = ( $executable, 3, $POSTDIRECTORY . "/" . $POSTSERIES . ".nii.gz", 
  $POSTDIRECTORY . "/" . $POSTSERIES . "_hot.mha",
  "hot" ); 
system( @args ) == 0 or die "system @args failed: $?";

$executable = "/home/tustison/Utilities/bin64/ConvertScalarImageToRGB";
@args = ( $executable, 3, $PREDIRECTORY . "/" . $PRESERIES . ".nii.gz", 
  $PREDIRECTORY . "/" . $PRESERIES . "_hot.mha",
  "hot" ); 
system( @args ) == 0 or die "system @args failed: $?";


#@args = ( $executable, 3,
#                     "-f", $fixedPointSet, 
#    																	"-m", $movingPointSet,
#    																	"-a", 2.0,                                # alpha [1.0, 2.0]      
#    																	"-c", 30,                                 # point-set sigma
#    																	"-d", 10,                                 # kernel sigma
#    																	"-r", 0.95,                               # annealing rate
#    																	"-v", 1,                                  # use anisotropic covariances  
#    																	"-t", 0,                                  # employ term 2
#    																	"-K", 10,                                 # covariance k-neighborhood 
#    																	"-k", 100,                                # evaulation k-neighborhood    
#    																	"-i", "30",                               # number of iterations at each level
#    																	"-B", 3,                                  # spline order
#    																	"-R", "2x2x1",
#    																	"-l", 3,  
#    																	"-L", "3x2",                              # point labels to use
#    																	"-W", "1.0x1.0",                          # gradient weights (for dmffd approach)
#    																	"-P", "1.0x1.0",                          # number of points used (in percentages)
#     																	"-o", $outputDirectory . $prefix1,
#    																	"-p", 0,
#    																	"-I", $domainImage, 
#    																	"-h", 1,                                  # use input as samples
#    																	"-G", 0,                                  # use Rueckert gradient
#                     "-s", 5000,
#                     "-j", 5000,
#                     "-x", 0
#          																	);
#system( @args ) == 0 or die "system @args failed: $?";

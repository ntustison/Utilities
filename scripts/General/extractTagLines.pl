#!/usr/bin/perl

$algorithm = "/Users/nick/pkg/Utilities/bin/3D/float/ExtractTagLinePoints";
#$algorithm = "/Users/nick/pkg/Utilities/bin/2D/float/RegisterImagesFFD";

# Usage: /Users/nick/pkg/Utilities/bin/3D/float/ExtractTagLinePoints 
#   imageFile maskFile outputFilePrefix tagSpacing [thresholdPercentage] 
#   [numberOfAngleSteps] [numberOfTagSpacingSteps] [angleOffset] [tagSpacingFactor]


$imageFile = 
$maskFile = 

$outputDirectory = 
$naming = "registered";

@args = ( $algorithm, $imageFile, 
                      $maskFile,
                      $outputDirectory . $naming,
                      28,
                      95
                      5
                      5
                      10
                      0.1 
																						);
system( @args ) == 0 or die "system @args failed: $?"


#  Cost/similarity Options: 
#     0: optical flow
#     1: n.a.
#     2: Mutual Information (not good)
#     3: MutualInformation (good - 2nd most recommended) 
#     4: Cross Correlation of radius 5 (recommended) 
#     5: Cross Correlation of radius 2 (similar to 4) 


#!/usr/bin/perl

use strict;

my $fixedImage = $ARGV[0];
my $movingImage = $ARGV[1];
my $outputDirectory = $ARGV[2];

my $antsDirectory = "/mnt/data1/tustison/PICSL/ANTS/bin/";
my $utilitiesDirectory = "/mnt/data1/tustison/Utilities/bin/";

my $ants = $antsDirectory . "ANTS";
my $antsScript = "/Users/nick/pkg/Utilities/scripts/antsMediaTest.sh";

my $dimension = 3;

## 
# Create otsu threshold fixed and moving images 
##

my $otsu = $antsDirectory . "ThresholdImage";
my $fixedOtsuImage = $outputDirectory . "fixedImageOtsu.nii.gz";
#`$otsu $dimension $fixedImage $fixedOtsuImage Otsu 4`;

my $otsu = $antsDirectory . "ThresholdImage";
my $movingOtsuImage = $outputDirectory . "movingImageOtsu.nii.gz";
#`$otsu $dimension $movingImage $otsuImage Otsu 4`;



my $count = 0;
my @gradientSteps;
for ( my $g = 0.4; $g <= 0.5; $g+=0.1 )
  {
  $gradientSteps[$count++] = $g;
  }

##
# Transformation/Regularizations
##

my @transformations = ( "SyN", "SyN", "Exp", "Elast", "Elast", "Elast", "Elast" );
my @regularizations = ( "Gauss", "Gauss", "Gauss", "Gauss", "DMFFD", "DMFFD", "DMFFD" );

my $gradientFieldSigma = 3.0;
my $totalFieldSigma = 1.0;

# greedy SyN
# time-varying SyN
my $meaninglessParam = 2;
my $integrationDelta = 0.01; 
# exponential
my $expTimeSteps = 8;
# deformable
# dmffd(order=2)
# dmffd(order=3)
# dmffd(order=4)

my @metrics = ( "MSQ", "PR", "PR", "PR", "PR" );
my @radii = ( 0, 2, 4, 6, 8 );
# mean-squared difference
# cross correlation (radius=1) 
# cross correlation (radius=2)
# cross correlation (radius=3)

my $iterations = "1x0x0";


##
# Run ANTS for all the different scenarios
##

my $numberOfRegistrations = 0;

for ( my $i = 0; $i < @gradientSteps; $i++ )
  {
  my $gradientStep = $gradientSteps[$i]; 
  for( my $j = 0; $j < @transformations; $j++ )
    {
    my $transformation = $transformations[$j] . "[" . $gradientStep;
    my $regularization = $regularizations[$j] . "[" . $gradientFieldSigma . "," . $totalFieldSigma;
    if( $j == 0 )      # greedy SyN
      {
      $transformation = $transformation . "]";
      $regularization .= "]";
      }
    elsif( $j == 1 )   # time-varying SyN 
      {
      $transformation = $transformation . "," . $meaninglessParam . "," . $integrationDelta . "]";      
      $regularization .= "]";
      }
    elsif( $j == 2 )   # exponential 
      {
      $transformation = $transformation . "]";      
      $regularization .= "]";
      }
    elsif( $j > 2 )   # elastic,dmffd 
      {
      $transformation = $transformation . "]";      
      if( $j == 3 )
        {
        $regularization .= "]";
        }
      else
        {
        my $order = $j - 2; 
        $regularization .= "," . $order . "]";
        }  
      }
    for( my $k = 0; $k < @metrics; $k++ )
      {
      my $metric = $metrics[$k] . "[" . $fixedImage . "," . $movingImage . ",1," . $radii[$k] . "]";  
       
      my $output = $outputDirectory . "ants_" . $i . "_" . $j . "_" . $k;
      my $deformedImage = $output . "deformed.nii";
       
      # perform the pair-wise registration 
      $numberOfRegistrations++;
      my @args = ( "sh", $antsScript, $dimension, 
                $output, $iterations, $transformation, $regularization, $metric
                  ); 
      system( @args ) == 0 or die "system @args failed: $?";
#      print $numberOfRegistrations . ": " . "@args\n";
      
      exit;

#       now perform registration assessment
#      # 
#         1. perform similarity metric assessment
#      $measureImageSimilarity = $antsDirectory . "MeasureImageSimilarity";
#      for( $sim = 0; $sim < 3; $sim++ )
#        {
#        $similarityOutput = $output . "similarityMetric" . $sim . ".txt";
#        print "    $measureImageSimilarity 3 $sim $fixedImage $deformedImage > $similarityOutput \n";
#        } 
#         2. create deformed otsu image
#      $warpImage = $antsDirectory . "WarpImageMultiTransform";
#      $deformedOtsuImage = $output . "deformedOtsu.nii.gz";
#      $vectorField = $output . "Warp.nii";
#      $affine = $output . "Affine.txt";
#      `$warpImage $dimension $movingOtsuImage $deformedOtsuImage $vectorField $affine -R $fixedImage --use-NN`;
#
#         3. calculate dice statistics
#      $imageMath = $antsDirectory . "ImageMath";
#      $outputDiceFile = $output . "diceStats.txt";
#      $outputMinDistImage = $output . "minDistSum.nii.gz";
#      `$imageMath $dimension $outputDiceFile DiceAndMinDistSum $fixedOtsuImage $movingOtsuImage $outputMinDistImage`; 
#            
#         4. calculate jacobian over labels
#      $compose = $antsDirectory . "ComposeMultiTransform";
#      $composedVectorField = $output . "TotalWarp.nii.gz";
#      $jacobian = $antsDirectory . "CreateJacobianDeterminantImage";
#      $jacobianImage = $output . "Jacobian.nii.gz";
#      $movingLabelsFile = $output . "movingLabelsStats.txt";
#      $jacobianLabelsFile = $output . "jacobianLabelsStats.txt";
#      
#      `$compose $dimension $composedVectorField -R $fixedImage $vectorField $affine`;
#      `$jacobian $dimension $composedVectorField $jacobianImage 0`;
#      `$imageMath $dimension $movingLabelsFile LabelStats $movingOtsuImage $movingOtsuImage`;
#      `$imageMath $dimension $jacobianLabelsFile LabelStats $deformedOtsuImage $jacobianImage`;              
      
      }
    }  
  }    



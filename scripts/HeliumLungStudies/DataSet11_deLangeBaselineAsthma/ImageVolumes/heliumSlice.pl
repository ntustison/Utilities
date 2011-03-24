#!/usr/bin/perl

$utilityDirectory = "/home/tustison/Utilities/bin64/";
$antsDirectory = "/home/tustison/ANTS/bin64/";

$prefix = $ARGV[0];
$inputImage = "${prefix}.nii.gz";
$anisotropicImage = "${prefix}_aniso.nii.gz";
$maskImage = "${prefix}_mask.nii.gz";
$closeMaskImage = "${prefix}_closing_mask.nii.gz";
$innerMaskImage = "${prefix}_inner_mask.nii.gz";
$distanceImage = "${prefix}_distance.nii.gz";
$mediastinumMaskImage = "${prefix}_mediastinum.nii.gz";
$featureImage = "${prefix}_hessian.nii.gz";
$featureImageN4 = "${prefix}_hessian_n4.nii.gz";
$vesselMaskImage = "${prefix}_vesselmask.nii.gz";
$otsuImage = "${prefix}_otsu.nii.gz";
$vesselImage = "${prefix}_vessels.nii.gz";
$priorLabelImage = "${prefix}_priorlabel.nii.gz";

# `${utilityDirectory}/GradientAnisotropicDiffusionImageFilter 2 ${inputImage} ${anisotropicImage} 5 0.05 2`;
`${utilityDirectory}/BinaryMorphology 2 ${maskImage} ${closeMaskImage} 2 100 1`;
`${utilityDirectory}/GenerateDistanceImage 2 ${closeMaskImage} ${distanceImage}`;
`${utilityDirectory}/ThresholdImage 2 ${distanceImage} ${innerMaskImage} -100000000000 -15 1 0`;
`${utilityDirectory}/BinaryOperateImages 2 ${innerMaskImage} x ${maskImage} ${innerMaskImage}`;
`${antsDirectory}/SmoothImage 2 ${innerMaskImage} 2.0 ${prefix}_mask_center_smooth.nii.gz`;

`${utilityDirectory}/HessianBasedFeatures 2 ${inputImage} ${featureImage} 1 0.01 3.0 20 5.0 0.5 0.5 0`;
`${utilityDirectory}/BinaryOperateImages 2 ${prefix}_mask_center_smooth.nii.gz x ${featureImage} ${featureImage}`;
# `${utilityDirectory}/RescaleImageIntensity 2 ${featureImage} ${featureImage} 100 1000`;
#
# @args = ( "${antsDirectory}/N4BiasFieldCorrection", 2,
#   '-i', $featureImage,
#   '-s', 1,
#   '-c', '[50x50x50x50,0.00001]',
#   '-b', '[100,3]',
#   '-t', '[0.15,0.01]',
#   '-o', $featureImageN4
#  );
# system( @args ) == 0 or die "system @args failed: $?";

`${utilityDirectory}/OtsuThresholdImage 2 ${featureImage} ${otsuImage} 1 200 ${innerMaskImage} 1`;
`${utilityDirectory}/ThresholdImage 2 ${otsuImage} ${vesselImage} 2 2 0 1`;
`${utilityDirectory}/BinaryOperateImages 2 ${vesselImage} x ${maskImage} ${vesselMaskImage}`;


$numberOfClasses = 7;
$threshold = 3;

@volume = `${utilityDirectory}/CalculateVolumeFromBinaryImage 2 $vesselMaskImage`;
if( $volume[1] > 0 )
  {
  # First get prior classification on whole image
  @args = ( "${utilityDirectory}/AtroposSegmentation", 2,
    "--initialization", "kmeans[${numberOfClasses}]",
    "--intensity-image", "${inputImage}",
    "--convergence", '[0]',
    "--output", "${priorLabelImage}" );
  system( @args ) == 0 or die "system @args failed: $?";

  `${utilityDirectory}/BinaryOperateImages 2 ${priorLabelImage} x ${vesselMaskImage} ${priorLabelImage}`;

  @args = ( "${utilityDirectory}/AtroposSegmentation", 2,
    "--initialization", "PriorLabelImage[${numberOfClasses},${priorLabelImage},0.0]",
    "--intensity-image", "${inputImage}",
    "--mask-image", "${vesselMaskImage}",
    "--convergence", '[5,0.00001]',
    "--likelihood-model", 'Gaussian',
    "--winsorize-outliers", 'BoxPlot[0.25,0.75,1.5]',
    "--mrf", "[5,2x2]",
    "--output", "${prefix}_atropos.nii.gz" );
  system( @args ) == 0 or die "system @args failed: $?";
  }
else
  {
  `/bin/cp ${vesselMaskImage} ${prefix}_atropos.nii.gz`;
  }

##
#   1. Get connected components and remove clusters which are < 40 total voxels
#   2. Remove clusters
#   3. Perform morphological opening and then get rid of clusters which are
#      < maxClusterSize and have elongation > minElongation
##

$maxClusterSize = 25;
$minElongation = 4;


`${utilityDirectory}/ThresholdImage 2 ${prefix}_atropos.nii.gz ${prefix}_atropos_2class.nii.gz 1 ${threshold} 1 0`;
`${utilityDirectory}/GetConnectedComponents 2 ${prefix}_atropos_2class.nii.gz ${prefix}_atropos_2class_cc.nii.gz 0`;

open( LGM, "${utilityDirectory}/LabelGeometryMeasures 2 ${prefix}_atropos_2class_cc.nii.gz |" );
$count = 1;
while( <LGM> )
  {
  if( $count > 1 )
    {
    ( $label, $volume, $eccentricity, $elongation, $orientation, $centroid, $axislength, $bbox ) = split;
    if( $volume < $maxClusterSize )
      {
      `${utilityDirectory}/UnaryOperateImage 2 ${prefix}_atropos_2class_cc.nii.gz r 0 ${prefix}_atropos_2class_cc.nii.gz ${label} 0`;
      }
    }
  $count++;
  }
close( LGM );

 `${utilityDirectory}/ThresholdImage 2 ${prefix}_atropos_2class.nii.gz ${prefix}_atropos_2class.nii.gz 0 0 0 1`;
 `${utilityDirectory}/BinaryMorphology 2 ${prefix}_atropos_2class.nii.gz ${prefix}_atropos_2class.nii.gz 3 1 1 1`;
 `${utilityDirectory}/GetConnectedComponents 2 ${prefix}_atropos_2class.nii.gz ${prefix}_atropos_2class_cc.nii.gz 0`;
 `${utilityDirectory}/GetConnectedComponents 2 ${prefix}_atropos_2class_cc.nii.gz ${prefix}_atropos_2class_cc.nii.gz 0 1 ${mediastinumMaskImage}`;

open( LGM, "${utilityDirectory}/LabelGeometryMeasures 2 ${prefix}_atropos_2class_cc.nii.gz |" );
$count = 1;
while( <LGM> )
  {
  if( $count > 1 )
    {
    ( $label, $volume, $eccentricity, $elongation, $orientation, $centroid, $axislength, $bbox ) = split;
    if( $elongation > $minElongation || $volume < $maxClusterSize )
      {
      print $label . " -> " . $volume . "\n";
      `${utilityDirectory}/UnaryOperateImage 2 ${prefix}_atropos_2class_cc.nii.gz r 0 ${prefix}_atropos_2class_cc.nii.gz ${label} 0`;
      }
    }
  $count++;
  }
close( LGM );
`${utilityDirectory}/GetConnectedComponents 2 ${prefix}_atropos_2class_cc.nii.gz ${prefix}_atropos_2class_cc.nii.gz 1`;

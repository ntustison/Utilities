#!/usr/bin/perl

$baseDirectory = "/Users/nick/pkg/Projects/HeliumLungStudies/NormalImages/PCALevelSets/PCA_2/";
$imageDirectory = "${baseDirectory}/AxialImages/";
$antsDirectory = "/Users/nick/pkg/PICSL/ANTS/bin/";

$template = "${baseDirectory}/OUT3template.nii.gz";

# opendir( DIR, $imageDirectory );
# while( defined( $file = readdir( DIR ) ) )
#   {
# 		if( $file =~ m/\.nii\.gz$/ )
# 				{
#     print $file . "\n";
#
#     $filePrefix = $file;
# 				$outputFile = "${imageDirectory}/warped_${file}";
#     $outputFilePrefix = $outputFile;
#     $outputFilePrefix =~ s/\..*//;
# 				$segmentedFile = "${imageDirectory}/warped_segmented_${file}";
# 				$distanceFile = "${imageDirectory}/warped_distance_${file}";
#
#     `${antsDirectory}/ANTS 3 -m PR[${template},${imageDirectory}/${file},1,4] -t SyN[0.5] -r Gauss[3,2] -i 0 -o $outputFile`;
#     `${antsDirectory}/WarpImageMultiTransform 3 ${imageDirectory}/${file} $outputFile -R ${template} ${outputFilePrefix}Affine.txt`;
#
#     `/Users/nick/pkg/Utilities/bin/SegmentHeliumLungs 3 $outputFile $segmentedFile`;
#     `/Users/nick/pkg/Utilities/bin/ThresholdImage 3 $segmentedFile $segmentedFile 2 3 1 0`;
#     `/Users/nick/pkg/Utilities/bin/GenerateDistanceImage 3 $segmentedFile $distanceFile`;
# 				}
#   }

for( $count = 0; $count < 157; $count++ )
  {
		if( $count < 10 )
				{
				$index = "00${count}";
				}
		elsif( $count < 100 )
				{
				$index = "0${count}";
				}
		else
				{
				$index = "${count}";
				}
  $file = "${imageDirectory}/warped_segmented_he3_${index}.nii.gz";
		$distanceFile = "${imageDirectory}/warped_distance_he3_${index}.nii.gz";
  print $file . "\n";
  `/Users/nick/pkg/Utilities/bin/GenerateDistanceImage 3 $file $distanceFile`;
  }

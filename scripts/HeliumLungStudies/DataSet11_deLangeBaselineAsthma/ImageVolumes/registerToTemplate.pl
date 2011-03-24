#!/user/bin/perl -w

$outputDirectory = $ARGV[0];
$uncorrectedHeliumImage = $ARGV[1];

@names = ( 413, 414, 416, 417, 418, 419, 420, 422, 423, 424, 429, 432, 438, 439,
441, 442, 443, 444, 446, 447, 448, 449, 451, 452, 453, 454, 456, 458, 459, 460,
461, 462, 466, 467, 468, 471, 475, 480, 484, 485, 487, 494, 496, 497, 498, 504,
505, 515, 531, 533, 538, 549, 550, 553, 568, 572, 620, 621, 635, 637, 638, 639,
645, 646, 648, 649, 660, 669, 670, 672, 674, 677, 678, 684, 699, 703 );

@diagnosis = ( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0,
 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0 );

$utilityDirectory = "/home/tustison/Utilities/bin64/";
$antsDirectory = "/home/tustison/ANTS/bin64/";
$shapeIntensityTemplateDirectory = "/home/tustison/Data/HeliumLungStudies/Templates/ANTSShapeIntensityTemplate/Axial/UnbiasedTemplate/";
$shapeIntensityTemplate = $shapeIntensityTemplateDirectory . "OUT3template_biascorrected.nii.gz";
$shapeIntensityTemplateSegmentation = $shapeIntensityTemplateDirectory . "OUT3template_segmentation2.nii.gz";



$movingImage = "${outputDirectory}/initial_segmentation_run00.nii.gz";

$output = "${outputDirectory}/antsToTemplate";

@args = ( "${antsDirectory}/ANTS", 3,
        "-m", "MSQ[${shapeIntensityTemplateSegmentation},${movingImage},1,4]",
        "-t", "SyN[0.5]",
        "-i", "100x100x100",
        "-r", "Gauss[3,1]",
        "-o", "${output}.nii.gz" );
system( @args ) == 0 or die "system @args failed: $?";

# $isAsthma = 0;
# $index = 0;
# for( $i = 0; $i < @diagnosis; $i++ )
#   {
#   if( $outputDirectory =~ /${names[$i]}/ )
#     {
#     $index = $i;
#     }
#
#   if( $outputDirectory =~ /${names[$i]}/ && $diagnosis[$i] == 1 )
#     {
#     $isAsthma = 1;
#     }
#   }
#
# $whichDirectory = '/home/tustison/Data/HeliumLungStudies/DataSet11_deLangeBaselineAsthma/ImageVolumes/WarpedNormalImages/';
# if( $isAsthma )
#   {
#   $whichDirectory = '/home/tustison/Data/HeliumLungStudies/DataSet11_deLangeBaselineAsthma/ImageVolumes/WarpedAsthmaImages/';
#   }
#
# $movingImage = "${outputDirectory}/apocrita_5class_both_run02.nii.gz";
# $warpedImage = "${whichDirectory}/warpedHe3#${names[$index]}.nii.gz";
# @args = ( "${antsDirectory}/WarpImageMultiTransform", 3,
#         $movingImage, $warpedImage, '-R', $shapeIntensityTemplate,
#         "${output}Warp.nii.gz", "${output}Affine.txt" );
# system( @args ) == 0 or die "system @args failed: $?";
#
# return;






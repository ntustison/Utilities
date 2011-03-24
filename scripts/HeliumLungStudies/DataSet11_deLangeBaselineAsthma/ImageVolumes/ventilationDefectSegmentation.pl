#!/user/bin/perl -w

use POSIX;

$outputDirectory = $ARGV[0];
$originalHeliumImage = "${outputDirectory}/${ARGV[1]}";
$inputImage = "${outputDirectory}/he3_bias_corrected_run00.nii.gz";
$initialMaskImage = "${outputDirectory}/initial_segmentation_run00.nii.gz";
$maskImage = "${outputDirectory}/refined_segmentation_run00.nii.gz";

$baseDirectory = "/home/tustison/Data/HeliumLungStudies/DataSet11_deLangeBaselineAsthma/ImageVolumes/";
$utilityDirectory = "/home/tustison/Utilities/bin64/";
$antsDirectory = "/home/tustison/ANTS/bin64/";

sleep( int( rand(5) ) );


# `/bin/cp ${initialMaskImage} ${maskImage}`;

#
# refine the initial segmentation mask
#
`${utilityDirectory}/UnaryOperateImage 3 ${initialMaskImage} + 1 ${maskImage}`;

@args = ( "${utilityDirectory}AtroposSegmentation", 3,
  '-i', "PriorLabelImage[2,${maskImage},0.9]",
  '-a', "${inputImage}",
  '-c', '[4,0.00001]',
  '-k', 'Gaussian',
  '-m', '[5,2x2x1]',
  '-w', 'BoxPlot[0.25,0.75,1.5]',
  '-o', "${maskImage}",
  '-l', "1[0.25,0.75]",
  '-l', "2[0.25,0.75]" );
system( @args ) == 0 or die "system @args failed: $?";

`${utilityDirectory}/ThresholdImage 3 ${maskImage} ${maskImage} 2 2 1 0`;

$shapeIntensityTemplateDirectory = "/home/tustison/Data/HeliumLungStudies/Templates/ANTSShapeIntensityTemplate/Axial/UnbiasedTemplate/";
# $shapeIntensityTemplate = $shapeIntensityTemplateDirectory . "OUT3template_segmentation_singleLabel.nii.gz";
$mediastinumTemplate = $shapeIntensityTemplateDirectory . "OUT3template_lung_mediastinum_segmentation.nii.gz";

#
# Warp the mediastinum mask from the template to the individual subject
#
$antsOutput = "${outputDirectory}/antsToTemplate";
$mediastinumMask = "${outputDirectory}/mediastinumMask.nii.gz";
@args = ( "${antsDirectory}/WarpImageMultiTransform", 3,
        $mediastinumTemplate, $mediastinumMask,
        "-R", $inputImage,
        "--use-NN",
        "-i", "${antsOutput}Affine.txt",
        "${antsOutput}InverseWarp.nii.gz" );
system( @args ) == 0 or die "system @args failed: $?";

#
# Process each subject slice by slice
#
@imageInfo = `${utilityDirectory}/GenerateStatisticsFromImage $inputImage`;
$begin = index( $imageInfo[1], '[' );
$end = index( $imageInfo[1], ']' );
$sizeString = substr( $imageInfo[1], $begin+1, $end - $begin - 1 );
@sizes = split( /, /, $sizeString );

$smallestSize = $sizes[0];
$smallestIndex = 0;
for( $i = 1; $i < @sizes; $i++ )
  {
  if( $sizes[$i] < $smallestSize )
    {
    $smallestSize = $sizes[$i];
    $smallestIndex = $i;
    }
  }
$smallestSize -= 1;

$defectsSliceFile = "${outputDirectory}/defectsSlice.txt";
open( DEF, ">${defectsSliceFile}" );

for( $i = 0; $i <= $smallestSize; $i++ )
  {
  $sliceImage = "${outputDirectory}/slice_${i}.nii.gz";
  $originalSliceImage = "${outputDirectory}/originalslice_${i}.nii.gz";
  $sliceMaskImage = "${outputDirectory}/slice_${i}_mask.nii.gz";
  $sliceMediastinumMaskImage = "${outputDirectory}/slice_${i}_mediastinum.nii.gz";
  `${utilityDirectory}/ExtractSliceFromImage $maskImage $sliceMaskImage $smallestIndex $i`;

  # erode the 2D mask slightly
  `${utilityDirectory}/GenerateDistanceImage 2 $sliceMaskImage $sliceMaskImage`;
  `${utilityDirectory}/ThresholdImage 2 $sliceMaskImage $sliceMaskImage -100000000 -2 1 0`;

  @volume = `${utilityDirectory}/CalculateVolumeFromBinaryImage 2 $sliceMaskImage`;
  if( $volume[1] > 0 )
    {

    `${utilityDirectory}/ExtractSliceFromImage $inputImage $sliceImage $smallestIndex $i`;
    `${utilityDirectory}/ExtractSliceFromImage $mediastinumMask $sliceMediastinumMaskImage $smallestIndex $i`;
    `${utilityDirectory}/ExtractSliceFromImage $originalHeliumImage $originalSliceImage $smallestIndex $i`;
    `${utilityDirectory}/ThresholdImage 2 $sliceMediastinumMaskImage $sliceMediastinumMaskImage 4 5 1 0`;
    `/usr/bin/perl ${baseDirectory}/heliumSlice.pl ${outputDirectory}/slice_${i}`;

    if( ! -e "${outputDirectory}/slice_${i}_otsu.nii.gz" )
      {
      `/bin/cp $sliceMaskImage ${outputDirectory}/slice_${i}_otsu.nii.gz`;
      }
    if( ! -e "${outputDirectory}/slice_${i}_atropos_2class_cc.nii.gz" )
      {
      `/bin/cp $sliceMaskImage ${outputDirectory}/slice_atropos_${i}.nii.gz`;
      }
    else
      {
      `/bin/mv ${outputDirectory}/slice_${i}_atropos_2class_cc.nii.gz ${outputDirectory}/slice_atropos_${i}.nii.gz`;
      }
    if( ! -e "${outputDirectory}/slice_${i}_atropos.nii.gz" )
      {
      `/bin/cp $sliceMaskImage ${outputDirectory}/slice_${i}_atropos.nii.gz`;
      }
    }
  else
    {
    `/bin/cp $sliceMaskImage ${outputDirectory}/slice_${i}_atropos.nii.gz`;
    `/bin/cp $sliceMaskImage ${outputDirectory}/slice_atropos_${i}.nii.gz`;
    `/bin/cp $sliceMaskImage ${outputDirectory}/slice_${i}_hessian.nii.gz`;
    `/bin/cp $sliceMaskImage ${outputDirectory}/slice_${i}_hessian_n4.nii.gz`;
    `/bin/cp $sliceMaskImage ${outputDirectory}/slice_${i}_otsu.nii.gz`;
    `/bin/cp $sliceMaskImage ${outputDirectory}/slice_${i}_inner_mask.nii.gz`;
   }
  @statsString = `${utilityDirectory}/CalculateFirstOrderStatisticsFromImage 2 ${outputDirectory}/slice_atropos_${i}.nii.gz`;
  $statsString[1] =~ s/\r|\n//g;
  @stats = split( / /, $statsString[1] );

  print DEF "${i} ${stats[11]}\n"
  }

close( DEF );

$atroposImage = "${outputDirectory}/atropos.nii.gz";
`${utilityDirectory}/ConvertImageSeries $outputDirectory slice_%d_atropos.nii.gz $atroposImage 0 $smallestSize 1`;
`${utilityDirectory}/ChangeImageInformation 3 ${atroposImage} ${atroposImage} 4 ${inputImage}`;

$defectImage = "${outputDirectory}/defects.nii.gz";
`${utilityDirectory}/ConvertImageSeries $outputDirectory slice_atropos_%d.nii.gz $defectImage 0 $smallestSize 1`;
`${utilityDirectory}/ChangeImageInformation 3 ${defectImage} ${defectImage} 4 ${inputImage}`;

$hessianImage = "${outputDirectory}/hessian.nii.gz";
`${utilityDirectory}/ConvertImageSeries $outputDirectory slice_%d_hessian.nii.gz $hessianImage 0 $smallestSize 1`;
`${utilityDirectory}/ChangeImageInformation 3 ${hessianImage} ${hessianImage} 4 ${inputImage}`;

$hessianImageN4 = "${outputDirectory}/hessian_n4.nii.gz";
`${utilityDirectory}/ConvertImageSeries $outputDirectory slice_%d_hessian_n4.nii.gz $hessianImageN4 0 $smallestSize 1`;
`${utilityDirectory}/ChangeImageInformation 3 ${hessianImageN4} ${hessianImageN4} 4 ${inputImage}`;

$otsuImage = "${outputDirectory}/otsu.nii.gz";
`${utilityDirectory}/ConvertImageSeries $outputDirectory slice_%d_otsu.nii.gz $otsuImage 0 $smallestSize 1`;
`${utilityDirectory}/ChangeImageInformation 3 ${otsuImage} ${otsuImage} 4 ${inputImage}`;

$innerMaskImage = "${outputDirectory}/innerMask.nii.gz";
`${utilityDirectory}/ConvertImageSeries $outputDirectory slice_%d_inner_mask.nii.gz $innerMaskImage 0 $smallestSize 1`;
`${utilityDirectory}/ChangeImageInformation 3 ${innerMaskImage} ${innerMaskImage} 4 ${inputImage}`;

$midSlice = int( 0.5 * ${smallestSize} + 0.5 );
`${utilityDirectory}/CreateColoredLabeledImage 2 ${outputDirectory}/slice_${midSlice}.nii.gz ${outputDirectory}/slice_atropos_${midSlice}.nii.gz ${outputDirectory}/latex_atropos_slice_${midSlice}.png hsv`;
`${utilityDirectory}/ConvertScalarImageToRGB 2 ${outputDirectory}/slice_${midSlice}.nii.gz ${outputDirectory}/latex_slice_${midSlice}.png nomask grey`;
`${utilityDirectory}/ConvertScalarImageToRGB 2 ${outputDirectory}/originalslice_${midSlice}.nii.gz ${outputDirectory}/latex_originalslice_${midSlice}.png nomask grey`;

unlink glob( "${outputDirectory}/slice*nii.gz" );
unlink glob( "${outputDirectory}/originalslice*nii.gz" );


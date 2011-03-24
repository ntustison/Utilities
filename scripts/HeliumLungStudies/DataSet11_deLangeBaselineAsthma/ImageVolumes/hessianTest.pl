#!/user/bin/perl -w

use File::Copy;

$outputDirectory = $ARGV[0];
$uncorrectedHeliumImage = $ARGV[1];
# $tmpDirectory = $ENV{ 'TMPDIR' };
# $tmpDirectory .= "/";
$tmpDirectory = "${outputDirectory}/";

$utilityDirectory = "/home/tustison/Utilities/bin64/";

$PID = "run00";

##
# Change image information to have 0x0x0 origin and identity direction
##

$changeInformation = $tmpDirectory . "change_image_information_${PID}.nii.gz";

@args = ( $utilityDirectory . "ChangeImageInformation", 3,
  $uncorrectedHeliumImage, $changeInformation, 0, '0x0x0' );
system( @args ) == 0 or die "system @args failed: $?";
@args = ( $utilityDirectory . "ChangeImageInformation", 3,
  $changeInformation, $changeInformation, 2 );
system( @args ) == 0 or die "system @args failed: $?";


##
# Bias field correction
##

$heliumImage = $tmpDirectory . "/he3_bias_corrected_" . $PID . ".nii.gz";
@args = ( $utilityDirectory . "InhomogeneityCorrectImage", 3,
           '-i', $changeInformation,
           '-s', 1,
           '-c', '[1,0.00001]',
           '-b', '[100,3,0,0]',
           '-t', '[0.1,0.1,200]',
           '-o', "[${heliumImage},${tmpDirectory}/bias_field_${PID}.nii.gz]" );
system( @args ) == 0 or die "system @args failed: $?";


#
# Do a slice by slice vesselness filter to remove vessels/airways
#
$hessianImage = $tmpDirectory . "hessian_mask_${PID}.nii.gz";

@imageInfo = `${utilityDirectory}/GenerateStatisticsFromImage $heliumImage`;
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

for( $i = 0; $i <= $smallestSize; $i++ )
  {
  $sliceImage = "${tmpDirectory}/slice_${i}.nii.gz";
  `${utilityDirectory}/ExtractSliceFromImage $heliumImage $sliceImage $smallestIndex $i`;
  `${utilityDirectory}/HessianBasedFeatures 2 $sliceImage $sliceImage 1 0.5 4 10 0.5 0.5 5 0`;
  `${utilityDirectory}/OtsuThresholdImage 2 $sliceImage $sliceImage 1 200`;
  }
`${utilityDirectory}/ConvertImageSeries $tmpDirectory slice_%d.nii.gz $hessianImage 0 $smallestSize 1`;
`${utilityDirectory}/ChangeImageInformation 3 $hessianImage $hessianImage 4 $heliumImage`;


#!/usr/bin/perl

use File::Basename;

$utilityDirectory = "/home/tustison/Utilities/bin64/";

$image = $ARGV[0];
$imagePrefix = $ARGV[0];
$imagePrefix =~ s/\..*//;

# ( $trash, $directory, $suffix ) = fileparse( $image, qr/\.[^.]*/ );
( $trash, $directory, $suffix ) = fileparse( $image, qr/\.nii\.gz/ );
$imagePrefix = $trash;

# make 2 directories
$inspirationTempDirectory = "${directory}/tmp_dir_inspiration/";
$expirationTempDirectory = "${directory}/tmp_dir_expiration/";

`/bin/mkdir $inspirationTempDirectory`;
`/bin/mkdir $expirationTempDirectory`;

$inspirationCount = 0;
$expirationCount = 0;

for( $i = 0; $i < 28; $i++ )
  {
  if( $i % 2 == 0 ) # inspiration
    {
    `${utilityDirectory}/ExtractSliceFromImage ${image} ${inspirationTempDirectory}/slice${inspirationCount}.nii.gz 2 $i`;
    ${inspirationCount}++;
    }
  else              # expiration
    {
    `${utilityDirectory}/ExtractSliceFromImage ${image} ${expirationTempDirectory}/slice${expirationCount}.nii.gz 2 $i`;
    ${expirationCount}++;
    }
  }
${inspirationCount}--;
${expirationCount}--;

`${utilityDirectory}/ConvertImageSeries $inspirationTempDirectory slice%d.nii.gz ${imagePrefix}_inspiration.nii.gz 0 $inspirationCount 1`;
`${utilityDirectory}/ChangeImageInformation 3 ${imagePrefix}_inspiration.nii.gz ${imagePrefix}_inspiration.nii.gz 1 1.7188x1.7188x15`;

`${utilityDirectory}/ConvertImageSeries $expirationTempDirectory slice%d.nii.gz ${imagePrefix}_expiration.nii.gz 0 $expirationCount 1`;
`${utilityDirectory}/ChangeImageInformation 3 ${imagePrefix}_expiration.nii.gz ${imagePrefix}_expiration.nii.gz 1 1.7188x1.7188x15`;
`${utilityDirectory}/PermuteAxesImage 3 ${imagePrefix}_expiration.nii.gz ${imagePrefix}_expiration.nii.gz 0x2x1`;
`${utilityDirectory}/FlipImage 3 ${imagePrefix}_expiration.nii.gz ${imagePrefix}_expiration.nii.gz 0x0x1`;
`${utilityDirectory}/ChangeImageInformation 3 ${imagePrefix}_expiration.nii.gz ${imagePrefix}_expiration.nii.gz 2 0x0x0`;


`/bin/rm -r $inspirationTempDirectory`;
`/bin/rm -r $expirationTempDirectory`;

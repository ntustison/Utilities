#!/usr/bin/perl -w

use strict;

use Cwd 'realpath';
use Switch;
use File::Find;
use File::Basename;
use File::Path;
use File::Spec;
use FindBin qw($Bin);

my $usage = qq{
  Usage: keyhole.pl <registeredLowResInput> <highResInput> <outputBaseDir> <outputPrefix> <filterDimensions>

      *******************            Overview             *********************

      *******************    Command Line Parameters      *********************

      *******************            References           *********************

  };

if (!($#ARGV + 1))
  {
  print "$usage\n";
  exit 0;
  }
#elsif( $#ARGV > 1 )
#  {
#  die "ERROR: Missing arguments, run without args to see usage\n\t";
#  }

my ( $highRes, $lowRes, $outputBaseDir, $outputPrefix, $dimensions ) = @ARGV;

# Convert I/O directories to absolute paths
$highRes = File::Spec->rel2abs( $highRes );
$lowRes = File::Spec->rel2abs( $lowRes );
$outputBaseDir = File::Spec->rel2abs( $outputBaseDir );

if( ! -d $outputBaseDir )
  {
  mkpath( $outputBaseDir, {verbose => 0, mode => 0755} ) or
    die "Can't create output directory $outputBaseDir\n\t";
  }

my $zerosImage = "${outputBaseDir}/zeros.nii.gz";
my $lowResFTReal = "${outputBaseDir}/lowResFTReal.nii.gz";
my $lowResFTImag = "${outputBaseDir}/lowResFTImag.nii.gz";
my $highResFTReal = "${outputBaseDir}/highResFTReal.nii.gz";
my $highResFTImag = "${outputBaseDir}/highResFTImag.nii.gz";
my $filter = "${outputBaseDir}/filter.nii.gz";

my @info = `GetImageInformation 3 $highRes`;
my $center = $info[4];
$center =~ s/ |Center|:|\[|\]|\r|\n//g;
$center =~ s/,/x/g;


`CreateZeroImage 2 $highRes $zerosImage 0`;
`CreateImageSource 2 $zerosImage $filter 2 ${center} ${dimensions}`;
`FourierTransformImage 2 $lowRes $zerosImage $lowResFTReal $lowResFTImag 0 0 1 0`;
`BinaryOperateImages 2 $lowResFTReal x $filter $lowResFTReal`;
`BinaryOperateImages 2 $lowResFTImag x $filter $lowResFTImag`;
#
`InvertImageIntensity 2 $filter $filter`;
`FourierTransformImage 2 $highRes $zerosImage $highResFTReal $highResFTImag 0 0 1 0`;
`BinaryOperateImages 2 $highResFTReal x $filter $highResFTReal`;
`BinaryOperateImages 2 $highResFTImag x $filter $highResFTImag`;
#
my $outputReal = "${outputBaseDir}/${outputPrefix}real.nii.gz";
my $outputImag = "${outputBaseDir}/${outputPrefix}imag.nii.gz";
#
`BinaryOperateImages 2 $lowResFTReal + $highResFTReal $outputReal`;
`BinaryOperateImages 2 $lowResFTImag + $highResFTImag $outputImag`;
#
`FourierTransformImage 2 $outputReal $outputImag $outputReal $outputImag 1 1 0 1`;

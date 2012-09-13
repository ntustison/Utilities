#! /usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use File::Path;
use File::Spec;

my $ANTSPATH = '/Users/ntustison/Pkg/ANTS/bin/bin/';
my $UTILITIESPATH = '/Users/ntustison/Pkg/Utilities/bin/';

my $usage = qq{

Usage: T1MappingMotionCorrection.pl <input_dir> <output_dir> <output_file_root>

  -------------------------------------------------------------------
  Command line arguments
  -------------------------------------------------------------------

  <input_dir> - input directory.

  <output_dir> - Output base directory.

  <output_file_root> - Root of output, e.g. subject_T1_corrected

  -------------------------------------------------------------------
  References
  -------------------------------------------------------------------

  H. Xue, S. Sha, A. Greiser, et al. "Motion Correction for Myocardial
  T1 Mapping Using Image Registration with Synthetic Image Estimation",
  Magnetic Resonance in Medicine,


Regarding the parameters, the

alpha = 5.0;
beta = 12.0;

There are few other parameters which may be relevant:

// the temporal step size for solving the PDE
deltaT = 0.1;

// maximal number of iterations for solving the PDE
maxPDEIter = 30;

// the threshold for the minimal changes of cost function values
tol = 1e-3;

// the derivative of synthetic images is computed by convolving with the derivative of Gaussian function
// the sigma for the guassian kernel in the unit of pixel
sigmaForDerivative = 1.2;


};


############################################################
############################################################
#
# Motion correction is performed with the following steps:
#  1.  The image sequence is ordered according to the inversion time (TI).
#  1a. Optionally, perform bias correction.
#  2.  The first and last images are registered and then used to
#      estimate the MOLLI signal curve.
#  3.  A signal image, S, for each TI of the image sequence is
#      constructed and made positive.
#  4.  Synthetic images, M, are created by minimizing equation [5] in
#      the reference above.
#  5.  The original images, I, are then pair-wise registered to the
#      synthetic images to create motion corrected images.
#  6.  The motion corrected images are used to fit the three-parameter
#      model.  However, this has to be done with a brute-force search as
#      it is not
#      known which images are negative.  Initially, all images are
#      considered to be positive and then fitted to the curve.  In the
#      next iteration, the first image is considered negative (with
#      all the rest positive) and fitting occurs.  This occurs for
#      all possibilities and the one with the minimum residual is
#      chosen.
#  7.  New signal images are then created.
#  8.  Go back to 4.
#
############################################################
############################################################

my ( $inputDir, $outputDir, $outputPrefix ) = @ARGV;


$inputDir = File::Spec->rel2abs( $inputDir );
$outputDir = File::Spec->rel2abs( $outputDir );

if( ! -d $outputDir )
  {
  mkpath($outputDir, {verbose => 0, mode => 0755}) or die
    "Can't make output directory $outputDir\n\t";
  }


my $alpha = 5.0;
my $beta = 12.0;
my $deltaK = 0.1;
my $numberOfIterations = 1;

my $performBiasCorrection = 0;
my $medianFilterSignalImages = 1;

my @args = ();
my @motionCorrectedImages = ();

############################################################
#
# Step 1: Order the image sequence is ordered according to the inversion
#         time (TI).
#
############################################################
print "** Step 1 **\n";

my $inversionTimeTagID = '0018\|0082';

my @files = <${inputDir}/*>;

my %hash = ();

for( my $i = 0; $i < @files; $i++ )
  {
  my $out = `GetDicomTagValue $files[$i] $inversionTimeTagID`;
  chomp( $out );
  my @outString = split( ',', $out );
  if( $outString[0] =~ m/^${inversionTimeTagID}$/g )
    {
    $hash{ $files[$i] } = scalar( $outString[1] );
    }
  }

# Sort the files based on the inversion time

my @originalImages = ();
my @inversionTimes = ();
for( my $i = 0; $i < @files; $i++ )
  {
  push( @originalImages, $files[$i] );
  push( @inversionTimes, $hash{$files[$i]} );

  my $motionImage = $originalImages[$i];
  $motionImage =~ s/${inputDir}/${outputDir}/;
  $motionImage .= "_motionCorrected.nii.gz";

  push( @motionCorrectedImages, $motionImage );
  }

print "--------------------------------------------------------\n";
print "Original image file name -> inversion time\n";
print "--------------------------------------------------------\n";
for( my $i = 0; $i < @originalImages; $i++ )
  {
  print "$originalImages[$i] -> $inversionTimes[$i]\n";
  }
print "========================================================\n\n";

############################################################
#
# Step 1: Order the image sequence is ordered according to the inversion
#         time (TI).
#
############################################################

my @n4CorrectedImages = ();

if( $performBiasCorrection )
  {
  print "** Step 1a **\n";

  print "--------------------------------------------------------\n";
  print "N4 bias correction\n";

  foreach my $image ( @originalImages )
    {
    my $corrected = $image;
    $corrected =~ s/$inputDir/$outputDir/;
    $corrected .= "_n4.nii.gz";

    push( @n4CorrectedImages, $corrected );

    @args = ( "${ANTSPATH}/N4BiasFieldCorrection",
      '--image-dimensionality', 2,
      '--input-image', $image,
      '--shrink-factor', 2,
      '--convergence', '[50x50x30x20,0.00000001]',
      '--bspline-fitting', '[200]',
      '--output', $corrected
      );
    if( ! -e $corrected )
      {
      system( @args ) == 0 || die "Error:  N4 on $image.\n";
      }
    else
      {
      print "$corrected already exists.\n";
      }
    }
  print "========================================================\n\n";
  }
else
  {
  foreach my $image ( @originalImages )
    {
    my $corrected = $image;
    $corrected =~ s/$inputDir/$outputDir/;
    $corrected .= "_original.nii.gz";

    `ConvertImage 2 $image $corrected 0`;
    push( @n4CorrectedImages, $corrected );
    }
  }

############################################################
#
# Step 2:  The first and last images are registered and then used to
#          estimate the MOLLI signal curve.
#
############################################################
print "** Step 2 **\n";

print "--------------------------------------------------------\n";
print "Registration of first image to last image\n";

my $lastToFirst = "${outputDir}/${outputPrefix}_lastToFirst";

@args = ( "${ANTSPATH}/antsRegistration",
            '-d', 2,
            '-o', "[${lastToFirst},${lastToFirst}Warped.nii.gz]",
            '-u', 0,
#             '-m', "MI[${n4CorrectedImages[0]},${n4CorrectedImages[-1]},1,20]",
#             '-t', 'Affine[0.25]',
#             '-i', '100x100x100',
#             '-s', '1x0.5x0',
#             '-f', '4x2x1',
            '-m', "CC[${n4CorrectedImages[0]},${n4CorrectedImages[-1]},1,6]",
            '-t', 'SyN[0.0001,3.0,0.0]',
            '-c', '0x1',
#             '-t', 'SyN[0.5,3.0,0.0]',
#             '-c', '100x100',
            '-s', '0x0',
            '-f', '2x1'
          );
system( @args ) == 0 || die "Error:  antsRegistration on first and last images.\n";


my $negativeFirstImage = "${outputDir}/${outputPrefix}_negativeFirstImage.nii.gz";

`${ANTSPATH}/ImageMath 2 $negativeFirstImage m ${n4CorrectedImages[0]} -1`;

@args = ( "${UTILITIESPATH}/SalernoFitVoxelwise3ParameterModel",
  "${outputDir}/${outputPrefix}",
  $negativeFirstImage,
  $inversionTimes[0],
  "${lastToFirst}Warped.nii.gz",
  $inversionTimes[-1]
  );
system( @args ) == 0 || die "Error: @args.\n";

my $Aimage = "${outputDir}/${outputPrefix}A.nii.gz";
my $Bimage = "${outputDir}/${outputPrefix}B.nii.gz";
my $T1image = "${outputDir}/${outputPrefix}T1.nii.gz";

print "========================================================\n\n";

############################################################
#
# Step 3. A signal image, S, for each TI of the image sequence is
#         constructed and made positive.
#
############################################################

print "** Step 2 **\n";

print "--------------------------------------------------------\n";
print "Generating signal images.\n";

my @signalImages = ();

for( my $i = 0; $i < @inversionTimes; $i++ )
  {
  my $idx = sprintf( "%05d", $i );

  print "  Generating signal image $idx (TI = $inversionTimes[$i])\n";

  push( @signalImages, "${outputDir}/${outputPrefix}_SignalImage${idx}.nii.gz" );
  @args = ( "${UTILITIESPATH}/SalernoGenerateSignalImage",
    $Aimage,
    $Bimage,
    $T1image,
    $inversionTimes[$i],
    $signalImages[$i]
    );
  system( @args ) == 0 || die "Error: Step 3.\n";

  if( $medianFilterSignalImages )
    {
    `SmoothImage 2 $signalImages[$i] 1 $signalImages[$i] 0 1`;
    }

#   `${ANTSPATH}/ImageMath 2 $signalImages[$i] abs $signalImages[$i]`;
  }

print "========================================================\n\n";


############################################################
#
# Begin loop
#
############################################################

my @inputImages = @n4CorrectedImages;

for( my $iteration = 0; $iteration <= $numberOfIterations; $iteration++ )
  {
  print "Iteration $iteration\n";

  ############################################################
  #
  # Step 4. Synthetic images, M, are created by minimizing equation [5] in
  #         the reference above.
  #
  ############################################################

  print "  ** Step 4 **\n";

  print "  --------------------------------------------------------\n";
  print "  Generating synthetic images.\n";

  my @parameters = ();
  for( my $i = 0; $i < @inputImages; $i++ )
    {
    push( @parameters, $inputImages[$i] );
    push( @parameters, $signalImages[$i] );
    }

  @args = ( "${UTILITIESPATH}/SalernoSyntheticImageEstimation",
    "${outputDir}/${outputPrefix}_SyntheticImage",
    $deltaK,
    $alpha,
    $beta,
    @parameters
    );
  system( @args ) == 0 || die "Error: Step 4.\n";

  my @syntheticImages = <${outputDir}/${outputPrefix}_SyntheticImage*nii.gz>;

  for( my $i = 0; $i < @syntheticImages; $i++ )
    {
    `${ANTSPATH}/ImageMath 2 $syntheticImages[$i] abs $syntheticImages[$i]`;
    }

  ############################################################
  #
  # 5.  The original images, I, are then pair-wise registered to the
  #     synthetic images to create motion corrected images.
  #
  ############################################################

  for( my $i = 0; $i < @n4CorrectedImages; $i++ )
    {
    my $motionXfrm = $originalImages[$i];
    $motionXfrm =~ s/\.nii\.gz//;
    $motionXfrm =~ s/$inputDir/$outputDir/;

    @args = ( "${ANTSPATH}/antsRegistration",
                '-d', 2,
                '-o', "[${motionXfrm},${motionCorrectedImages[$i]}]",
                '-u', 0,
#                 '-w', '[0.025,0.0975]',
#                 '-m', "MI[${syntheticImages[$i]},${n4CorrectedImages[$i]},1,20]",
#                 '-t', 'Affine[0.25]',
#                 '-c', '100x100x100',
#                 '-s', '1x0.5x0',
#                 '-f', '4x2x1',
                '-m', "CC[${syntheticImages[$i]},${n4CorrectedImages[$i]},1,6]",
#                 '-t', 'SyN[0.25,3.0,0.0]',
#                 '-c', '100x100',
            '-t', 'SyN[0.0001,3.0,0.0]',
            '-c', '0x1',
                '-s', '0x0',
                '-f', '2x1'
              );
    system( @args ) == 0 || die "Error:  antsRegistration on pairs.\n";

    $inputImages[$i] = $motionCorrectedImages[$i];
    }

  ############################################################
  #
  #  6.  The motion corrected images are used to fit the three-parameter
  #      model.  However, this has to be done with a brute-force search as
  #      it is not
  #      known which images are negative.  Initially, all images are
  #      considered to be positive and then fitted to the curve.  In the
  #      next iteration, the first image is considered negative (with
  #      all the rest positive) and fitting occurs.  This occurs for
  #      all possibilities and the one with the minimum residual is
  #      chosen.
  #
  ############################################################

  my $minResidual = 1e10;
  my $minResidualIndex = -1;

  # only try the first 5 images;
  my $numberOfImagesForSearch = 5;
  if( $numberOfImagesForSearch > @inputImages )
    {
    $numberOfImagesForSearch = @inputImages;
    }

  for( my $i = 0; $i < $numberOfImagesForSearch; $i++ )
    {
    my @signedInputImages = @inputImages;

    for( my $j = 0; $j < $i; $j++ )
      {
      my $negativeImage = $inputImages[$j];
      $negativeImage =~ s/\.nii\.gz/Negative\.nii\.gz/;
      if( ! -e $negativeImage )
        {
        `${ANTSPATH}/ImageMath 2 $negativeImage m ${inputImages[$j]} -1`;
        }
      $signedInputImages[$j] = $negativeImage;
      }

    my @parameters = ();
    for( my $j = 0; $j < @inputImages; $j++ )
      {
      push( @parameters, $signedInputImages[$j] );
      push( @parameters, $inversionTimes[$j] );
      }

    @args = ( "${UTILITIESPATH}/SalernoFitVoxelwise3ParameterModel",
      "${outputDir}/${outputPrefix}",
      @parameters
      );
    system( @args ) == 0 || die "Error: @args.\n";

    my $currentResidual = 0;

    for( my $j = 0; $j < @inversionTimes; $j++ )
      {
      my $idx = sprintf( "%05d", $j );

      @args = ( "${UTILITIESPATH}/SalernoGenerateSignalImage",
        $Aimage,
        $Bimage,
        $T1image,
        $inversionTimes[$j],
        $signalImages[$j]
        );
      system( @args ) == 0 || die "Error: @args.\n";

      if( $medianFilterSignalImages )
        {
        `SmoothImage 2 $signalImages[$j] 1 $signalImages[$j] 0 1`;
        }

      my $tmpResidualImage = "${outputDir}/${outputPrefix}_tmpResidual.nii.gz";

      `${UTILITIESPATH}/BinaryOperateImages 2 $signalImages[$j] - $signedInputImages[$j] $tmpResidualImage`;
      `${UTILITIESPATH}/UnaryOperateImage 2 $tmpResidualImage ^ 2 $tmpResidualImage`;

      my @out = `${UTILITIESPATH}/CalculateFirstOrderStatisticsFromImage 2 $tmpResidualImage`;
      my @stats = split( ' ', $out[1] );

      $currentResidual += $stats[0];

      unlink( $tmpResidualImage );
      }



    if( $currentResidual < $minResidual )
      {
      $minResidual = $currentResidual;
      $minResidualIndex = $i;
      }
    }

  print "$minResidualIndex: $minResidual\n";

  ############################################################
  #
  #  7.  New signal images are created.
  #
  ############################################################

  my @signedInputImages = @inputImages;

  for( my $i = 0; $i <= $minResidualIndex; $i++ )
    {
    my $negativeImage = $inputImages[$i];
    $negativeImage =~ s/\.nii\.gz/Negative\.nii\.gz/;
    if( ! -e $negativeImage )
      {
      `${ANTSPATH}/ImageMath 2 $negativeImage m ${inputImages[$i]} -1`;
      }
    $signedInputImages[$i] = $negativeImage;
    }

  @parameters = ();
  for( my $i = 0; $i < @inputImages; $i++ )
    {
    push( @parameters, $signedInputImages[$i] );
    push( @parameters, $inversionTimes[$i] );
    }

  print "@signedInputImages\n";

  @args = ( "${UTILITIESPATH}/SalernoFitVoxelwise3ParameterModel",
    "${outputDir}/${outputPrefix}",
    @parameters
    );
  system( @args ) == 0 || die "Error: @args.\n";

  for( my $i = 0; $i < @inversionTimes; $i++ )
    {
    my $idx = sprintf( "%05d", $i );

    @args = ( "${UTILITIESPATH}/SalernoGenerateSignalImage",
      $Aimage,
      $Bimage,
      $T1image,
      $inversionTimes[$i],
      $signalImages[$i]
      );
    system( @args ) == 0 || die "Error: @args.\n";

    if( $medianFilterSignalImages )
      {
      `SmoothImage 2 $signalImages[$i] 1 $signalImages[$i] 0 1`;
      }
    }
  print "  ========================================================\n\n";
  }


#!/user/bin/perl -w

$outputDirectory = "/home/tustison/Data/HeliumLungStudies/DataSet11_deLangeBaselineAsthma/ImageVolumes/";
$outputFile = "${outputDirectory}/defectsMNDS.csv";

$baseDirectory = "/home/tustison/Data/HeliumLungStudies/DataSet11_deLangeBaselineAsthma/ImageVolumes/";

@names = ( 413, 414, 416, 417, 418, 419, 420, #421,
422, 423, 424, 429, 432, #436, 437,
438, 439, 441, 442, 443, 444, 446, 447, 448, #449, 450,
451, 452, 453, 454, 456, 458, 459, 460, 461, 462, 466, 467, 468, #470,
471, #473, 474, 475,
480, 484, 485, #486,
487, #489, 491,
494, #495,
496, 497, 498 #, 504,
# 505, 515, 531, 533, 538, 549, 550, 553, 568, 572, 620, 621, 635, 637, 638, 639,
# 645, 646, 648, 649, 660, 669, 670, 672, 674, 677, 678, 684, 699, 703
);


# @diagnosis = ( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#  0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0,
#  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#  0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0 );

@begSlices = ( 3, 3, 5, 2, 5, 5, 6, #4,
3, 3, 3, 5, 5, #3, 7,
4, 3, 4, 3, 4, 4, 5, 4, 4, #4,
4, 4, 6, 4, 4, 5, 5, 3, 4, 4, 3, 5, 5, #2,
5, #3, 4,
4, 4, 5,
4, #5, 4,
4, #5,
3, 4, 3 #, 0
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
);

@endSlices = ( 23, 38, 24, 25, 25, 25, 25, #27,
28, 27, 25, 28, 28, #28, 25,
25, 29, 25, 26, 22, 25, 28, 26, 28, #24,
25, 27, 23, 24, 27, 23, 26, 25, 26, 24, 23, 25, 27, #28,
24, #28, 27,
28, 27, 30,
30, #24, 25,
25, #28,
26, 23, 28 #,0
#0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
);

# for( $i = 0; $i < @names; $i++ )
#   {
#   print "${names[$i]} ${begSlices[$i]} ${endSlices[$i]}\n";
#   }
# exit( 0 );

# @slices = ( 20, 28, 21, 24, 22, 23, 24, 23, 22, 22, 24, 25, 23, 26, 22, 25, 19,
# 22, 25, 24, 25, 25, 22, 25, 20, 22, 24, 21, 23, 24, 24, 20, 21, 23, 24, 19, 24,
# 25, 23, 27, 24, 21, 21, 21, 26, 0,
# 20,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );


open( DEF, ">${outputFile}" );

for( $i = 0; $i < @names; $i++ )
  {
  if( ! -e "${baseDirectory}/He3\#${names[${i}]}" )
    {
    print DEF "${names[$i]},\n";
    }
  else
    {
    @directories = <${baseDirectory}/He3\#${names[${i}]}/Series*>;

    open( DEF_IN, "${directories[0]}/defectsSlice.txt" );
    @defects = <DEF_IN>;
    close( DEF_IN );

    $begin = $begSlices[$i] - 1;
    $end = $endSlices[$i] - 1;

#     $begin = int( 0.5 * ( @defects - $slices[$i] ) );
#     $end = $begin + $slices[$i];
#
#     if( $begin < 0 )
#       {
#       print "ERROR: ${names[$i]}\n";
#       print DEF "${names[$i]},\n";
#       next;
#       }

    $sum = 0;
    for( $s = $begin; $s < $end; $s++ )
      {
      @ndefects = split( / /, $defects[$s] );
      $sum += $ndefects[1];
      }
    $mnds = $sum / ( $end - $begin );
    print DEF "${names[$i]},${mnds}\n";
    }
  }

close( DEF );



# for( $i = 0; $i < @names; $i++ )
#   {
#
#
#
# $movingImage = "${outputDirectory}/initial_segmentation_run02.nii.gz";
#
# $output = "${outputDirectory}/antsToTemplate.nii.gz";
#
# @args = ( "${antsDirectory}/ANTS", 3,
#         "-m", "PR[${shapeIntensityTemplateSegmentation},${movingImage},1,4]",
#         "-t", "SyN[0.5]",
#         "-i", "100x100x100",
#         "-r", "Gauss[3,1]",
#         "-o", "${output}.nii.gz" );
# system( @args ) == 0 or die "system @args failed: $?";
#
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
#
#
#
#
#

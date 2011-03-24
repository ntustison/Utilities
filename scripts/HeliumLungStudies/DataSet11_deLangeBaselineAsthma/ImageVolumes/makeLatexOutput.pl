#!/user/bin/perl -w

$outputDirectory = "/home/tustison/Data/HeliumLungStudies/DataSet11_deLangeBaselineAsthma/ImageVolumes/LatexOutput/";
$outputFile = "${outputDirectory}/output.tex";
$figureDirectory = "${outputDirectory}/Figures/";

$baseDirectory = "/home/tustison/Data/HeliumLungStudies/DataSet11_deLangeBaselineAsthma/ImageVolumes/";

@names = ( 413, 414, 416, 417, 418, 419, 420, 422, 423, 424, 429, 432, 438, 439,
441, 442, 443, 444, 446, 447, 448, 449, 451, 452, 453, 454, 456, 458, 459, 460,
461, 462, 466, 467, 468, 471, 475, 480, 484, 485, 487, 494, 496, 497, 498, 504,
505, 515, 531, 533, 538, 549, 550, 553, 568, 572, 620, 621, 635, 637, 638, 639,
645, 646, 648, 649, 660, 669, 670, 672, 674, 677, 678, 684, 699, 703 );

@diagnosis = ( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0,
 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0 );


open( FILE, ">${outputFile}" );

print FILE "\\documentclass[11pt]{article}\n";
print FILE "\\usepackage{geometry}\n";
print FILE "\\geometry{letterpaper,top=0.5in,left=0.5in,bottom=0.5in,top=0.5in,headsep=6pt,footskip=18pt}\n";
print FILE "\\usepackage{fancyhdr}\n";
print FILE "\\pagestyle{empty}\n";
print FILE "\\renewcommand{\\headrulewidth}{0.75pt}\n";
print FILE "\\renewcommand{\\footrulewidth}{0.75pt}\n";
print FILE "\\usepackage{helvet}\n";
print FILE "\\renewcommand{\\familydefault}{\\sfdefault}\n";
print FILE "\\renewcommand{\\baselinestretch}{.9}\n\n";
print FILE "\\usepackage{graphicx}\n";
print FILE "\\graphicspath{{/Users/nick/Desktop/LatexOutput/Figures/}}\n\n";
print FILE "\\renewcommand{\\topfraction}{0.85}\n";
print FILE "\\renewcommand{\\textfraction}{0.1}\n";
print FILE "\\renewcommand{\\floatpagefraction}{0.75}\n\n";
print FILE "\\begin{document}\n\n";

for( $i = 0; $i < @names; $i++ )
  {
  if( $i % 4 == 0 )
    {
    if( $i > 0 )
      {
      print FILE "\\end{tabular}\n\\end{center}\n\\end{figure}\n\n";
      print FILE "\\clearpage\n\n";
      }
    print FILE "\n\\begin{figure}\n\\begin{center}\n\\begin{tabular}{ccc}\n";
    }
  if( $diagnosis[$i] == 1 )
    {
    print FILE "\\rotatebox{90}{\\centering\\parbox{47.5mm}{\\centering \\bf \$\^3\$He ${names[${i}]} (Asthmatic)}} & \n";
    }
  else
    {
    print FILE "\\rotatebox{90}{\\centering\\parbox{47.5mm}{\\centering \\bf \$\^3\$He ${names[${i}]} (Normal)}} & \n";
    }

  @directories = <${baseDirectory}/He3\#${names[${i}]}/Series*>;
  @originalpaths = <${directories[0]}/latex_originalslice*.png>;
  @slicepaths = <${directories[0]}/latex_slice*.png>;
  @atropospaths = <${directories[0]}/latex_atropos*.png>;

  @originalcomponents = split( /\//, $originalpaths[0] );
  @slicecomponents = split( /\//, $slicepaths[0] );
  @atroposcomponents = split( /\//, $atropospaths[0] );

  print $originalcomponents[@originalcomponents-1] . "\n";
  print $atroposcomponents[@atroposcomponents-1] . "\n";
  print $slicecomponents[@slicecomponents-1] . "\n";

  $originalfile = "He3_${names[${i}]}$originalcomponents[@slicecomponents-1]";
  $slicefile = "He3_${names[${i}]}$slicecomponents[@slicecomponents-1]";
  $atroposfile = "He3_${names[${i}]}$atroposcomponents[@atroposcomponents-1]";

  `/bin/cp ${slicepaths[0]} ${figureDirectory}/${slicefile}`;
  `/bin/cp ${originalpaths[0]} ${figureDirectory}/${originalfile}`;
  `/bin/cp ${atropospaths[0]} ${figureDirectory}/${atroposfile}`;

  print FILE "\\includegraphics[width=90mm]{${slicefile}} & \n";
  print FILE "\\includegraphics[width=90mm]{${atroposfile}} \\\\ \n";
  }

print FILE "\\end{tabular}\n\\end{center}\n\\end{figure}\n\n";
print FILE "\n\n\\end{document}\n";

close( FILE );



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

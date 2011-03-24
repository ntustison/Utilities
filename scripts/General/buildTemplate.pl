#!/usr/bin/perl

$utilitiesDirectory = "/mnt/data2/Avants/bin/3D/";

$imageDirectory = "/mnt/data2/PUBLIC/Data/Input/TedCareyUltrasound/high/G3 high/";
$outputDirectory = "/mnt/data1/tustison/Projects/UltrasoundRegistration/";

for $i ( 0..4 )
  {  
  opendir( DIR, $imageDirectory );
  @pngfiles = grep( /\.png$/, readdir(DIR) );
  closedir( DIR );
  
  $count = 0;
  $template = "";
  foreach $file( @pngfiles )
    {
    print "(" . $i . "): handling file " . $file . "\n";    
    
    if ( $count == 0 )
      {
      $template = $outputDirectory . "template.nii";
      
      @args = ( "/mnt/data1/tustison/Utilities/bin/2D/float/ConvertImage", 
             $imageDirectory . $file, $template ); 
      system( @args ) == 0 or die "system @args =  failed: $?"
      }
    else
      {  
      $movingFile = $outputDirectory . "moving" . $count . ".nii";
      @args = ( "/mnt/data1/tustison/Utilities/bin/3D/float/ConvertImage", 
             $imageDirectory . $file, $movingFile ); 
      system( @args ) == 0 or die "system @args =  failed: $?";

      $algorithm = $utilitiesDirectory . "SyMMN";
      @args = ( $algorithm, "-f", $template, 
                         "-m", $movingFile, 
                         "-l", 0.25, 
                         "-s", 3,
                         "-c", 3,
                         "-o", $outputDirectory . "output" . $count,
                         "-n", 5,
                         "-i", "41x41x41x31x7",
                         "-a", 0.001
  		 																				);
      system( @args ) == 0 or die "system @args =  failed: $?";
      } 
    
    $count = $count + 1;
    }
  @args = ( "rm", $outputDirectory . "*inversewarp*" );
  system( @args ) == 0 or die "system @args =  failed: $?";
  @args = ( "rm", $outputDirectory . "templatewarp*" );
  system( @args ) == 0 or die "system @args =  failed: $?";
    
  $algorithm = $utilitiesDirectory . "AverageImages";
  @args = ( $algorithm, $outputDirectory . "templatewarpxvec.nii", 0, 
         $outputDirectory . "*warpxvec.nii" );
  system( @args ) == 0 or die "system @args =  failed: $?";
  $algorithm = $utilitiesDirectory . "AverageImages";
  @args = ( $algorithm, $outputDirectory . "templatewarpyvec.nii", 0, 
         $outputDirectory . "*warpyvec.nii" );
  system( @args ) == 0 or die "system @args =  failed: $?";
  $algorithm = $utilitiesDirectory . "DiffeomorphicFitToDeformationField";
  @args = ( $algorithm, 0, $outputDirectory . "templatewarp", 
                        $outputDirectory . "template.nii", 
                        $outputDirectory . "templateOut", 
                        $outputDirectory . "template.nii", 15 );
  system( @args ) == 0 or die "system @args =  failed: $?";

  }
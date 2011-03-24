#!/user/bin/perl -w

use File::Find;
use File::Basename;

$directory = $ARGV[0];

find( \&wanted, $directory );

$count = 0;

sub wanted
  {

  ( $filename, $directories, $suffix ) = fileparse( $File::Find::name );
  @current_vars = split( /\//, $directories ); 

  if( $filename =~ m/^Series/ && $filename =~ m/.nii.gz$/  )
    {
    print $directories . $filename . "\n"; 
    if( $count < 10 )
      {
      @args = ( "cp",  $directories . $filename, 
        "/home/tustison/Data/HeliumLungStudies/PrePostAlbuterol/ImageVolumes/PCA/CF/he3_0${count}.nii.gz" );
      system( @args ) == 0 or die "system @args failed: $?";
      }
    else
      {
      @args = ( "cp",  $directories . $filename, 
        "/home/tustison/Data/HeliumLungStudies/PrePostAlbuterol/ImageVolumes/PCA/CF/he3_${count}.nii.gz" );
      system( @args ) == 0 or die "system @args failed: $?";
      }
    $count++;  
      
#      unlink glob( $directories . $filename . '/*run0*' );
#      unlink glob( $directories . $filename . '/*Warp*' );
#      unlink glob( $directories . $filename . '/*Affine*' );
#      unlink glob( $directories . $filename . '/registrationToTemplate*' );
#      unlink glob( $directories . $filename . '/points.vtk' );
#      unlink glob( $directories . $filename . '/*_hot.mha' );
     
#      $image = $directories . $filename . "/" . $filename . ".nii.gz"; 
#      $segmentation = $directories . $filename . "/segmentation.nii.gz"; 
#        
#      @args = ( "snap", "-g", $image, "-s", $segmentation );
#      system( @args ) == 0 or die "system @args failed: $?";
#      print "@args\n";
    }
  }


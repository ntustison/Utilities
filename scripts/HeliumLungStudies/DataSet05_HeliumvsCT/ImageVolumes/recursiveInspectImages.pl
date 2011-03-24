#!/user/bin/perl -w

use File::Find;
use File::Basename;

$directory = $ARGV[0];

find( \&wanted, $directory );


sub wanted
  {

  ( $filename, $directories, $suffix ) = fileparse( $File::Find::name );
  @current_vars = split( /\//, $directories );

  if( -d && $filename =~ m/^Series/ && $filename !~ m/\.nii\.gz$/ )
    {
    eval
      {
      print $directories . $filename . "/" . "\n";
#       unlink glob( $directories . $filename . '/*run0*' );
#       unlink glob( $directories . $filename . '/*Warp*' );
#       unlink glob( $directories . $filename . '/*Affine*' );
#       unlink glob( $directories . $filename . '/registrationToTemplate*' );
#       unlink glob( $directories . $filename . '/points.vtk' );
#       unlink glob( $directories . $filename . '/*_hot.mha' );

       $directory = "/home/tustison/Data/HeliumLungStudies/DataSet05_HeliumvsCT/ImageVolumes/";

       $which = @current_vars - 1;

 ############################
 #  Do two class
 ###########################

       open( TWOCLASS_FILE, ">>${directory}/TwoClassData.cvs" );

       $name = $current_vars[$which];
       print TWOCLASS_FILE "${name},";

       ## Do whole lung
       open( NONVENTILATED_FILE, "${directories}/${filename}/apocrita_2class_ventilation_ratio_both_1_run03.txt" );
       @fileContent = <NONVENTILATED_FILE>;
       close( NONVENTILATED_FILE );

       chomp( $fileContent[1] );
       print TWOCLASS_FILE "${fileContent[1]},";

       open( VENTILATED_FILE, "${directories}/${filename}/apocrita_2class_ventilation_ratio_both_2_run03.txt" );
       @fileContent = <VENTILATED_FILE>;
       close( VENTILATED_FILE );

       chomp( $fileContent[1] );
       print TWOCLASS_FILE "${fileContent[1]},";

       ## Do left lung
       open( NONVENTILATED_FILE, "${directories}/${filename}/apocrita_2class_ventilation_ratio_both_1_left_run03.txt" );
       @fileContent = <NONVENTILATED_FILE>;
       close( NONVENTILATED_FILE );

       chomp( $fileContent[1] );
       print TWOCLASS_FILE "${fileContent[1]},";

       open( VENTILATED_FILE, "${directories}/${filename}/apocrita_2class_ventilation_ratio_both_2_left_run03.txt" );
       @fileContent = <VENTILATED_FILE>;
       close( VENTILATED_FILE );

       chomp( $fileContent[1] );
       print TWOCLASS_FILE "${fileContent[1]},";

       ## Do right lung
       open( NONVENTILATED_FILE, "${directories}/${filename}/apocrita_2class_ventilation_ratio_both_1_right_run03.txt" );
       @fileContent = <NONVENTILATED_FILE>;
       close( NONVENTILATED_FILE );

       chomp( $fileContent[1] );
       print TWOCLASS_FILE "${fileContent[1]},";

       open( VENTILATED_FILE, "${directories}/${filename}/apocrita_2class_ventilation_ratio_both_2_right_run03.txt" );
       @fileContent = <VENTILATED_FILE>;
       close( VENTILATED_FILE );

       chomp( $fileContent[1] );
       print TWOCLASS_FILE "${fileContent[1]}\n";

       close( TWOCLASS_FILE );



 ############################
 #  Do three class
 ###########################

       open( THREECLASS_FILE, ">>${directory}/ThreeClassData.cvs" );

       $name = $current_vars[$which];
       print THREECLASS_FILE "${name},";


       ## Do whole lung
       open( NONVENTILATED_FILE, "${directories}/${filename}/apocrita_3class_ventilation_ratio_both_1_run03.txt" );
       @fileContent = <NONVENTILATED_FILE>;
       close( NONVENTILATED_FILE );

       chomp( $fileContent[1] );
       print THREECLASS_FILE "${fileContent[1]},";

       open( VENTILATED_FILE, "${directories}/${filename}/apocrita_3class_ventilation_ratio_both_2_run03.txt" );
       @fileContent = <VENTILATED_FILE>;
       close( VENTILATED_FILE );

       chomp( $fileContent[1] );
       print THREECLASS_FILE "${fileContent[1]},";

       open( VENTILATED_FILE, "${directories}/${filename}/apocrita_3class_ventilation_ratio_both_3_run03.txt" );
       @fileContent = <VENTILATED_FILE>;
       close( VENTILATED_FILE );

       chomp( $fileContent[1] );
       print THREECLASS_FILE "${fileContent[1]},";

       ## Do left lung
       open( NONVENTILATED_FILE, "${directories}/${filename}/apocrita_3class_ventilation_ratio_both_1_left_run03.txt" );
       @fileContent = <NONVENTILATED_FILE>;
       close( NONVENTILATED_FILE );

       chomp( $fileContent[1] );
       print THREECLASS_FILE "${fileContent[1]},";

       open( VENTILATED_FILE, "${directories}/${filename}/apocrita_3class_ventilation_ratio_both_2_left_run03.txt" );
       @fileContent = <VENTILATED_FILE>;
       close( VENTILATED_FILE );

       chomp( $fileContent[1] );
       print THREECLASS_FILE "${fileContent[1]},";

       open( VENTILATED_FILE, "${directories}/${filename}/apocrita_3class_ventilation_ratio_both_3_left_run03.txt" );
       @fileContent = <VENTILATED_FILE>;
       close( VENTILATED_FILE );

       chomp( $fileContent[1] );
       print THREECLASS_FILE "${fileContent[1]},";

       ## Do right lung
       open( NONVENTILATED_FILE, "${directories}/${filename}/apocrita_3class_ventilation_ratio_both_1_right_run03.txt" );
       @fileContent = <NONVENTILATED_FILE>;
       close( NONVENTILATED_FILE );

       chomp( $fileContent[1] );
       print THREECLASS_FILE "${fileContent[1]},";

       open( VENTILATED_FILE, "${directories}/${filename}/apocrita_3class_ventilation_ratio_both_2_right_run03.txt" );
       @fileContent = <VENTILATED_FILE>;
       close( VENTILATED_FILE );

       chomp( $fileContent[1] );
       print THREECLASS_FILE "${fileContent[1]},";

       open( VENTILATED_FILE, "${directories}/${filename}/apocrita_3class_ventilation_ratio_both_3_right_run03.txt" );
       @fileContent = <VENTILATED_FILE>;
       close( VENTILATED_FILE );

       chomp( $fileContent[1] );
       print THREECLASS_FILE "${fileContent[1]}\n";

       close( THREECLASS_FILE );













#      $image = $directories . $filename . "/" . $filename . ".nii.gz";
#      $segmentation = $directories . $filename . "/segmentation.nii.gz";
#
#      @args = ( "snap", "-g", $image, "-s", $segmentation );
#      system( @args ) == 0 or die "system @args failed: $?";
#      print "@args\n";
      }
    }
  }


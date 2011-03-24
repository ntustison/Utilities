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
#     `/home/tustison/Utilities/bin64/CreateZeroImage 3 ${directories}/${filename}/${filename}.nii.gz ${directories}/${filename}/ones_tmp.nii.gz 1`;
#
#     @fileContent = `/home/tustison/Utilities/bin64/CalculateFirstOrderStatisticsFromImage 3 ${directories}/${filename}/${filename}.nii.gz ${directories}/${filename}/ones_tmp.nii.gz 1 200 ${directories}/${filename}/hist_tmp.txt`;
#     @statistics = split( / /, $fileContent[1] );
#    print $statistics[9] . "\n";

#     unlink( "${directories}/${filename}/ones_tmp.nii.gz" );
#     unlink( "${directories}/${filename}/hist_tmp.txt" );

#     if( $statistics[9] > 200 )
#       {
      print $directories . $filename . "/" . "\n";

      $directory = "/home/tustison/Data/HeliumLungStudies/DataSet11_deLangeBaselineAsthma/ImageVolumes/";

      $which = @current_vars - 1;

 ############################
 #  Do two class
 ###########################

      open( TWOCLASS_FILE, ">>${directory}/TwoClassData.csv" );

      $name = $current_vars[$which];
      print TWOCLASS_FILE "${name},";

      ## Do whole lung
      $file = "${directories}/${filename}/apocrita_2class_ventilation_ratio_both_1_run02.txt";
      open( NONVENTILATED_FILE, $file ) or die "Unable to open ${file}: $!";
      @fileContent = <NONVENTILATED_FILE>;
      close( NONVENTILATED_FILE );

      chomp( $fileContent[1] );
      print TWOCLASS_FILE "${fileContent[1]},";

      $file = "${directories}/${filename}/apocrita_2class_ventilation_ratio_both_2_run02.txt";
      open( VENTILATED_FILE, $file ) or die "Unable to open ${file}: $!";
      @fileContent = <VENTILATED_FILE>;
      close( VENTILATED_FILE );

      chomp( $fileContent[1] );
      print TWOCLASS_FILE "${fileContent[1]},";

      ## Do left lung
      $file = "${directories}/${filename}/apocrita_2class_ventilation_ratio_both_1_left_run02.txt";
      open( NONVENTILATED_FILE, $file ) or die "Unable to open ${file}: $!";
      @fileContent = <NONVENTILATED_FILE>;
      close( NONVENTILATED_FILE );

      chomp( $fileContent[1] );
      print TWOCLASS_FILE "${fileContent[1]},";

      $file = "${directories}/${filename}/apocrita_2class_ventilation_ratio_both_2_left_run02.txt";
      open( VENTILATED_FILE, $file ) or die "Unable to open ${file}: $!";
      @fileContent = <VENTILATED_FILE>;
      close( VENTILATED_FILE );

      chomp( $fileContent[1] );
      print TWOCLASS_FILE "${fileContent[1]},";

      ## Do right lung
      $file = "${directories}/${filename}/apocrita_2class_ventilation_ratio_both_1_right_run02.txt";
      open( NONVENTILATED_FILE, $file ) or die "Unable to open ${file}: $!";
      @fileContent = <NONVENTILATED_FILE>;
      close( NONVENTILATED_FILE );

      chomp( $fileContent[1] );
      print TWOCLASS_FILE "${fileContent[1]},";

      $file = "${directories}/${filename}/apocrita_2class_ventilation_ratio_both_2_right_run02.txt";
      open( VENTILATED_FILE, $file ) or die "Unable to open ${file}: $!";
      @fileContent = <VENTILATED_FILE>;
      close( VENTILATED_FILE );

      chomp( $fileContent[1] );
      print TWOCLASS_FILE "${fileContent[1]}\n";

      close( TWOCLASS_FILE );



 ############################
 #  Do three class
 ###########################

#       open( THREECLASS_FILE, ">>${directory}/ThreeClassData.csv" );
#
#       $name = $current_vars[$which];
#       print THREECLASS_FILE "${name},";
#
#
#       ## Do whole lung
#       open( NONVENTILATED_FILE, "${directories}/${filename}/apocrita_3class_ventilation_ratio_both_1_run02.txt" );
#       @fileContent = <NONVENTILATED_FILE>;
#       close( NONVENTILATED_FILE );
#
#       chomp( $fileContent[1] );
#       print THREECLASS_FILE "${fileContent[1]},";
#
#       open( VENTILATED_FILE, "${directories}/${filename}/apocrita_3class_ventilation_ratio_both_2_run02.txt" );
#       @fileContent = <VENTILATED_FILE>;
#       close( VENTILATED_FILE );
#
#       chomp( $fileContent[1] );
#       print THREECLASS_FILE "${fileContent[1]},";
#
#       open( VENTILATED_FILE, "${directories}/${filename}/apocrita_3class_ventilation_ratio_both_3_run02.txt" );
#       @fileContent = <VENTILATED_FILE>;
#       close( VENTILATED_FILE );
#
#       chomp( $fileContent[1] );
#       print THREECLASS_FILE "${fileContent[1]},";
#
#       ## Do left lung
#       open( NONVENTILATED_FILE, "${directories}/${filename}/apocrita_3class_ventilation_ratio_both_1_left_run02.txt" );
#       @fileContent = <NONVENTILATED_FILE>;
#       close( NONVENTILATED_FILE );
#
#       chomp( $fileContent[1] );
#       print THREECLASS_FILE "${fileContent[1]},";
#
#       open( VENTILATED_FILE, "${directories}/${filename}/apocrita_3class_ventilation_ratio_both_2_left_run02.txt" );
#       @fileContent = <VENTILATED_FILE>;
#       close( VENTILATED_FILE );
#
#       chomp( $fileContent[1] );
#       print THREECLASS_FILE "${fileContent[1]},";
#
#       open( VENTILATED_FILE, "${directories}/${filename}/apocrita_3class_ventilation_ratio_both_3_left_run02.txt" );
#       @fileContent = <VENTILATED_FILE>;
#       close( VENTILATED_FILE );
#
#       chomp( $fileContent[1] );
#       print THREECLASS_FILE "${fileContent[1]},";
#
#       ## Do right lung
#       open( NONVENTILATED_FILE, "${directories}/${filename}/apocrita_3class_ventilation_ratio_both_1_right_run02.txt" );
#       @fileContent = <NONVENTILATED_FILE>;
#       close( NONVENTILATED_FILE );
#
#       chomp( $fileContent[1] );
#       print THREECLASS_FILE "${fileContent[1]},";
#
#       open( VENTILATED_FILE, "${directories}/${filename}/apocrita_3class_ventilation_ratio_both_2_right_run02.txt" );
#       @fileContent = <VENTILATED_FILE>;
#       close( VENTILATED_FILE );
#
#       chomp( $fileContent[1] );
#       print THREECLASS_FILE "${fileContent[1]},";
#
#       open( VENTILATED_FILE, "${directories}/${filename}/apocrita_3class_ventilation_ratio_both_3_right_run02.txt" );
#       @fileContent = <VENTILATED_FILE>;
#       close( VENTILATED_FILE );
#
#       chomp( $fileContent[1] );
#       print THREECLASS_FILE "${fileContent[1]}\n";
#
#       close( THREECLASS_FILE );
#
#      $image = $directories . $filename . "/" . $filename . ".nii.gz";
#      $segmentation = $directories . $filename . "/segmentation.nii.gz";
#
#      @args = ( "snap", "-g", $image, "-s", $segmentation );
#      system( @args ) == 0 or die "system @args failed: $?";
#      print "@args\n";
#       }
    }
  }


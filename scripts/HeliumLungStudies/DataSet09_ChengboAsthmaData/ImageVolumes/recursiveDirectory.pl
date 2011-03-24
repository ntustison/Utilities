#!/user/bin/perl -w

use File::Find;
use File::Basename;

$directory = $ARGV[0];
find( \&wanted, $directory );

sub wanted
  {
  ( $filename, $directories, $suffix ) = fileparse( $File::Find::name );
  @current_vars = split( /\//, $directories );

  ##
  # do ct alone
  ##

  $ct_m[0] = "stats_bothlungs_run00.txt";
  $ct_m[1] = "stats_bothlungs_run00.txt";
  $ct_m[2] = "rlm_bothlungs_run00.txt";
  $ct_m[3] = "rlm_bothlungs_run00.txt";
  $ct_m[4] = "stats_bothlungs_run00.txt";

  $ct_w[0] = 6;
  $ct_w[1] = 4;
  $ct_w[2] = 5;
  $ct_w[3] = 3;
  $ct_w[4] = 8;

  # inspiration = SR3, SR4
  # expiration = SR6, SR7

  $ct_p[0] = "expiration";
  $ct_p[1] = "inspiration";
  $ct_p[2] = "inspiration";
  $ct_p[3] = "expiration";
  $ct_p[4] = "inspiration";

  ##
  # do helium-3 alone
  ##

  $he3_m[0] = "otsu_ventilation_ratio_1_both_run02.txt";
  $he3_m[1] = "otsu_ventilation_ratio_1_both_run02.txt";
  $he3_m[2] = "rlm_bothlungs_run02.txt";
  $he3_m[3] = "rlm_bothlungs_run02.txt";
  $he3_m[4] = "cm_bothlungs_run02.txt";

  $he3_w[0] = 0;
  $he3_w[1] = 0;
  $he3_w[2] = 4;
  $he3_w[3] = 7;
  $he3_w[4] = 1;

  $he3_p[0] = "pre";
  $he3_p[1] = "post";
  $he3_p[2] = "post";
  $he3_p[3] = "post";
  $he3_p[4] = "post";

  ##
  # do helium-3 and ct registration
  ##

  $ants_m[0] = "statsFromImagePair_LungPre1_2.txt";
  $ants_m[1] = "statsFromImagePair_LungPost1_2.txt";
  $ants_m[2] = "statsFromImagePair_LungPre2_2.txt";
  $ants_m[3] = "statsFromImagePair_LungPost2_2.txt";



  open CT1FILE, ">>/home/tustison/Data/HeliumLungStudies/DataSet09_Chengbo/ImageVolumes/CT1measurements.txt" or die $!;
  open CT2FILE, ">>/home/tustison/Data/HeliumLungStudies/DataSet09_Chengbo/ImageVolumes/CT2measurements.txt" or die $!;

  open SET1FILE, ">>/home/tustison/Data/HeliumLungStudies/DataSet09_Chengbo/ImageVolumes/He3Set1measurements.txt" or die $!;
  open SET2FILE, ">>/home/tustison/Data/HeliumLungStudies/DataSet09_Chengbo/ImageVolumes/He3Set2measurements.txt" or die $!;

  open CORR_SR3_PRE1, ">>/home/tustison/Data/HeliumLungStudies/DataSet09_Chengbo/ImageVolumes/CorrSR3Pre1measurements.txt" or die $!;
  open CORR_SR3_PRE2, ">>/home/tustison/Data/HeliumLungStudies/DataSet09_Chengbo/ImageVolumes/CorrSR3Pre2measurements.txt" or die $!;
  open CORR_SR3_POST1, ">>/home/tustison/Data/HeliumLungStudies/DataSet09_Chengbo/ImageVolumes/CorrSR3Post1measurements.txt" or die $!;
  open CORR_SR3_POST2, ">>/home/tustison/Data/HeliumLungStudies/DataSet09_Chengbo/ImageVolumes/CorrSR3Post2measurements.txt" or die $!;

  open CORR_SR4_PRE1, ">>/home/tustison/Data/HeliumLungStudies/DataSet09_Chengbo/ImageVolumes/CorrSR4Pre1measurements.txt" or die $!;
  open CORR_SR4_PRE2, ">>/home/tustison/Data/HeliumLungStudies/DataSet09_Chengbo/ImageVolumes/CorrSR4Pre2measurements.txt" or die $!;
  open CORR_SR4_POST1, ">>/home/tustison/Data/HeliumLungStudies/DataSet09_Chengbo/ImageVolumes/CorrSR4Post1measurements.txt" or die $!;
  open CORR_SR4_POST2, ">>/home/tustison/Data/HeliumLungStudies/DataSet09_Chengbo/ImageVolumes/CorrSR4Post2measurements.txt" or die $!;

  if( -d $filename && $filename =~ m/^Set/ )
    {

    print CT1FILE "${filename},";
    print CT2FILE "${filename},";
    print SET1FILE "${filename},";
    print SET2FILE "${filename},";

    for( $i = 0; $i < @he3_p; $i++ )
      {
      my @measurements1;
      my @measurements2;
      if( $he3_p[$i] =~ m/pre/ )
        {
        $file = "${directories}/${filename}/PreSteroid/Pre1/${he3_m[$i]}";
        if( -e $file )
          {
          open FILE, $file;
          @lines = <FILE>;
          chomp( @lines );
          @measurements1 = split( / /, $lines[1] );
          close( FILE );
          print SET1FILE $measurements1[$he3_w[$i]] . ",";
          }
        else
          {
          print SET1FILE "0,";
          }

        $file = "${directories}/${filename}/PreSteroid/Pre2/${he3_m[$i]}";
        if( -e $file )
          {
          open FILE, $file;
          @lines = <FILE>;
          chomp( @lines );
          @measurements2 = split( / /, $lines[1] );
          close( FILE );
          print SET2FILE $measurements2[$he3_w[$i]] . ",";
          }
        else
          {
          print SET2FILE "0,";
          }
        }
      else
        {
        $file = "${directories}/${filename}/PostSteroid/Post1/${he3_m[$i]}";
        if( -e $file )
          {
          open FILE, $file;
          @lines = <FILE>;
          chomp( @lines );
          @measurements1 = split( / /, $lines[1] );
          close( FILE );
          print SET1FILE $measurements1[$he3_w[$i]] . ",";
          }
        else
          {
          print SET1FILE "0,";
          }

        $file = "${directories}/${filename}/PostSteroid/Post2/${he3_m[$i]}";
        if( -e $file )
          {
          open FILE, $file;
          @lines = <FILE>;
          chomp( @lines );
          @measurements2 = split( / /, $lines[1] );
          close( FILE );
          print SET2FILE $measurements2[$he3_w[$i]] . ",";
          }
        else
          {
          print SET2FILE "0,";
          }
        }
      }
    print SET1FILE "\n";
    print SET2FILE "\n";


    for( $i = 0; $i < @ct_p; $i++ )
      {
      my @measurements1;
      my @measurements2;
      if( $ct_p[$i] =~ m/inspiration/ )
        {
        $file = "${directories}/${filename}/CT/SR3/${ct_m[$i]}";
        if( -e $file )
          {
          open FILE, $file;
          @lines = <FILE>;
          chomp( @lines );
          @measurements1 = split( / /, $lines[1] );
          close( FILE );
          print CT1FILE $measurements1[$ct_w[$i]] . ",";
          }
        else
          {
#          print $file . "\n";
          print CT1FILE "0,";
          }

        $file = "${directories}/${filename}/CT/SR4/${ct_m[$i]}";
        if( -e $file )
          {
          open FILE, $file;
          @lines = <FILE>;
          chomp( @lines );
          @measurements2 = split( / /, $lines[1] );
          close( FILE );
          print CT2FILE $measurements2[$ct_w[$i]] . ",";
          }
        else
          {
#          print $file . "\n";
          print CT2FILE "0,";
          }
        }
      else
        {
        $file = "${directories}/${filename}/CT/SR6/${ct_m[$i]}";
        if( -e $file )
          {
          open FILE, $file;
          @lines = <FILE>;
          chomp( @lines );
          @measurements1 = split( / /, $lines[1] );
          close( FILE );
          print CT1FILE $measurements1[$ct_w[$i]] . ",";
          }
        else
          {
#          print $file . "\n";
          print CT1FILE "0,";
          }

        $file = "${directories}/${filename}/CT/SR7/${ct_m[$i]}";
        if( -e $file )
          {
          open FILE, $file;
          @lines = <FILE>;
          chomp( @lines );
          @measurements2 = split( / /, $lines[1] );
          close( FILE );
          print CT2FILE $measurements2[$ct_w[$i]] . ",";
          }
        else
          {
#          print $file . "\n";
          print CT2FILE "0,";
          }
        }
      }
    print CT1FILE "\n";
    print CT2FILE "\n";
    }
  elsif( ( ( $filename =~ m/^SR3$/ || $filename =~ m/^SR4$/ ) &&
    !( $current_vars[@current_vars-1] =~ m/RegistrationResults/ ) ) )
    {
#     print $filename . "\n";
#     print $directories . "\n";
#     print @current_vars . "\n";
#     print "    " . $current_vars[@current_vars-1] . "\n";

    if( $filename =~ m/SR3/ )
      {
      print CORR_SR3_PRE1 "${current_vars[@current_vars-2]} ";
      print CORR_SR3_PRE2 "${current_vars[@current_vars-2]} ";
      print CORR_SR3_POST1 "${current_vars[@current_vars-2]} ";
      print CORR_SR3_POST2 "${current_vars[@current_vars-2]} ";
      }
    elsif( $filename =~ m/SR4/ )
      {
      print CORR_SR4_PRE1 "${current_vars[@current_vars-2]} ";
      print CORR_SR4_PRE2 "${current_vars[@current_vars-2]} ";
      print CORR_SR4_POST1 "${current_vars[@current_vars-2]} ";
      print CORR_SR4_POST2 "${current_vars[@current_vars-2]} ";
      }

    for( $i = 0; $i < @ants_m; $i++ )
      {
      my @measurements = (0,0,0,0,0,0,0,0,0);

      $file = "${directories}/${filename}/${ants_m[$i]}";
      if( -e $file )
        {
        open FILE, $file;
        @lines = <FILE>;
        chomp( @lines );
        @measurements = split( / /, $lines[1] );
        close( FILE );
        }

      if( $filename =~ m/SR3/ )
        {
        if( $ants_m[$i] =~ m/Pre1/ )
          {
          print CORR_SR3_PRE1 "@measurements" . "\n";
          }
        if( $ants_m[$i] =~ m/Pre2/ )
          {
          print CORR_SR3_PRE2 "@measurements" . "\n";
          }
        if( $ants_m[$i] =~ m/Post1/ )
          {
          print CORR_SR3_POST1 "@measurements" . "\n";
          }
        if( $ants_m[$i] =~ m/Post2/ )
          {
          print CORR_SR3_POST2 "@measurements" . "\n";
          }
        }
      else
        {
        if( $ants_m[$i] =~ m/Pre1/ )
          {
          print CORR_SR4_PRE1 "@measurements" . "\n";
          }
        if( $ants_m[$i] =~ m/Pre2/ )
          {
          print CORR_SR4_PRE2 "@measurements" . "\n";
          }
        if( $ants_m[$i] =~ m/Post1/ )
          {
          print CORR_SR4_POST1 "@measurements" . "\n";
          }
        if( $ants_m[$i] =~ m/Post2/ )
          {
          print CORR_SR4_POST2 "@measurements" . "\n";
          }
        }
      }
    }

  close( CT1FILE );
  close( CT2FILE );
  close( SET1FILE );
  close( SET2FILE );

  close( CORR_SR3_PRE1 );
  close( CORR_SR3_PRE2 );
  close( CORR_SR3_POST1 );
  close( CORR_SR3_POST2 );

  close( CORR_SR4_PRE1 );
  close( CORR_SR4_PRE2 );
  close( CORR_SR4_POST1 );
  close( CORR_SR4_POST2 );

  }

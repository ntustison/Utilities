#! /usr/bin/perl -w

use File::Find;
use File::Basename;
use File::Path;
use File::Spec;

my ( $inputDataDir, $outputDataDir, $outputScriptDir ) = @ARGV;

$inputDataDir = File::Spec->rel2abs( $inputDataDir );
$outputDataDir = File::Spec->rel2abs( $outputDataDir );
$outputScriptDir = File::Spec->rel2abs( $outputScriptDir );
if( ! -e $outputScriptDir )
  {
  mkpath( $outputScriptDir, {verbose => 0, mode => 0755} ) or
    die "Can't create output directory $outputScriptDir\n\t";
  }


my $registrationTemplate = "/home/njt4n/share/Data/Public/fsl/MNI152_t1_1mm.nii.gz";

find( \&wanted, $inputDataDir );

sub wanted
  {
  my ( $filename, $directories, $suffix ) = fileparse( $File::Find::name );

  my @comps = split( '/', $directories );

  if( $filename =~ m/WarpedToTemplate\.nii\.gz$/ )
    {

    my $baseScript = "/home/njt4n/share/Pkg/Utilities/scripts/pbs_test.sh";
    open( BASEFILE, "<$baseScript" );
    my @baseScriptContents = <BASEFILE>;
    close( BASEFILE );

    my $subjectname = $filename;
    $subjectname =~ s/WarpedToTemplate\.nii\.gz//;

    my $localOutputDir = "${outputDataDir}/${subjectname}/";
    if( !-e $localOutputDir )
      {
      mkpath( $localOutputDir, {verbose => 0, mode => 0755} ) or
        die "Can't create output directory $localOutputDir\n\t";
      }

    $baseScriptContents[9] = "T1=${directories}/${filename}\n";
    $baseScriptContents[10] = "REGISTRATION_TEMPLATE=${registrationTemplate}\n";
    $baseScriptContents[11] = "OUTPUT_DATA_DIR=${localOutputDir}\n";

    my $pbsFile = "${outputScriptDir}/${filename}";
    $pbsFile =~ s/\.nii\.gz$/_pbs\.sh/;
    open( FILE, ">$pbsFile" ) or die "Couldn't open $pbsFile";
    print FILE @baseScriptContents;
    close( FILE );

    system( "qsub $pbsFile" );

    }
  }

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
  Usage: convertFLIRTtoITK.pl <FLIRTfile> <ITKfile>

      *******************            Overview             *********************

      *******************    Command Line Parameters      *********************

      *******************            References           *********************

  };

if (!($#ARGV + 1))
  {
  print "$usage\n";
  exit 0;
  }
# elsif( $#ARGV > 1 )
#  {
#  die "ERROR: Missing arguments, run without args to see usage\n\t";
#  }

my ( $flirt, $itk ) = @ARGV;

# Convert I/O directories to absolute paths
# $flirt = File::Spec->rel2abs( $flirt );
# $itk = File::Spec->rel2abs( $itk );

open( FLIRT, "${flirt}" ) || die "FLIRT file does not exist.";
my @contents = <FLIRT>;
close( FLIRT );



open( ITK, ">${itk}" );

print( ITK "# Insight Transform File V1.0\n" );
print( ITK "# Transform 0\n" );
print( ITK "Transform: MatrixOffsetTransformBase_double_3_3\n" );
print( ITK "Parameters:  " );
for( my $i = 1; $i < 4; $i++ )
  {
  $contents[$i] =~ s/\n|\r//g;
  my @row = split( ' ', $contents[$i] );
  for( my $j = 0; $j < 3; $j++ )
    {
    print ITK $row[$j] . ' ';
    }
  }
for( my $i = 1; $i < 4; $i++ )
  {
  $contents[$i] =~ s/\n|\r//g;
  my @row = split( ' ', $contents[$i] );
  print ITK $row[3] . ' ';
  }
print( ITK "Parameters:  " );

print( ITK "FixedParameters: 0 0 0" );
close( ITK );

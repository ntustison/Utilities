#! /usr/bin/perl -w

use File::Find;
use File::Basename;
use File::Path;
use File::Spec;

my ( $first, $last ) = @ARGV;

for( my $i = 0; $i <= $last; $i++ )
  {
  system( "qdel ${i}.lc5" );
  }

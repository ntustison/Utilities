#!/usr/bin/perl -w

use strict;

use Cwd 'realpath';
use File::Find;
use File::Basename;
use File::Path;
use File::Spec;
use FindBin qw($Bin);


my $inputBaseDir = $ARGV[0];
$inputBaseDir = File::Spec->rel2abs( $inputBaseDir );

find( \&segmentHead, $inputBaseDir );

sub segmentHead
  {
  ( my $filename, my $directories, my $suffix ) = fileparse( $File::Find::name );

  my @comps = split( '/', $directories );

  $directories =~ s/\(/\\\(/g;
  $directories =~ s/\)/\\\)/g;
  $directories =~ s/\ /\\\ /g;

  if( $filename =~ m/^MR/ && $filename !~ /dcm$/ )
    {
    my $file = "${directories}/${filename}";
    my $dcmfile = "${file}.dcm";

    print "$file -> $dcmfile\n";

    `mv $file $dcmfile`;
    }
  }

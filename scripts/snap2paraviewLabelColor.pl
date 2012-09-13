#!/usr/bin/perl -w

# Usage:
#
#  snap2paraviewLabelColor.pl snap_labels.txt my_snap_labels > paraview_labels.xml
#
# Code:


use strict;

my $snapFile = $ARGV[0];
my $name = $ARGV[1];

# Assume snap labels in format 0     0    0    0        0  0  0    "Label Name"

open SNAP, "<$snapFile";

print "<ColorMap name=\"$name\" space=\"RGB\">\n";

LINE: while (my $line = <SNAP>) {

 next LINE unless $line =~ m/\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+\d(.\d+)?\s+\d\s+\d\s+".*"/;

 my $label = $1;

 # Labels will be off in Paraview unless you rescale, or omit label 0
 next if ($label == 0);

 my $red = $2 / 255.0;
 my $green = $3 / 255.0;
 my $blue = $4 / 255.0;

 print "<Point x=\"${label}.000\" o=\"1\" r=\"${red}.0\" g=\"${green}.0\" b=\"${blue}.0\"/>\n";

}

print "</ColorMap>\n";

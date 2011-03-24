#!/bin/bash

cputime(){
a=$(/usr/bin/time -f "\ncputimesonggang  %U+%S" $*  2>&1 | grep cputimesonggang | awk '{print $2}' | bc)
echo cpu time elapsed: $a seconds
}


cputime $*



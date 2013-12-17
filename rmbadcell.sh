#!/bin/sh

file=$1
fileout=$2
grep -v " 0 0 0" $file \
 | grep -v " 1 0 0" \
 | grep -v " 1 0.1 0" \
 | grep -v " 1.0 0.1 0 0"  \
&> $fileout

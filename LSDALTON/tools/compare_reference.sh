#!/bin/sh
# Shell script for comparing a new profiling summary file with a reference summary file
# Written by Simen Reine January 13th 2012
#
# Syntax:
#
#      sh compare_reference.sh <new summary file> <reference summary file> <cpu-fraction tolerace> <wall-fraction tolerance> <memory tolerance>
#
# where for example a <cpu-fraction tolerace> of 20 indicate that we will not allow an increase of CPU time by more than a maximum of 20 %
#
if [ $# -ne 5 ]
  then
  echo "Usage $0 \<new-sumary file\> \<ref-sumary file\> \<cpu-fraction tolerace\> \<wall-fraction tolerance\> \<memory tolerance\>" 
  exit 2
fi

new=$1
ref=$2
tolCPU=`echo "scale=20; 1 + $3/100" | bc`
tolWall=`echo "scale=20; 1 + $4/100" | bc`
tolMem=`echo "scale=20; 1 + $5/100" | bc`

echo "\n"
echo "Comparing profiling summary files $1 and $2 \n\n"
echo "The CPU-time fraction specified to be below:     $tolCPU"
echo "The wall-time fraction below specified to be:    $tolWall"
echo "The memory usage fraction below specified to be: $tolMem"
sh ./compare_summary.sh $1 $2 | awk 'BEGIN {Ident="T";CPU="OK";Wall="OK";Mem="OK";Pass="OK";print "\n","Profiles not included in the comparison:\n"} {if ( $1 != "The") {if ( $2 >= '$tolCPU' ) CPU="NOT ok"; if ( $3 >= '$tolWall' ) Wall="NOT ok"; if ( $4 >= '$tolMem' ) Mem="NOT ok"; if ( $5 == "different" ) Pass="NOT ok"} else {print $0; Ident="F"}} END {{if ( Ident == "T" ) print "none \n"; else print "";} print "Summary of comparison:","\n\n CPU time:                ",CPU,"\n Wall time:               ",Wall,"\n Memory usage:            ",Mem,"\n Passing profiling cases: ",Pass,"\n\n"}'


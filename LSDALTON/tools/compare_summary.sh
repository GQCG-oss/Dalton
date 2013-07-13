#!/bin/sh
# Shell script for comparing two profiling summary files
# Written by Simen Reine January 13th 2012
#
# Syntax:
#
#      sh compare_summary.sh <new summary file> <reference summary file>
#
if [ $# -ne 2 ]
  then
  echo "Usage $0 \<file1\> \<file2\>"
  exit 2
fi


new=$1
ref=$2
nnew=`grep -ic '^' $new`
nref=`grep -ic '^' $ref`
basenew=`echo "$new" | cut -d'.' -f1`
baseref=`echo "$ref" | cut -d'.' -f1`

# Find number of items per line in both new and reference summary
maxnumnew=0
for (( i=1; $i <= $nnew; i++ ))
do
  newline=`head -${i} $new | tail -1`
  numitems=0
  for item in ${newline}
  do
    numitems=`expr $numitems \+ 1`
  done
  if [ $numitems -gt $maxnumnew ]; then
    maxnumnew=$numitems
  fi
done
maxnumref=0
for (( i=1; $i <= $nref; i++ ))
do
  newline=`head -${i} $ref | tail -1`
  numitems=0
  for item in ${newline}
  do
    numitems=`expr $numitems \+ 1`
  done
  if [ $numitems -gt $maxnumref ]; then
    maxnumref=$numitems
  fi
done

if [ $maxnumnew -ne $maxnumref ]; then
  echo "Mismatching number of items between new summary (${maxnumnew}) and reference (${maxnumref}). Something is wrong!"
  exit 1
fi

nitems=$maxnumnew

# Extract new results
for (( i=1; $i <= $nnew; i++ ))
do
  newline=`head -${i} $new | tail -1`
  numitems=0
  for item in ${newline}
  do
    numitems=`expr $numitems \+ 1`
    let "index = ($i - 1) * $nitems + $numitems"
    newresults[$index]=$item
  done
  if [ $numitems -ne $nitems ]; then
    echo "Wrong number of items (${numitems} instead of ${nitems}) for new results line number ${i}: ${newline}"
    exit 1
  fi
done

# Extract ref results
for (( i=1; $i <= $nref; i++ ))
do
  newline=`head -${i} $ref | tail -1`
  numitems=0
  for item in ${newline}
  do
    numitems=`expr $numitems \+ 1`
    let "index = ($i - 1) * $nitems + $numitems"
    refresults[$index]=$item
  done
  if [ $numitems -ne $nitems ]; then
    echo "Wrong number of items (${numitems} instead of ${nitems}) for ref results line number ${i}: ${newline}"
    exit 1
  fi
done

for (( i=1; $i <= $nnew; i++ ))
do
  let "indI = ($i - 1) * $nitems + 1"
  found[$i]=0
  for (( j=1; $j <= $nref; j++ ))
  do
    let "indJ = ($j - 1) * $nitems + 1"
    if [ ${newresults[$indI]} = ${refresults[$indJ]} ]; then
      found[$i]=1
      printresult[1]=${newresults[$indI]}
      for (( k=2; $k <= $nitems; k++ ))
      do
        let "indI = $indI + 1"
        let "indJ = $indJ + 1"
#       Only the two options logical or number can be compared
        if [ ${newresults[$indI]} = "T" ] || [ ${newresults[$indI]} = "F" ];then
          if [ ${newresults[$indI]} = ${refresults[$indJ]} ];then
            fraction="same"
          else
            fraction="different"
          fi
        else
          fraction=`echo "scale=20; ${newresults[$indI]}/${refresults[$indJ]}" | bc`
        fi
        printresult[$k]=$fraction
      done
      echo "${printresult[*]}"
#      echo "${newresults[$indI]} ${refresults[$indJ]} found"
    fi
  done
done

# Print new profiling cases not in reference
for (( i=1; $i <= $nnew; i++ ))
do
  if [ ${found[$i]} -eq 0 ]; then
    let "indI = ($i - 1) * $nitems + 1"
    echo "The profiling case ${newresults[$indI]} only in the new summary file"
  fi
done

# Print missing profiling cases (only in reference not in new summary)
for (( j=1; $j <= $nref; j++ ))
do
  let "indJ = ($j - 1) * $nitems + 1"
  found[$j]=0
  for (( i=1; $i <= $nnew; i++ ))
  do
    let "indI = ($i - 1) * $nitems + 1"
    if [ ${newresults[$indI]} = ${refresults[$indJ]} ]; then
      found[$j]=1
    fi
  done
  if [ ${found[$j]} -eq 0 ]; then
    let "indJ = ($j - 1) * $nitems + 1"
    echo "The profiling case ${refresults[$indJ]} only in the reference summary file"
  fi
done


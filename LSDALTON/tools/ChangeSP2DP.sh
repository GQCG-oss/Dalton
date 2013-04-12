#!/bin/bash 

dirlist=`ls -d */`
test=2
test2=2
echo $dirlist
for k in $dirlist
do
echo $k
if [ $k != "pdpack/" ]
then
    cd $k
    list=`ls *.f90`
    echo $list
    for i in $list
    do
	echo $i
	tail -n +937 $i > tmp
	mv tmp $i
    done
    list=`ls *.F`
    echo $list
    for i in $list
    do
	echo $i
	tail -n +937 $i > tmp
	mv tmp $i
    done
    list=`ls *.c`
    echo $list
    for i in $list
    do
	echo $i
	tail -n +937 $i > tmp
	mv tmp $i
    done
    cd ..
fi
done



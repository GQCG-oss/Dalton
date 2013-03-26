#!/bin/bash 

dirlist=`ls -d */`
test=2
test2=2
#echo $dirlist
for k in $dirlist
do
if [ $k != "pdpack/" ]
then
    cd $k
    list=`ls *.f90`
#echo $list
    for i in $list
    do
	cat /home/tkjaer/daltondev/SinglePrecision/tools/single.h > tmp
	cat $i >> tmp
	mv tmp $i
    done
    list=`ls *.F`
#echo $list
    for i in $list
    do
	cat /home/tkjaer/daltondev/SinglePrecision/tools/single.h > tmp
	cat $i >> tmp
	mv tmp $i
    done
    list=`ls *.c`
#echo $list
    for i in $list
    do
	cat /home/tkjaer/daltondev/SinglePrecision/tools/single.h > tmp
	cat $i >> tmp
	mv tmp $i
    done
    cd ..
fi
done



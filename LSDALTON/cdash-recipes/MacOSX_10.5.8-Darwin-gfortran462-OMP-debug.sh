#!/bin/bash
source /Users/patrime/.bash_profile
WRK=/Users/patrime/NightlyTestingCdash
CC=/sw/bin/gcc-4
FC=/sw/bin/gfortran
TMP_DIR=/Users/patrime/NightlyTestingCdash/tmp/Dalton-gfortran-OMP-debug
SCRATCH_DIR=$TMP_DIR/scratch
BUILD_NAME=build-debug-omp
BUILDNAME="MacOSX_10.5.8-Darwin-gfortran462-OMP-debug"

cd $WRK
export OMP_NUM_THREADS=2

## Mac OS X does not include /usr/local/bin in its default PATH
export PATH=/usr/local/bin:$PATH
export PATH=/sw/bin/:$PATH
export BLAS_LIB="-framework vecLib"
MKL_NUM_THREADS=1
export MKL_NUM_THREADS

git config --global user.name "Patrick Merlot"
git config --global user.email "patrick.merlot@gmail.com"


rm -rf $TMP_DIR 
rm -rf $SCRATCH_DIR

/sw/bin/git clone git@repo.ctcc.no:dalton.git $TMP_DIR
cd $TMP_DIR
/sw/bin/git checkout -b linsca-develop  origin/linsca-develop

./setup --fc=$FC --cc=$CC --debug --omp  --build=$BUILD_NAME --scratch=$SCRATCH_DIR -D BUILDNAME=$BUILDNAME
cd $TMP_DIR/$BUILD_NAME;Pwd;

#make -j 2;make test VERBOSE=1
#make Experimental
make Nightly

rm -rf $TMP_DIR 
rm -rf $SCRATCH_DIR


exit 0
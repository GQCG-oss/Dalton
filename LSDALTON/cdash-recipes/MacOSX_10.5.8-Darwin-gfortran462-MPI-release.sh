#!/bin/bash
source /Users/patrime/.bash_profile
WRK=/Users/patrime/NightlyTestingCdash
CC=/opt/openmpi/bin/mpicc
FC=/opt/openmpi/bin/mpif90
TMP_DIR=/Users/patrime/NightlyTestingCdash/tmp/Dalton-gfortran-MPI-release
SCRATCH_DIR=$TMP_DIR/scratch
BUILD_NAME=build-release-mpi
BUILDNAME="MacOSX_10.5.8-Darwin-gfortran462-MPI-release"

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

# running testcases
rm -rf $TMP_DIR 
rm -rf $SCRATCH_DIR

/sw/bin/git clone git@repo.ctcc.no:dalton.git $TMP_DIR
cd $TMP_DIR
/sw/bin/git checkout -b linsca-develop  origin/linsca-develop

./setup --fc=$FC --cc=$CC --mpi  --build=$BUILD_NAME --scratch=$SCRATCH_DIR -D BUILDNAME=$BUILDNAME
cd $TMP_DIR/$BUILD_NAME;Pwd;

#make -j 2;make test
#make Experimental
make Nightly


exit 0

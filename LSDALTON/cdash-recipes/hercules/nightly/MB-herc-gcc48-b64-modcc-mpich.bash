bname=`basename "$0"`
##########

source /opt/mpich-3.1.4-gnu-4.8/source.bash

export DALTON_NUM_MPI_PROCS=3
export OMP_NUM_THREADS=2
export CTEST_PROJECT_NAME=LSDALTON
export DALTON_TMPDIR=$wrk/$lib/$bname

root=$(pwd)
wrk=$1
lib=lsdalton_nightly
if [ ! -d $wrk/$lib ]
then
   cd $wrk
   git clone --recursive git@gitlab.com:dalton/dalton.git $lib
fi
cd $wrk/$lib
#
if [ ! -d $wrk/$lib/$bname ]
then
   ./setup --fc=mpif90 --cc=mpicc --cxx=mpic++ --mpi --omp --type=debug --coverage -DENABLE_DEC=ON -DENABLE_TENSORS=ON -DENABLE_RSP=OFF -DENABLE_XCFUN=OFF -DBUILDNAME="$bname" $bname
fi
if [ ! -d $DALTON_TMPDIR ]
then
   mkdir $DALTON_TMPDIR
fi
#
cd $bname
#
ctest -D Nightly
#
cd $root
rm -rf $wrk/$lib

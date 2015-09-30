bname=$(basename "$0")
bname=${bname/.bash/}
##########
#
source /opt/mpich-3.1.4-gnu-4.8/source.bash
#
export DALTON_NUM_MPI_PROCS=3
export OMP_NUM_THREADS=2
export CTEST_PROJECT_NAME=LSDALTON
#
wrk=$1
lib=lsdalton_$bname
export DALTON_TMPDIR=$wrk/$lib/$bname
if [ ! -d $wrk/$lib ]
then
   cd $wrk
   git clone --recursive git@gitlab.com:dalton/dalton.git $lib
fi
cd $wrk/$lib
#
if [ ! -d $wrk/$lib/$bname ]
then
   ./setup --fc=mpif90 --cc=mpicc --cxx=mpic++ --mpi --omp --type=release --coverage -DENABLE_DEC=ON -DENABLE_TENSORS=ON -DENABLE_RSP=OFF -DENABLE_XCFUN=OFF -DBUILDNAME="$bname" $bname
fi
if [ ! -d $DALTON_TMPDIR ]
then
   mkdir $DALTON_TMPDIR
fi
#
cd $bname
#
ctest -D Nightly -L linsca
#

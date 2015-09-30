bname=$(basename "$0")
bname=${bname/.bash/}
##########
export LM_LICENSE_FILE=/opt/pgi/license.dat
source /opt/pgi/linux86-64/15.7/pgi.sh
source /opt/pgi/linux86-64/15.7/mpi.sh
export LD_LIBRARY_PATH=/opt/pgi/15.7/share_objects/lib64:$LD_LIBRARY_PATH
export DALTON_NUM_MPI_PROCS=6
export OMP_NUM_THREADS=2
export CTEST_PROJECT_NAME=LSDALTON
#
wrk=$1
lib=lsdalton_$bname
export DALTON_TMPDIR=$wrk/$lib/$bname/tmp
#
if [ ! -d $wrk/$lib ]
then
   cd $wrk
   git clone --recursive git@gitlab.com:dalton/dalton.git $lib
fi
#
cd $wrk/$lib
#
if [ ! -d $wrk/$lib/$bname ]
then
   ./setup --fc=mpif90 --cc=mpicc --cxx=mpic++ --mpi --omp --int64 --type=release -DENABLE_DEC=ON -DENABLE_TENSORS=ON -DENABLE_RSP=OFF -DENABLE_XCFUN=OFF -DBUILDNAME="$bname" $bname
fi
if [ ! -d $DALTON_TMPDIR ]
then
   mkdir $DALTON_TMPDIR
fi
#
cd $bname
#
ctest -D Nightly -L linsca

#!/bin/bash

BRANCH='master'
COMPILER='GNU'
NUM_CORES=1
MAKE_RELEASE=0
BUILDNAME='undefined'
INT64=0
TRACK='Nightly'

while getopts b:m:r:c:n:i:t: opt; do
  case $opt in
  b)
      BRANCH=$OPTARG
      ;;
  m)
      NUM_CORES=$OPTARG
      ;;
  r)
      MAKE_RELEASE=$OPTARG
      ;;
  c)
      COMPILER=$OPTARG
      ;;
  n)
      BUILDNAME=$OPTARG
      ;;
  i)
      INT64=$OPTARG
      ;;
  t)
      TRACK=$OPTARG
      ;;
  esac
done

shift $((OPTIND - 1))

CTEST_TEMP_DIR=/tmp/cdash-scratch
mkdir -p $CTEST_TEMP_DIR

export DALTON_TMPDIR=$CTEST_TEMP_DIR/run

git clone --recursive git@gitlab.com:dalton/dalton.git $CTEST_TEMP_DIR/compile

cd $CTEST_TEMP_DIR/compile
git checkout $BRANCH &> /dev/null

if [ $MAKE_RELEASE -eq 1 ]; then
    mkdir -p build
    cd build
    cmake ..
    make release
    tar xvzf DALTON-Source.tar.gz
    cd DALTON-Source
fi

if [ "$COMPILER" == 'Intel' ]; then
    source /opt/intel/bin/compilervars.sh intel64
    SETUP_FLAGS="--fc=ifort --cc=icc --cxx=icpc --mkl=sequential"
else
    export MATH_ROOT=/opt/intel/mkl
    SETUP_FLAGS="--fc=gfortran --cc=gcc --cxx=g++"
fi

if [ $INT64 -eq 1 ]; then
    SETUP_FLAGS="$SETUP_FLAGS --int64"
fi

./setup $SETUP_FLAGS -D BUILDNAME=$BUILDNAME

cd build
make -j$NUM_CORES
ctest -D Nightly -j$NUM_CORES --track $TRACK

rm -rf $CTEST_TEMP_DIR

exit 0

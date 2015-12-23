#!/bin/bash

BRANCH='master'
COMPILER='GNU'
MAKE_RELEASE=0
BUILDNAME='undefined'
INT64=0
TRACK='Nightly'
DEPLOY_KEY='undefined'
NUM_JOBS=1

while getopts b:r:c:n:i:t:k:j: opt; do
  case $opt in
  b)
      BRANCH=$OPTARG
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
  k)
      DEPLOY_KEY=$OPTARG
      ;;
  j)
      NUM_JOBS=$OPTARG
      ;;
  esac
done

shift $((OPTIND - 1))

CTEST_TEMP_DIR=/tmp/cdash-scratch
mkdir -p $CTEST_TEMP_DIR

export DALTON_TMPDIR=$CTEST_TEMP_DIR/run

ssh-agent bash -c "ssh-add ${DEPLOY_KEY}; git clone --recursive -b $BRANCH git@gitlab.com:dalton/dalton.git $CTEST_TEMP_DIR/compile"

cd $CTEST_TEMP_DIR/compile

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
make -j$NUM_JOBS
ctest -D Nightly -j$NUM_JOBS --track $TRACK

rm -rf $CTEST_TEMP_DIR

exit 0

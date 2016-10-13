#!/bin/sh

#SBATCH --job-name=prop_exci_aosoppa_trip_h2o_STO2G
#SBATCH -o /home/zth702/source/dalton/DALTON/test/prop_exci_aosoppa_trip/prop_exci_aosoppa_trip_h2o_STO2G-%j.err
#SBATCH -e /home/zth702/source/dalton/DALTON/test/prop_exci_aosoppa_trip/prop_exci_aosoppa_trip_h2o_STO2G-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1"gb"
#SBATCH -t 8000:00:00
#SBATCH -p sauer


date

# set scratch dir
SCRATCH=/scratch/$SLURM_JOB_ID
# make scratch dir
mkdir -p $SCRATCH || exit 0
export DALTON_TMPDIR=$SCRATCH
echo $DALTON_TMPDIR
export DALTON_NUM_MPI_PROCS=1

#
# Define backup function
#
backup () {
  if [ "$1" = "-v" ]; then ekko=1; shift
  else ekko=0;
  fi
  for i in $* ; do
     if [ -f "${i}" ]; then
        for j in 6 5 4 3 2 1 0 ; do
           jp=$(($j+1))
           if [ -f "${i}.${j}" ]; then
              [ $ekko -eq 1 ] && echo "Backup: renaming ${i}.${j} to ${i}.${jp}"
              mv -f "${i}.${j}" "${i}.${jp}"
           fi
        done
        if [ $ekko -eq 1 ]; then echo "Backup: renaming ${i} to ${i}.0"; fi
        mv -f "${i}" "${i}.0"
     fi
  done
}



backup prop_exci_aosoppa_trip_h2o_STO2G.out
/home/zth702/source/dalton/build/dalton -noarch -M 1"000" -t $SCRATCH -d -o prop_exci_aosoppa_trip_h2o_STO2G.out prop_exci_aosoppa_trip h2o_STO2G 

# remove scratch dir
rm -rf $SCRATCH || exit 0

echo ========= Job finished ===================


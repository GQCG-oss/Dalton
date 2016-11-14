#!/usr/bin/env bash
#SBATCH --account=sdujk_fat
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 24
#SBATCH --job-name test
#SBATCH --time=24:00:00
#
# only send emails when the job fails
#SBATCH --mail-type=FAIL
#SBATCH --output dsred_mep.stdout
#SBATCH --error dsred_mep.stderr

module purge
module use /work/sdujk/nhlist/Modules/
module load dalton/pe-efffield
export OMP_NUM_THREADS=1
export DALTON_LAUNCHER="srun"
export DALTON_TMPDIR=$LOCALSCRATCH
srun mkdir -p $DALTON_TMPDIR

dalton -mb 2048 -noappend -t ${LOCALSCRATCH} -o pehf_2pa_lf.out pehf_2pa_lf.dal pehf_2pa_lf pehf_2pa_lf 


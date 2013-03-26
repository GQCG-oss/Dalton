NUM_NODES="1 8"
TASKS_PER_NODE=1
WALL_TIME=8
MAX_MEM=16
PROFS="df-J J XC"
DATE=`date +'%Y_%m_%dT%H%M'`
ACCOUNT="nn4654k"
MAIN_DIR="/usit/titan/u1/simensr/mpi"
BUILD="mpi-build"
PROF_DIR="profile"
MODULES="intel cmake openmpi/1.4.2.intel"
OUTPUTDIR=`pwd`

for NUM_NODE in $NUM_NODES
do
for PROF in $PROFS
do
IDENTIFIER=${PROF}_${NUM_NODE}_${TASKS_PER_NODE}_${DATE}
echo "#!/bin/sh"                                                                      > runjob.prof_${IDENTIFIER}
echo "##"                                                                            >> runjob.prof_${IDENTIFIER}
echo "## Jobname"                                                                    >> runjob.prof_${IDENTIFIER}
echo "#SBATCH --job-name=${IDENTIFIER}"                                              >> runjob.prof_${IDENTIFIER}
echo "##"                                                                            >> runjob.prof_${IDENTIFIER}
echo "## Project"                                                                    >> runjob.prof_${IDENTIFIER}
echo "#SBATCH --account=${ACCOUNT}"                                                  >> runjob.prof_${IDENTIFIER}
echo "##"                                                                            >> runjob.prof_${IDENTIFIER}
echo "## Walltime"                                                                   >> runjob.prof_${IDENTIFIER}
echo "#SBATCH --time=${WALL_TIME}:0:0"                                               >> runjob.prof_${IDENTIFIER}
echo "##"                                                                            >> runjob.prof_${IDENTIFIER}
echo "## Number of nodes"                                                            >> runjob.prof_${IDENTIFIER}
echo "#SBATCH --nodes=${NUM_NODE}"                                                   >> runjob.prof_${IDENTIFIER}
echo "##"                                                                            >> runjob.prof_${IDENTIFIER}
echo "## Number of cores per node"                                                   >> runjob.prof_${IDENTIFIER}
echo "#SBATCH --ntasks-per-node=${TASKS_PER_NODE}"                                   >> runjob.prof_${IDENTIFIER}
echo "##"                                                                            >> runjob.prof_${IDENTIFIER}
echo "## Memory needed"                                                              >> runjob.prof_${IDENTIFIER}
echo "#SBATCH --mem-per-cpu=${MAX_MEM}G"                                             >> runjob.prof_${IDENTIFIER}
echo ""                                                                              >> runjob.prof_${IDENTIFIER}
echo "source /site/bin/jobsetup"                                                     >> runjob.prof_${IDENTIFIER}
echo "export OMP_NUM_THREADS=${TASKS_PER_NODE}"                                      >> runjob.prof_${IDENTIFIER}
echo ""                                                                              >> runjob.prof_${IDENTIFIER}
for MODULE in ${MODULES}
do
echo "module load ${MODULE}"                                                         >> runjob.prof_${IDENTIFIER}
done #MODULE
echo ""                                                                              >> runjob.prof_${IDENTIFIER}
echo "## copy data file to work area"                                                >> runjob.prof_${IDENTIFIER}
echo ""                                                                              >> runjob.prof_${IDENTIFIER}
echo "sbcast ${MAIN_DIR}/${BUILD}/dalton \$SCRATCH/dalton"                           >> runjob.prof_${IDENTIFIER}
echo "sbcast ${MAIN_DIR}/${BUILD}/dalton.x \$SCRATCH/dalton.x"                       >> runjob.prof_${IDENTIFIER}
echo "sbcast ${MAIN_DIR}/${PROF_DIR}/PROF \$SCRATCH/PROF"                            >> runjob.prof_${IDENTIFIER}
echo "cp -r ${MAIN_DIR}/${PROF_DIR}/LSint \$SCRATCH"                                 >> runjob.prof_${IDENTIFIER}
echo "cp -r ${MAIN_DIR}/${PROF_DIR}/DECtest \$SCRATCH"                               >> runjob.prof_${IDENTIFIER}
echo ""                                                                              >> runjob.prof_${IDENTIFIER}
echo "chkfile ${IDENTIFIER}.log"                                                     >> runjob.prof_${IDENTIFIER}
echo "chkfile ${IDENTIFIER}.alllog"                                                  >> runjob.prof_${IDENTIFIER}
echo ""                                                                              >> runjob.prof_${IDENTIFIER}
echo "cd \$SCRATCH"                                                                  >> runjob.prof_${IDENTIFIER}
echo "./PROF -log ${IDENTIFIER}.log -dalton ./dalton -param \"-N ${NUM_NODE} -exe \$SCRATCH/dalton.x -t \$SCRATCH/tmp\" ${PROF} &> ${IDENTIFIER}.allout" >> runjob.prof_${IDENTIFIER}
echo ""                                                                              >> runjob.prof_${IDENTIFIER}
echo "cp ${IDENTIFIER}.log ${OUTPUTDIR}/"                                           >> runjob.prof_${IDENTIFIER}
echo "cp ${IDENTIFIER}.allout ${OUTPUTDIR}/"                                        >> runjob.prof_${IDENTIFIER}
echo ""                                                                              >> runjob.prof_${IDENTIFIER}
echo ""                                                                              >> runjob.prof_${IDENTIFIER}
echo ""                                                                              >> runjob.prof_${IDENTIFIER}
echo ""                                                                              >> runjob.prof_${IDENTIFIER}
echo ""                                                                              >> runjob.prof_${IDENTIFIER}
echo ""                                                                              >> runjob.prof_${IDENTIFIER}
echo "## end of script"                                                              >> runjob.prof_${IDENTIFIER}

sbatch runjob.prof_${IDENTIFIER}
done #PROFS
done #NUM_NODES

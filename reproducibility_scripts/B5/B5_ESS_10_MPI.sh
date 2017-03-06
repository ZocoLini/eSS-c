SACESS_HOME=/home/davidp/bitbucket
INPUT_PATH=inputs/reproducibility/TEMPLATE_B5_eSSm.xml
OUTPUT_PATH=output/B5_ESS_10_MPI
NUM_PROC=10 # NUMBRE OF MPI PARALLEL PROCESSORS
OMP_NUM_THREADS=1 # NUMBER OF OPENMP THREAD PER MPI PROCESSOR
INIT_RUN=1 # ID OF THE FIRST RUN
NUMBER_OF_RUNS=20 # ID OF THE LAST RUN
MPIRUN=mpirun

##################################################################################################################
CURRENT_PWD=$(pwd)
cd $SACESS_HOME
module load intel/2.144
module load impi/5.0.0.028
export OMP_NUM_THREADS=$OMP_NUM_THREADS

for i in `seq $INIT_RUN $NUMBER_OF_RUNS`; do
   echo RUN$i
   if [ "$i" -lt "10" ]; then
	$MPIRUN -np $NUM_PROC bin/paralleltestbed $INPUT_PATH "$OUTPUT_PATH"_0$i
   else
        $MPIRUN -np $NUM_PROC bin/paralleltestbed $INPUT_PATH "$OUTPUT_PATH"_$i
   fi 
done
cd $CURRENT_PWD
##################################################################################################################
echo END OF REPRODUCIBILITY SCRIPT

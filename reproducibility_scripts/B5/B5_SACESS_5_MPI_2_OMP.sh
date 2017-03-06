SACESS_HOME=/home/davidp/bitbucket
INPUT_PATH=inputs/reproducibility/TEMPLATE_B5_saCeSS.xml
OUTPUT_PATH=output/B5_SACESS_6_MPI_2_OMP
NUM_PROC=6 # NUMBRE OF MPI PARALLEL PROCESSORS
OMP_NUM_THREADS=2 # NUMBER OF OPENMP THREAD PER MPI PROCESSOR
INIT_RUN=1 # ID OF THE FIRST RUN
NUMBER_OF_RUNS=20 # ID OF THE LAST RUN
MPIRUN=mpirun

module load intel/2.144
module load impi/5.0.0.028

##################################################################################################################
CURRENT_PWD=$(pwd)
cd $SACESS_HOME
export OMP_NUM_THREADS=$OMP_NUM_THREADS
export I_MPI_PIN_DOMAIN=omp

for i in `seq $INIT_RUN $NUMBER_OF_RUNS`; do
   echo RUN$i
   if [ "$i" -lt "10" ]; then
       $MPIRUN --map-by slot:pe=$OMP_NUM_THREADS -np $NUM_PROC bin/paralleltestbed $INPUT_PATH "$OUTPUT_PATH"_0$i
   else
       $MPIRUN --map-by slot:pe=$OMP_NUM_THREADS -np $NUM_PROC bin/paralleltestbed $INPUT_PATH "$OUTPUT_PATH"_$i
   fi 
done
cd $CURRENT_PWD
##################################################################################################################
echo END OF REPRODUCIBILITY SCRIPT

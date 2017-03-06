SACESS_HOME=/home/david.penas/DavidPenas-sacess-library-b6f08b860add 
INPUT_PATH=inputs/reproducibility/TEMPLATE_B2_CeSS_TIME_800.xml
OUTPUT_PATH=output/B2_CESS_10_MPI_TIME_800
NUM_PROC=11 # NUMBRE OF MPI PARALLEL PROCESSORS
OMP_NUM_THREADS=1 # NUMBER OF OPENMP THREAD PER MPI PROCESSOR
INIT_RUN=1 # ID OF THE FIRST RUN
NUMBER_OF_RUNS=20 # ID OF THE LAST RUN
MPIRUN=/home/david.penas/openmpi/bin/mpirun

##################################################################################################################
CURRENT_PWD=$(pwd)
cd $SACESS_HOME
#module load intel/64/compiler/13.1_up1
export OMP_NUM_THREADS=$OMP_NUM_THREADS

        #!/bin/bash
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

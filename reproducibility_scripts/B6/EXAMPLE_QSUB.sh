# Examples of QSUB commmands to send the reproducibility scripts to the queue in SGE

# MPI RUNs
qsub -cwd -N B6_ESS_10MPI    -l excl=true,s_rt=36:00:00,mem_free=1G -q compute-0-x.q -pe mpi_fu 10 B6_ESS_10_MPI.sh
qsub -cwd -N B6_CESS_10MPI_TIME1000  -l excl=true,s_rt=36:00:00,mem_free=1G -q compute-0-x.q -pe mpi_fu 11 B6_CESS_10_MPI_TIME_1000.sh
qsub -cwd -N B6_CESS_10MPI_TIME2000  -l excl=true,s_rt=36:00:00,mem_free=1G -q compute-0-x.q -pe mpi_fu 11 B6_CESS_10_MPI_TIME_2000.sh
qsub -cwd -N B6_SACESS_10MPI -l excl=true,s_rt=36:00:00,mem_free=1G -q compute-0-x.q -pe mpi_fu 11 B6_SACESS_10_MPI.sh




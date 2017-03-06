# Examples of QSUB commmands to send the reproducibility scripts to the queue in SGE

# MPI RUNs
qsub -cwd -N B1_ESS_10MPI    -l excl=true,s_rt=36:00:00,mem_free=1G -q compute-0-x.q -pe mpi_fu 10 B1_ESS_10_MPI.sh
qsub -cwd -N B1_CESS_10MPI   -l excl=true,s_rt=36:00:00,mem_free=1G -q compute-0-x.q -pe mpi_fu 11 B1_CESS_10_MPI.sh
qsub -cwd -N B1_SACESS_10MPI -l excl=true,s_rt=36:00:00,mem_free=1G -q compute-0-x.q -pe mpi_fu 11 B1_SACESS_10_MPI.sh
qsub -cwd -N B1_ESS_20MPI    -l excl=true,s_rt=36:00:00,mem_free=1G -q compute-0-x.q -pe mpi_fu 20 B1_ESS_20_MPI.sh
qsub -cwd -N B1_SACESS_20MPI -l excl=true,s_rt=36:00:00,mem_free=1G -q compute-0-x.q -pe mpi_fu 21 B1_SACESS_20_MPI.sh
qsub -cwd -N B1_ESS_40MPI    -l excl=true,s_rt=36:00:00,mem_free=1G -q compute-0-x.q -pe mpi_fu 40 B1_ESS_40_MPI.sh
qsub -cwd -N B1_SACESS_40MPI -l excl=true,s_rt=36:00:00,mem_free=1G -q compute-0-x.q -pe mpi_fu 41 B1_SACESS_40_MPI.sh
# hybrid MPI+openMP RUNs
qsub -cwd -N B1_SACESS_5MPI_2_OMP  -l excl=true,s_rt=36:00:00,mem_free=1G -q compute-0-x.q -pe mpi_6p   12 B1_SACESS_5_MPI_2_OMP.sh
qsub -cwd -N B1_SACESS_5MPI_4_OMP  -l excl=true,s_rt=36:00:00,mem_free=1G -q compute-0-x.q -pe mpi_6p   24 B1_SACESS_5_MPI_4_OMP.sh
qsub -cwd -N B1_SACESS_5MPI_8_OMP  -l excl=true,s_rt=36:00:00,mem_free=1G -q compute-0-x.q -pe mpi_12p  48 B1_SACESS_5_MPI_8_OMP.sh
qsub -cwd -N B1_SACESS_10MPI_2_OMP -l excl=true,s_rt=36:00:00,mem_free=1G -q compute-0-x.q -pe mpi_11p  22 B1_SACESS_10_MPI_2_OMP.sh
qsub -cwd -N B1_SACESS_10MPI_4_OMP -l excl=true,s_rt=36:00:00,mem_free=1G -q compute-0-x.q -pe mpi_11p  44 B1_SACESS_10_MPI_4_OMP.sh
qsub -cwd -N B1_SACESS_20MPI_2_OMP -l excl=true,s_rt=36:00:00,mem_free=1G -q compute-0-x.q -pe mpi_6p   42 B1_SACESS_20_MPI_2_OMP.sh


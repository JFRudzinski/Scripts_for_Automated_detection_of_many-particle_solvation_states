#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob.out.%j
#SBATCH -e ./tjob.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J Get_HMM
# Queue (Partition):
#SBATCH --partition=small
# Number of MPI Tasks, e.g. 8:
#SBATCH --ntasks=8
#SBATCH --ntasks-per-core=1
# Memory usage of the job [MB], 3800X2 MB per task:
#SBATCH --mem=15200
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rudzinski@mpip-mainz.mpg.de
#
# Wall clock limit:
#SBATCH --time=24:00:00

# Run the program:
#srun -n 4 /u/jrudz/pkg/miniconda2/bin/mpiexec -n 4 /u/jrudz/pkg/miniconda2/bin/python Get_HMM_AAss0.py &> log.dat
#srun -n 4 /u/jrudz/pkg/miniconda2/bin/python Get_HMM_AAss0.py &> log.dat
#srun -n 4 /u/jrudz/pkg/miniconda2/bin/mpiexec /u/jrudz/pkg/miniconda2/bin/python Get_HMM_AAss0.py &> log.dat
srun -n 8 --mpi=pmi2 /u/jrudz/pkg/miniconda2/bin/python driver_mpi.py &> log.dat



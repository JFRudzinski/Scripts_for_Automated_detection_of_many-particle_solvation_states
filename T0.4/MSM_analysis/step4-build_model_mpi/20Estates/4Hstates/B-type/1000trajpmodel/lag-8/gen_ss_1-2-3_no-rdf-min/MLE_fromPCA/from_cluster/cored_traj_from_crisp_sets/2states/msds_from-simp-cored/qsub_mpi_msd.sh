#! /bin/bash
#This script is generated by rudzinski@pckr177
#on Mon Jan 12 17:02:31 CET 2015
#by the command "q2start -8p --walltime 00:15:00 mdrun_mpi"
#using q2start, version 0.5.5

#sge options
#$ -cwd
#$ -j y
#$ -N Get_mss_props
#$ -e $JOB_ID.log
#$ -o $JOB_ID_2.log
#$ -l h_rt=36:00:00
#$ -pe PE_16 16
#$ -m bes
#$ -M rudzinski@mpip-mainz.mpg.de
#$ -S /bin/bash

echo ""
echo $JOB_ID
echo $HOME
echo $HOSTNAME 
echo $PWD

#module load openmpi

#BEGIN COMPILER and MPI VARIABLES for gnu
#export LD_LIBRARY_PATH="/sw/linux/mpi/gcc/openmpi/lib:$LD_LIBRARY_PATH";
#export PATH="/sw/linux/mpi/gcc/openmpi/bin:$PATH";
#export CC="gcc";
#export F77="gfortran";
#export CXX="g++";
#export MPIEXEC="/sw/linux/mpi/gcc/openmpi/bin/mpirun --prefix /sw/linux/mpi/gcc/openmpi/";
#END  COMPILER and MPI VARIABLES for gnu

#BEGIN useful script variables
#walltime=24:00:00
#wallsecs=900
#wallhours=0
#ncpus=8
#END useful script variables

echo Hi, I am job $JOB_ID on $HOSTNAME in $PWD

echo Starting simulation in $PWD
#MPIEXEC -n 8 mdrun_mpi
# JFR - use my own script
#/sw/linux/mpi/gcc/openmpi/bin/mpiexec 
/home/theorie/rudzinski/soft/anaconda/envs/PyEmma-new/bin/mpiexec -n 16 /home/theorie/rudzinski/soft/anaconda/envs/PyEmma-new/bin/python calc_mss_msds_mpi_cored.py

result=$?
[ $result -ne 0 ] && echo "$JOB_ID finished unhappy!"



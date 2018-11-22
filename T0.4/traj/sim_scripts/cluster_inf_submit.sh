#!/bin/bash
#
#sge options
#$ -cwd
#$ -j y
#$ -N kob-andersen_t0.2
#$ -o $JOB_NAME.out
#$ -e $JOB_NAME.err
#$ -l h_rt=36:00:00
#$ -pe PE_8 8
#$ -m n
#$ -S /bin/bash

# source ~/.bashrc
source /people/thnfs/homes/radumarc/Software/Espressopp_Cluster/espressopp-1.9.2/ESPRC

# WORKDIRECTORY
workdir=`pwd`
echo "Workdir is $workdir"

# WORK SCRATCH
if [ ! -d "/usr/scratch/radumarc" ]
then
  mkdir /usr/scratch/radumarc
fi

jno=0
while [ -d "/usr/scratch/radumarc/job_$jno" ]
do
  jno=`expr $jno + 1`
done

jobdir="/usr/scratch/radumarc/job_$jno"
if mkdir $jobdir
then
  jdir_made=0
else
  jdir_made=-1
fi

if [ $jdir_made -ne -1 ]
then
  rm -rf $jobdir/*
  echo "Jobdir is $jobdir"

  # COPY FOLDER
  rsync -ar $workdir/*.py $jobdir
  if [ -a restart.in.gz ]
  then
    rsync -ar $workdir/restart.in.gz $jobdir
    if [ -a config_series_eq.xyz.gz ]
    then
      if [ ! -a config_series_sim.xyz.gz ]
      then
        rsync -ar $workdir/config_series_eq.xyz.gz $jobdir
      fi
    fi
    if [ -a config_series_sim.xyz.gz ]
    then
      rsync -ar $workdir/config_series_sim.xyz.gz $jobdir
    fi
    if [ -a output_eq.dat ]
    then
      if [ ! -a output_sim.dat ]
      then
        rsync -ar $workdir/output_eq.dat $jobdir
      fi
    fi
    if [ -a output_sim.dat ]
    then
      rsync -ar $workdir/output_sim.dat $jobdir
    fi
  fi
  if [ -a endofsim.txt ]
  then
    rsync -ar $workdir/endofsim.txt $jobdir
  fi

  cd $jobdir
  echo $JOB_ID > MYJOBID.txt
  rsync -a $jobdir/MYJOBID.txt $workdir

  # echo Starting simulation in $PWD
  #/sw/linux/mpi/gcc/openmpi/bin/mpirun --prefix /sw/linux/mpi/gcc/openmpi/ -np 8 python LJ_PARRINELLO.py
  #/sw/linux/mpi/gcc/openmpi/bin/mpirun --prefix /sw/linux/mpi/gcc/openmpi/ -np 8 python LJ_POLYMER.py
  /sw/linux/mpi/gcc/openmpi/bin/mpirun --prefix /sw/linux/mpi/gcc/openmpi/ -np 8 python LJ.py
  #/sw/linux/mpi/gcc/openmpi/bin/mpirun --prefix /sw/linux/mpi/gcc/openmpi/ -np 8 python LJ_CMUMD_POLYMER.py

  if [ $? -eq 0 ]
  then
    # SYNCHRONIZE BACK & CLEAN
    rsync -ar $jobdir/* $workdir --exclude "*.out"

    if [ ! -f endofsim.txt ]
    then
      rm -rf $jobdir
      cd $workdir
      /people/thnfs/homes/radumarc/bin/qsub -pe PE_8 8 -wd $workdir cluster_inf_submit.sh
    else
      rm -rf $jobdir
    fi
  else
    rm -rf $jobdir
    echo "Job was cancelled due to failed simulation!"
    exit 1
  fi

fi

exit 0

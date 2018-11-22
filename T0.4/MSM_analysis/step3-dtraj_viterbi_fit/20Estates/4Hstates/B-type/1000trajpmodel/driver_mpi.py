import os
import numpy as np
# Add the MPI stuff
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
print 'hello from proc '+str(rank)
size = comm.Get_size()

lag = XX
sys_nm = XX

if ( rank == 0 ):

    script_path = './'+sys_nm+'/hmsm_lag-'+str(lag)+'/Dtraj_convert_HMM_'+sys_nm+'.py'
    os.system('/u/jrudz/pkg/miniconda2/bin/python '+script_path)

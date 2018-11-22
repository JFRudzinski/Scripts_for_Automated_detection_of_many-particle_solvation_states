import os
import numpy as np
# Add the MPI stuff
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
print 'hello from proc '+str(rank)
size = comm.Get_size()

lag = np.array( [16, 16, 16, 16, 16, 16, 16, 16] )

sys_nm = 'AAss'+str(rank)
script_path = './'+sys_nm+'/hmsm_lag-'+str(lag[rank])+'/Dtraj_convert_HMM_'+sys_nm+'.py'
os.system('/u/jrudz/pkg/miniconda2/bin/python '+script_path)

lag = np.array( [16, 16, 16, 16, 16, 16, 16, 16] )

sys_nm = 'AB-Ass'+str(rank)
script_path = './'+sys_nm+'/hmsm_lag-'+str(lag[rank])+'/Dtraj_convert_HMM_'+sys_nm+'.py'  
os.system('/u/jrudz/pkg/miniconda2/bin/python '+script_path)

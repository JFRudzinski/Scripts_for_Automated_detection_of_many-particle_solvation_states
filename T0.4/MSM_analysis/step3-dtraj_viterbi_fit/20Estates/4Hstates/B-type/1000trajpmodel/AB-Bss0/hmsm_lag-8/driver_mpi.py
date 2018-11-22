import os
import numpy as np
# Add the MPI stuff
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
print 'hello from proc '+str(rank)
size = comm.Get_size()

lag = 8
sys_nm = 'AB-Bss0'

if ( rank == 0 ):

    script_path = './Dtraj_convert_HMM_'+sys_nm+'.py'
    os.system('/u/jrudz/pkg/miniconda2/bin/python '+script_path)

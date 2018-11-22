import numpy as np
import mdtraj as md
import pyemma
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
from collections import Counter
import mdtraj as md
#from CN_functions import *
#from gen_dtraj_functions import *
import os
from copy import deepcopy
# Add the MPI stuff
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
print 'hello from proc '+str(rank)
size = comm.Get_size()

# some values
Temp = 0.4
nstates = 2

tica_lag = np.load('../../../TICA/tica_lag.npy')
tica_dim = np.load('../../../TICA/tica_dim.npy')
clust_dim = np.load('../../../TICA/tica_dim.npy')
nclust_max = np.load('../../../TICA/nclust_max.npy')
nclust = np.load('../../../TICA/nclust.npy')
Nprune = np.load('../../../TICA/Nprune.npy')
fr_max = 30006

indir = '/data/isilon/rudzinski/cluster_tmp/LJ/Kobb_Andersen/DiffCoeff_red-dt_30ns/T'+str(Temp)
topfile =  indir+'/LJ_AK_5000.pdb'

feat = coor.featurizer(topfile)

ptypes = ['A','B']
part_type = ptypes[1]
other_type = ptypes[0]

self_list = feat.topology.select('name '+part_type)
other_list = feat.topology.select('name '+other_type)

tau_CG = np.load('../tau_CG.npy')
lag = tau_CG

n_mol = feat.topology.n_residues
n_sites_p_mol = feat.topology.n_atoms / n_mol
n_excl = n_sites_p_mol
n_excl

# get the full traj
traj_fnm = indir+'/traj/config_series_sim_T'+str(Temp)+'_nojump_0-30ns.xtc'
# get the dtraj
dtrajs = np.load('../../../TICA/dtrajs_regspaceB_nclust-'+str(nclust)+'_ticadim-'+str(tica_dim)+'_ticalag-'+str(tica_lag)+'_clustdim-'+str(clust_dim)+'.npy')
# get the mss sets
mss_sets = np.load('../mss_sets_'+str(nstates)+'.npy')
Nmss = len(mss_sets)

# Use the lag time of the MSM
tica_lag = lag

# Calculate the distribution of msds for transitions between mss
msd = [[[] for x in range(Nmss)] for y in range(Nmss)]
wait_time = [[] for x in range(Nmss)]

#dtrajs = dtrajs[1:8] # DEBUG
# split the trajs up by process
trajs = np.split(np.arange(size*(len(dtrajs)/size+1)),size)[rank]
trajs = trajs[np.where(trajs<len(dtrajs))[0]]
for traj in trajs:
    # find the first frame in a mss
    fr0 = 0
    flag_first_fr = False
    while ( not flag_first_fr ):
        for mss in range(Nmss):
            if ( len(np.where( mss_sets[mss] == dtrajs[traj][fr0/Nprune] )[0]) != 0 ):
                mss_ind_0 = mss
                flag_first_fr = True
        if ( not flag_first_fr ):
            fr0 += tica_lag*Nprune
    # get some variables ready
    r0 = np.array( md.load_frame(traj_fnm, fr0, top=topfile, atom_indices=[self_list[traj]]).xyz )
    t0 = -1
    shift = 0
    for fr in range(fr0+tica_lag*Nprune,fr_max,tica_lag*Nprune): # dtrajs[traj].shape[0]):
        shift = int(fr/5001) # need to account for the repeat first frame in the split trajectories
        flag_jump = False
        flag_mss = False
        for mss in range(Nmss):
            if ( len(np.where( mss_sets[mss] == dtrajs[traj][fr/Nprune] )[0]) != 0 ):
                mss_ind = mss
                flag_mss = True
        if ( not flag_mss ): # wait for the next fr within some mss state
            continue
        # get the new position
        r = np.array( md.load_frame(traj_fnm, fr-shift, top=topfile, atom_indices=[self_list[traj]]).xyz )
        # get the msd
        msd_inst = np.linalg.norm(r-r0)**2
        msd[mss_ind_0][mss_ind].append(msd_inst)
        # check if there was a jump between mss
        if ( mss_ind_0 != mss_ind ):
            flag_jump = True
        # get the waiting time in the case of a jump
        if ( flag_jump ):
            if (t0 ==-1): # first jump
                t0 = deepcopy(fr)
            else: # store the jump properties
                wait_time[mss_ind_0].append(fr-t0)
                # update the initial time
                t0 = deepcopy(fr)
        # get ready for the next step
        r0 = deepcopy(r)
        mss_ind_0 = deepcopy(mss_ind)
    print 'proc '+str(rank)+' done with traj '+str(traj)

#np.save('msd_mss_trajs_'+str(trajs[0])+'-'+str(trajs[-1]), msd)
#np.save('wait_time_trajs_'+str(trajs[0])+'-'+str(trajs[-1]), wait_time)

# get the data from all the processes
msd_tot = comm.gather(msd, root=0)
wait_time_tot = comm.gather(wait_time, root=0)

if ( rank == 0 ):
    np.save('msd_mss_trajs_all', msd_tot)
    np.save('wait_time_trajs_all', wait_time_tot)



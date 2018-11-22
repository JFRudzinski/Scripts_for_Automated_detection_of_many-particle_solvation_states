
# coding: utf-8

# TICA and clustering with CN- of 80:20 KA sims
# ====

# In[ ]:

import pyemma
pyemma.__version__


# In[ ]:

import os
#get_ipython().magic(u'pylab inline')
#matplotlib.rcParams.update({'font.size': 12})


# In[ ]:

import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt


# Add the MPI stuff
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
print 'hello from proc '+str(rank)
size = comm.Get_size()

# Some useful inhouse functions
# ------

# In[ ]:

#from plot_functions import *
#import viterbi
import numpy as np


import pickle
def save_object(filename, obj):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


# Read in the dtrajs
# ------

# input params
pair_type = 'AB'
cent_type = 'B'
ss = 0
#
n_Estates = 20
dmin = 0.23
n_Hstates = 4
#
Nprune = 1
Nparam_traj = 1000
Ntraj_sets = 1
Ntraj_0 = 0
Ntraj_f = 0
lags = [8]
lag0 = np.zeros(Ntraj_sets).astype(int)
Nign = 7 # number of the processors to leave out, for memory conservation
#
#Nlags_ck_tot = 100
#
sys_nm = pair_type+'-'+cent_type+'ss'+str(ss)+'_'+str(n_Estates)+'Estates'+'_'+str(n_Hstates)+'Hstates'

# get the CN dtrajs
indir_CN = '/data/isilon/rudzinski/cluster_tmp/LJ/Kobb_Andersen/DiffCoeff_red-dt_30ns/T0.4/MSM_analysis/2016_12-Dec_13/step1-dtraj_mpi/'

if ( rank == 0 ):
    dtraj_CNs = np.load(indir_CN+'/dtraj_CNs_'+pair_type+'_0-5ns.npz')
    N_solshel = dtraj_CNs['N_solshel']
    Xre = dtraj_CNs['Xre_'+cent_type]
    dtraj_CNs = np.load(indir_CN+'/dtraj_CNs_'+pair_type+'_5-10ns.npz')
    Xre = np.hstack((Xre,dtraj_CNs['Xre_'+cent_type]))
    dtraj_CNs = np.load(indir_CN+'/dtraj_CNs_'+pair_type+'_10-15ns.npz')
    Xre = np.hstack((Xre,dtraj_CNs['Xre_'+cent_type]))
    dtraj_CNs = np.load(indir_CN+'/dtraj_CNs_'+pair_type+'_15-20ns.npz')
    Xre = np.hstack((Xre,dtraj_CNs['Xre_'+cent_type]))
    dtraj_CNs = np.load(indir_CN+'/dtraj_CNs_'+pair_type+'_20-25ns.npz')
    Xre = np.hstack((Xre,dtraj_CNs['Xre_'+cent_type]))
    dtraj_CNs = np.load(indir_CN+'/dtraj_CNs_'+pair_type+'_25-30ns.npz')
    Xre = np.hstack((Xre,dtraj_CNs['Xre_'+cent_type]))

    # get only the relevant solvation shell
    dtraj_CN = []
    for i in range(len(Xre)):
        dtraj_CN.append(Xre[i][::Nprune,ss]) # prune the data to lighten the load


# calculate the HMSM on subsets of the trajectories
for traj_frac in range(Ntraj_0,Ntraj_f+1):

    if ( rank == 0 ):
        print 'Starting trajfrac '+str(traj_frac)+' of '+str(Ntraj_sets)
        # get the subset
        dtraj_CN_act = dtraj_CN[traj_frac*Nparam_traj/Ntraj_sets:(traj_frac+1)*Nparam_traj/Ntraj_sets]

        # clustering
        n_clusters = n_Estates      # number of clusters
        clustering = coor.cluster_regspace(dtraj_CN_act,max_centers=n_clusters,dmin=dmin)
        save_object('clustering'+sys_nm+'_trajfrac-'+str(traj_frac)+'.pkl', clustering)
        # already did this, read it in
        #with open('clustering'+sys_nm+'_trajfrac-'+str(traj_frac)+'.pkl', 'rb') as f:
        #    clustering = pickle.load(f)
        dtrajs = clustering.dtrajs
        cc = clustering.clustercenters[:,0]
        print 'n_clusters = '+str(len(cc))
    else:
        dtrajs = None

    # send the dtraj info
    dtrajs = comm.bcast(dtrajs,root=0)


    # HMSM
    if ( rank < size-Nign ): # leave out some processors
        nstates = n_Hstates
        for lag in range(lag0[traj_frac],len(lags)):

            print 'Starting lag '+str(lag)+' of '+str(len(lags))
            hmsm = msm.estimate_hidden_markov_model(dtrajs, nstates, lags[lag], reversible=True, stationary=False,stride=1)
            #hmsm = msm.bayesian_hidden_markov_model(dtrajs, nstates, lags[lag], nsamples=Nsamples, reversible=True, stationary=False, stride=1, conf=Iconf)
            save_object('HMSM_'+sys_nm+'_trajfrac-'+str(traj_frac)+'_lag-'+str(lags[lag])+'.pkl', hmsm)
            # already did this, read it in
            #with open('HMSM_'+sys_nm+'_trajfrac-'+str(traj_frac)+'_lag-'+str(lags[lag])+'.pkl', 'rb') as f:
            #    hmsm = pickle.load(f)


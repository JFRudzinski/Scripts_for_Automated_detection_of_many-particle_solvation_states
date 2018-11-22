
# coding: utf-8

# TICA and clustering with CN- of 80:20 KA sims
# ====

# In[1]:

import pyemma
pyemma.__version__


# In[2]:

import os

# In[3]:

import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
from collections import Counter
import mdtraj as md
# from CN_functions import *
# from gen_dtraj_functions import *


# Some useful inhouse functions
# ------

# In[4]:

#from plot_functions import *
import pickle
import numpy as np


# Read in the data
# ------

# In[5]:

N_mss = [2,3,4,5,6,7]
Temp = '0.4'

tica_lag = np.load('../../../PCA/tica_lag.npy')
tica_lag = 1
# tica_dim = np.load('../../../PCA/tica_dim.npy')
# clust_dim = np.load('../../../PCA/clust_dim.npy')
# nclust_max = np.load('../../../PCA/nclust_max.npy')
# nclust = 50 # np.load('../../../TICA/nclust.npy')
Nprune = np.load('../../../PCA/Nprune.npy')
print Nprune
# tica_lag


# In[6]:

#  clustering
Y = []
Y.append(np.load('../../../PCA/Y.npy'))
with open('../../../PCA/clustering_kmeans_nclust-50_clustdim-2.pkl', 'rb') as f:
    clustering = pickle.load(f)


# In[7]:

dtraj = np.load('../../../PCA/dtrajs_regspaceB_nclust-50_ticadim-8_ticalag-1_clustdim-3.npy')


# In[8]:

# convert to 
dtrajs = []
for traj in range(dtraj.shape[0]):
    dtrajs.append(dtraj[traj])


# get the mss definitions
mss_sets = []
for nstate in range(len(N_mss)):
    mss_sets.append( np.load('mss_sets_'+str(N_mss[nstate])+'states.npy') )

lag = tica_lag

# load the relevant rdf data
rdfs_nm = ['BB','AB']
cent_nm = [ '', 'B']
rdf_dir = '/data/isilon/rudzinski/cluster_tmp/LJ/Kobb_Andersen/DiffCoeff_red-dt_30ns/T'+str(Temp)+'/MSM_analysis/2016_12-Dec_13/step1-dtraj_mpi/'
data_rdf = []
sel = []
pairs_excl = []
r = []
gr = []
N_solshel = []
solshel_max = []
sig = []
rcut = []
weight = []
for rdf in range(len(rdfs_nm)):
    sys = rdfs_nm[rdf]
    data_rdf = np.load(rdf_dir+'data_'+sys+'_rdf.npz')
    sel.append(data_rdf['sel'])
    pairs_excl.append(data_rdf['pairs_excl'])
    r.append(data_rdf['r'])
    gr.append(data_rdf['gr'])
    N_solshel.append(data_rdf['N_solshel'])
    solshel_max.append(data_rdf['solshel_max'])
    sig.append(data_rdf['sig'])
    rcut.append(data_rdf['rcut'])
    weight.append(data_rdf['weight'])


# In[60]:

# sort the dtraj by mss
indir_CN = '/data/isilon/rudzinski/cluster_tmp/LJ/Kobb_Andersen/DiffCoeff_red-dt_30ns/T'+str(Temp)+'/MSM_analysis/2016_12-Dec_13/step1-dtraj_mpi/'

# first get the dtrajs
Xre = []
for rdf in range(len(rdfs_nm)):
    Xre.append([])
    pair_type = rdfs_nm[rdf]
    if (pair_type == 'AB'):
        dtraj_CNs = np.load(indir_CN+'/dtraj_CNs_'+pair_type+'_0-5ns.npz')
        Xre[rdf] = dtraj_CNs['Xre_'+cent_nm[rdf]]
        dtraj_CNs = np.load(indir_CN+'/dtraj_CNs_'+pair_type+'_5-10ns.npz')
        Xre[rdf] = np.hstack((Xre[rdf],dtraj_CNs['Xre_'+cent_nm[rdf]]))
        dtraj_CNs = np.load(indir_CN+'/dtraj_CNs_'+pair_type+'_10-15ns.npz')
        Xre[rdf] = np.hstack((Xre[rdf],dtraj_CNs['Xre_'+cent_nm[rdf]]))
        dtraj_CNs = np.load(indir_CN+'/dtraj_CNs_'+pair_type+'_15-20ns.npz')
        Xre[rdf] = np.hstack((Xre[rdf],dtraj_CNs['Xre_'+cent_nm[rdf]]))
        dtraj_CNs = np.load(indir_CN+'/dtraj_CNs_'+pair_type+'_20-25ns.npz')
        Xre[rdf] = np.hstack((Xre[rdf],dtraj_CNs['Xre_'+cent_nm[rdf]]))
        dtraj_CNs = np.load(indir_CN+'/dtraj_CNs_'+pair_type+'_25-30ns.npz')
        Xre[rdf] = np.hstack((Xre[rdf],dtraj_CNs['Xre_'+cent_nm[rdf]]))
    else:
        dtraj_CNs = np.load(indir_CN+'/dtraj_CNs_'+pair_type+'_0-5ns.npz')
        Xre[rdf] = dtraj_CNs['Xre']
        dtraj_CNs = np.load(indir_CN+'/dtraj_CNs_'+pair_type+'_5-10ns.npz')
        Xre[rdf] = np.hstack((Xre[rdf],dtraj_CNs['Xre']))
        dtraj_CNs = np.load(indir_CN+'/dtraj_CNs_'+pair_type+'_10-15ns.npz')
        Xre[rdf] = np.hstack((Xre[rdf],dtraj_CNs['Xre']))
        dtraj_CNs = np.load(indir_CN+'/dtraj_CNs_'+pair_type+'_15-20ns.npz')
        Xre[rdf] = np.hstack((Xre[rdf],dtraj_CNs['Xre']))
        dtraj_CNs = np.load(indir_CN+'/dtraj_CNs_'+pair_type+'_20-25ns.npz')
        Xre[rdf] = np.hstack((Xre[rdf],dtraj_CNs['Xre']))
        dtraj_CNs = np.load(indir_CN+'/dtraj_CNs_'+pair_type+'_25-30ns.npz')
        Xre[rdf] = np.hstack((Xre[rdf],dtraj_CNs['Xre']))
        
# now organize per mss
CN_rdf = []
avg_CN = []
for nstate in range(len(N_mss)):
    CN_rdf.append([])
    avg_CN.append([])
    for rdf in range(len(rdfs_nm)):
        CN_rdf[nstate].append([])
        avg_CN[nstate].append([])
        for shell in range(N_solshel[rdf]):
            CN_rdf[nstate][rdf].append([[] for x in range(N_mss[nstate])])
            avg_CN[nstate][rdf].append(np.zeros(N_mss[nstate]))
            for traj in range(len(Xre[rdf])):
                CN_tmp = Xre[rdf][traj][::Nprune*lag,shell]
                for state in range( N_mss[nstate] ):
                    state_frs = np.array([])
                    for micro in mss_sets[nstate][state]:
                        state_frs = np.hstack( (state_frs,np.where(dtrajs[traj] == micro)[0]) )
                    if ( len(state_frs) != 0 ):
                        CN_rdf[nstate][rdf][shell][state] = np.hstack((CN_rdf[nstate][rdf][shell][state],CN_tmp[state_frs.astype(int)]))
            for state in range( N_mss[nstate] ):
                avg_CN[nstate][rdf][shell][state] = np.mean( CN_rdf[nstate][rdf][shell][state] )    
                
np.save('CN_rdf',CN_rdf)
np.save('avg_CN', avg_CN)



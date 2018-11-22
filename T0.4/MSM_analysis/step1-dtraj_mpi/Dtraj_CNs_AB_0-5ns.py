
# coding: utf-8

# Generate CN-dtrajs of 80:20 KA sims
# ====

# In[1]:

import pyemma
pyemma.__version__


# In[2]:

import os
#get_ipython().magic(u'pylab inline')
#matplotlib.rcParams.update({'font.size': 12})


# In[3]:

import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt


# Some useful inhouse functions
# ------

# In[4]:

from CN_functions import *
from gen_dtraj_functions import *

# Add the MPI stuff
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
print 'hello from proc '+str(rank)
size = comm.Get_size()

from copy import copy, deepcopy

# Now, load filenames and topology
# ------

# In[5]:

traj_nm = '0-5ns'
indir = '/data/isilon/rudzinski/cluster_tmp/LJ/Kobb_Andersen/DiffCoeff_red-dt_30ns/T0.4'
topfile =  indir+'/LJ_1frame.pdb'
indir += '/traj'
traj_list = []
Nf = 1
t0 = 0
dt = 5
#for i in range(Nf):
#    traj_fnm = indir+'/config_series_sim_T0.4_wV_'+str(t0+i*dt)+'-'+str(t0+(i+1)*dt)+'ns.xtc' 
#    traj_list.append(traj_fnm)
traj_fnm = indir+'/config_series_sim_T0.4_'+traj_nm+'.xtc'
traj_list.append(traj_fnm)


# In[6]:

feat = coor.featurizer(topfile)


# In[7]:

# select all the indices with a particular site name, nb - site indices are 0-indexed
#print feat.topology.select("name B")
# select a particular molecule, nb - residue numbers are 1-indexed
#print feat.topology.select("residue 1")
# total number of sites
#print feat.topology.n_atoms
# total number of molecules
#print feat.topology.n_residues


# In[8]:

#traj_dt = 1 # in ps
#n_frames_p_traj = 1001
n_mol = feat.topology.n_residues
n_sites_p_mol = feat.topology.n_atoms / n_mol


# Load the rdf data
# ------

# In[9]:

# AA
sys = 'AB'
data_rdf = np.load('data_'+sys+'_rdf.npz')
sel = data_rdf['sel']
pairs_excl = data_rdf['pairs_excl']
r = data_rdf['r']
gr = data_rdf['gr']
N_solshel = data_rdf['N_solshel']
solshel_max = data_rdf['solshel_max']
sig = data_rdf['sig']
rcut = data_rdf['rcut']
weight = data_rdf['weight']
# plus one specific parameter
n_AA_mol = int(n_mol*0.8)
n_BB_mol = int(n_mol*0.2)
#n_mol_type = n_AA_mol
#n_feat_p_mol = 1
n_traj = Nf
chunk = 312

# In[10]:


# In[11]:

feat_CN, frames_p_chunk = calc_chunkwise_noavg_mpi( lambda x: calc_wghtd_CN_multss( x, (str(sel), pairs_excl, rcut, r, weight) ), traj_list, topfile, size, rank, chunk_size=chunk, dim=1, stride=1, skip=0)

# get the results from each processor
feat_CN_global = comm.gather(feat_CN, root=0)
frames_p_chunk_global = comm.gather(frames_p_chunk, root=0)
chunk_cumsum = np.hstack( (0,np.cumsum(frames_p_chunk)) )
chunk_cumsum_global = comm.gather(chunk_cumsum, root=0)
# Now back together the results in the proper order
if ( rank == 0 ):

#    feat_CN_re = []
#    # get the total number of frames
#    ctr_frames = 0
#    for proc in range(size):
#        ctr_frames += feat_CN_global[proc][0].shape[0]
#    for shell in range(N_solshel):
#        #tmp = [None]*(len(feat_CN_global[:][shell]))
#        tmp = np.zeros(shape=(ctr_frames,feat_CN_global[0][0].shape[1]))
#        for proc in range(size):
#            tmp[proc::size] = feat_CN_global[proc][shell]
#        feat_CN_re.append( tmp )

    chunk_ind = np.zeros(size,dtype='int')
    for i in range(n_traj):
        proc = 0
        if ( i == 0 ):
            feat_CN_re = feat_CN_global[proc][chunk_cumsum_global[proc][chunk_ind[proc]]:chunk_cumsum_global[proc][chunk_ind[proc]+1],:,:]
        else:
            feat_CN_re = np.vstack( (feat_CN_re,feat_CN_global[proc][chunk_cumsum_global[proc][chunk_ind[proc]]:chunk_cumsum_global[proc][chunk_ind[proc]+1],:,:]) )
        chunk_ind[proc] += 1
        proc += 1
        for chunks in range(1,len(np.concatenate(frames_p_chunk_global))/n_traj): 
            feat_CN_re = np.vstack( (feat_CN_re,feat_CN_global[proc][chunk_cumsum_global[proc][chunk_ind[proc]]:chunk_cumsum_global[proc][chunk_ind[proc]+1],:,:]) )
            chunk_ind[proc] += 1
            proc += 1
            proc = proc%size


# In[12]:

    feat_CN_re = np.array(feat_CN_re)
    n_frames_p_traj = feat_CN_re.shape[0]
    n_feat_p_mol = feat_CN_re.shape[1]
    n_mol_type = feat_CN_re.shape[2]
    Xre = feat_CN_per_mol_reonly(feat_CN_re, len(feat_CN_re), n_frames_p_traj, n_traj, n_feat_p_mol, n_mol_type)

    Xre_A = Xre[1000:]
    Xre_B = Xre[0:1000]

# In[ ]:

# save all the dtraj data
    np.savez('dtraj_CNs_'+sys+'_'+traj_nm, N_solshel=N_solshel, Xre_A=Xre_A, Xre_B=Xre_B)







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


# Now, load filenames and topology
# ------

# In[5]:

indir = '/data/isilon/rudzinski/cluster_tmp/LJ/Kobb_Andersen/DiffCoeff_red-dt_30ns/T0.4'
topfile =  indir+'/LJ_AK_5000.pdb'
indir += '/traj'
traj_list = []
Nf = 6
t0 = 0
dt = 5
for i in range(Nf):
    traj_fnm = indir+'/config_series_sim_T0.4_wV_'+str(t0+i*dt)+'-'+str(t0+(i+1)*dt)+'ns.xtc' 
    traj_list.append(traj_fnm)
#traj_fnm = indir+'/config_series_sim_T0.4_wV_test.xtc'
#traj_list.append(traj_fnm)


# In[6]:

feat = coor.featurizer(topfile)


# In[7]:

# select all the indices with a particular site name, nb - site indices are 0-indexed
#print feat.topology.select("name B")
# select a particular molecule, nb - residue numbers are 1-indexed
print feat.topology.select("residue 1")
# total number of sites
print feat.topology.n_atoms
# total number of molecules
print feat.topology.n_residues


# In[8]:

#traj_dt = 1 # in ps
#n_frames_p_traj = 1001
n_mol = feat.topology.n_residues
n_sites_p_mol = feat.topology.n_atoms / n_mol


# Load the rdf data
# ------

# In[9]:

# BB
data_rdf = np.load('data_BB_rdf.npz')
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
n_mol_type = n_BB_mol
n_feat_p_mol = 1
n_traj = Nf


# In[10]:

feat_CN = []


# In[11]:

for i in range(N_solshel):
    feat_CN_tmp, count = calc_chunkwise_noavg( lambda x: calc_wghtd_CN( x, (str(sel), pairs_excl, rcut[i], r, weight[i]) ), traj_list, topfile, chunk_size=50, dim=1)
    feat_CN.append(feat_CN_tmp)
    print 'done with solvation shell '+str(i)+' of '+str(N_solshel)


# In[12]:

n_frames_p_traj = len(feat_CN[0])
Xre = []


# In[13]:

for i in range(N_solshel):
    Xre_tmp = feat_CN_per_mol(feat_CN[i], len(feat_CN[i]), n_frames_p_traj, n_traj, n_feat_p_mol, n_mol_type)
    Xre.append(Xre_tmp)


# In[ ]:

# save all the dtraj data
np.savez('dtraj_CNs_BB', N_solshel=N_solshel, Xre=Xre)





# In[ ]:




# In[ ]:




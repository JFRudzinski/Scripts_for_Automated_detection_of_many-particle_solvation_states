
# coding: utf-8

# TICA and clustering with CN- of 80:20 KA sims
# ====

# In[1]:

import pyemma
pyemma.__version__


# In[2]:

import os
import numpy as np
#get_ipython().magic(u'pylab inline')
#matplotlib.rcParams.update({'font.size': 12})


# In[3]:

import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt
from collections import Counter
import mdtraj as md
#from CN_functions import *
#from gen_dtraj_functions import *


# Some useful inhouse functions
# ------

# In[4]:

#from plot_functions import *
import pickle


# Read in the data
# ------

# In[5]:

tica_lags = [1] 
tica_dim = 5
clust_dim = 5
nclust_max = 200
Y = []
tica_corr = []
for lag in tica_lags:
    Y.append(np.load('Y_B_ticadim-'+str(tica_dim)+'_ticalag-'+str(lag)+'.npy'))
    tica_corr.append(np.load('tica_corr_lag-'+str(lag)+'.npy'))
    #  clustering
    with open('clustering_regspaceB_ticadim-'+str(tica_dim)+'_ticalag-'+str(lag)+'_nclust-'+str(nclust_max)+'_clustdim-'+str(clust_dim)+'.pkl', 'rb') as f:
        clustering = pickle.load(f)
#Y_pca = np.load('Y_B_nclust-200_pcadim-10_clustdim-3.npy')


# In[6]:

np.save('tica_lag',tica_lags[0])
np.save('tica_dim',tica_dim)
np.save('clust_dim',clust_dim)
np.save('nclust_max',nclust_max)

nclust = clustering.n_clusters
nclust

np.save('nclust',nclust)

dtrajs = clustering.dtrajs

np.save('dtrajs_regspaceB_nclust-'+str(nclust)+'_ticadim-'+str(tica_dim)+'_ticalag-'+str(tica_lags[0])+'_clustdim-'+str(clust_dim),dtrajs)


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




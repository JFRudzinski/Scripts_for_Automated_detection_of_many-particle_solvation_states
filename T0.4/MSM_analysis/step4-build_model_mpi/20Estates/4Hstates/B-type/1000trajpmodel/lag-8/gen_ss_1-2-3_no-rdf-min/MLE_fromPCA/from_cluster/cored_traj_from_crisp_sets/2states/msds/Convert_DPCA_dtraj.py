
# coding: utf-8

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
import msmbuilder
from msmbuilder.msm.ratematrix import ContinuousTimeMSM
import scipy
from msmtools.analysis.dense.decomposition import eigenvectors, eigenvalues
import operator


# In[4]:

Nstates=2

# get the dtraj from DPCA analysis
dtraj_DPCA = np.genfromtxt('../coring/dtraj_mss_concat_'+str(Nstates)+'states_cored.dat')
dtraj_DPCA = dtraj_DPCA.astype(int)
traj_len = np.genfromtxt('../../dtraj_mss_traj_len.dat').astype(int)


# In[5]:

# just put back into traj form and save the mapping


# In[6]:

n_traj = len(dtraj_DPCA)/traj_len
dtraj_DPCA_renum = []
for traj in range(n_traj):
    dtraj_DPCA_renum.append(dtraj_DPCA[traj*traj_len:(traj+1)*traj_len])

np.save('dtraj_mss_'+str(Nstates)+'states_cored', dtraj_DPCA_renum)



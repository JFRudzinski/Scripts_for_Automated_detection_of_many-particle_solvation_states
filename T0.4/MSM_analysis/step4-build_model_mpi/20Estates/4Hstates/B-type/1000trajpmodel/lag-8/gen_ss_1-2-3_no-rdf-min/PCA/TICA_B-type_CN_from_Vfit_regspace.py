
# coding: utf-8

# TICA and clustering with CN- of 80:20 KA sims
# ====

# In[1]:

import pyemma
pyemma.__version__


# In[2]:

#import os
#get_ipython().magic(u'pylab inline')
#matplotlib.rcParams.update({'font.size': 12})


# In[3]:

import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as mplt


import numpy as np
# Some useful inhouse functions
# ------

# In[4]:

#from plot_functions import *

import pickle
def save_object(filename, obj):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

# Inputs
# Clustering
n_clusters = 50
dmin=0.55
clust_dim = 2
#
# General
n_Estates = 20
n_Hstates = 4
part_type = 'B'
Temp = 0.4
# TICA
tica_dim = 5
tica_lag = [16,20,24]
Nprune = 1
fr_max = 30006

np.save('tica_dim',tica_dim)
np.save('tica_lag',tica_lag[0])
np.save('Nprune',Nprune)
# Read in the dtrajs
# ------

# In[5]:
indir = '/data/isilon/rudzinski/cluster_tmp/LJ/Kobb_Andersen/DiffCoeff_red-dt_30ns/T'+str(Temp)+'/MSM_analysis/2016_12-Dec_13/step3-dtraj_viterbi_fit/'+str(n_Estates)+'Estates_fromMPIP/'+str(n_Hstates)+'Hstates/'+part_type+'-type/1000trajpmodel/'

N_BB_ss = [0,1,2,4,6]
lag = [ 8, 8, 8, 8, 8, 8, 8, 8, 8 ]
pair_type = 'BB'
dtraj_Vfit_ind = []
for ss in N_BB_ss:
    BB_sys_nm = pair_type+'ss'+str(ss)+'_'+str(n_Estates)+'Estates'+'_'+str(n_Hstates)+'Hstates'
    dtraj_Vfit_ind.append( np.load(indir+pair_type+'ss'+str(ss)+'/hmsm_lag-'+str(lag[ss])+'/dtraj_Vfit_opt_apruned_'+BB_sys_nm+'_lag-'+str(lag[ss])+'.npy')[:,:fr_max:Nprune] )

N_AB_ss = [0,2,3,5]
pair_type = 'AB'
cent_type = 'B'
lag = [ 8, 8, 8, 8, 8, 8, 8, 8 ]
for ss in N_AB_ss:
    AB_sys_nm = pair_type+'-'+cent_type+'ss'+str(ss)+'_'+str(n_Estates)+'Estates'+'_'+str(n_Hstates)+'Hstates'
    dtraj_Vfit_ind.append( np.load(indir+pair_type+'-'+part_type+'ss'+str(ss)+'/hmsm_lag-'+str(lag[ss])+'/dtraj_Vfit_opt_apruned_'+AB_sys_nm+'_lag-'+str(lag[ss])+'.npy')[:,:fr_max:Nprune] )

dtraj_Vfit_ind = np.array(dtraj_Vfit_ind)

dtraj_Vfit = []
for traj in range(dtraj_Vfit_ind.shape[1]):
    tmp = dtraj_Vfit_ind[0,traj]
    for ss in range(1,dtraj_Vfit_ind.shape[0]):
        tmp = np.vstack( (tmp, dtraj_Vfit_ind[ss,traj]) )
    dtraj_Vfit.append(tmp.T)

# In[6]:

np.array(dtraj_Vfit).shape

# **TICA**

pca_obj = coor.pca(dtraj_Vfit, var_cutoff=0.95)

save_object('pca_obj.pkl', pca_obj)

#plt.plot(tica_obj.eigenvalues,marker='x')
#plt.xlim([-1,20])
#plt.ylim([0.5,1])

# here we do a little trick to ensure that eigenvectors always have the same sign structure. 
# That's irrelevant to the analysis and just nicer plots - you can ignore it.
#for i in range(2):
#    if tica_obj.eigenvectors[0, i] > 0: 
#        tica_obj.eigenvectors[:, i] *= -1

Y = pca_obj.get_output() # get tica coordinates
np.save('Y.npy',Y)

# Now, do the clustering
Y_clust = []
for i in range(len(Y)):
    Y_clust.append(Y[i][:,0:clust_dim])
clustering = coor.cluster_kmeans(data=Y_clust, k=n_clusters, max_iter=50, tolerance=1e-05, stride=1)
save_object('clustering_kmeans_nclust-'+str(clustering.n_clusters)+'_clustdim-'+str(clust_dim)+'.pkl', clustering)
#clustering = coor.cluster_regspace(Y_clust,max_centers=n_clusters,dmin=dmin)
#save_object('clustering_regspace_nclust-'+str(clustering.n_clusters)+'_clustdim-'+str(clust_dim)+'.pkl', clustering)
#print 'n_clusters = '+str(clustering.n_clusters)


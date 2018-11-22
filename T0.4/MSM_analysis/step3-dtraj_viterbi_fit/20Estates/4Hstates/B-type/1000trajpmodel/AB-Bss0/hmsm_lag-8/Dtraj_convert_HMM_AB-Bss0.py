# coding: utf-8

# TICA and clustering with CN- of 80:20 KA sims
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

import numpy as np
# Some useful inhouse functions
# ------

# In[4]:

#from plot_functions import *
import viterbi_logop
import pickle
import operator
from copy import deepcopy

# Get the HMSM
# ------

# input params
pair_type = 'AB'
cent_type = 'B'
ss = 0
# HMSM params
n_Estates = 20
n_Hstates = 4
Nprune = 1
Nparam_traj = 1000
N_traj_sets = 1
Ntraj_0 = 0
Ntraj_f = 0
lag = 8
#
indir = '/data/isilon/rudzinski/cluster_tmp/LJ/Kobb_Andersen/DiffCoeff_red-dt_30ns/T0.4/MSM_analysis/2016_12-Dec_13/step2-dtraj_get_HMM/20Estates/4Hstates/B-type/1000trajpmodel/AB-Bss0/hmsm_lag-8/'

# Read in the dtrajs (i.e., the sequence data)
# ------

indir_CN = '/data/isilon/rudzinski/cluster_tmp/LJ/Kobb_Andersen/DiffCoeff_red-dt_30ns/T0.4/MSM_analysis/2016_12-Dec_13/step1-dtraj_mpi/'

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

dtraj_CN = []
for i in range(len(Xre)):
    dtraj_CN.append(Xre[i][::Nprune*lag,ss])

# find the most prob seq from each model
nstates = n_Hstates
#sys_nm = pair_type+'ss'+str(ss)+'_'+str(n_Estates)+'Estates'+'_'+str(n_Hstates)+'Hstates'
sys_nm = pair_type+'-'+cent_type+'ss'+str(ss)+'_'+str(n_Estates)+'Estates'+'_'+str(n_Hstates)+'Hstates'
fit_full = []
map_Hstate = []
for traj_frac in range(Ntraj_0,Ntraj_f+1):

    fit_full.append([])
    sys_nm_long = sys_nm+'_trajfrac-'+str(traj_frac)
    sys_nm_full = sys_nm_long+'_lag-'+str(lag)

    # get the model
    with open(indir+'HMSM_'+sys_nm_full+'.pkl', 'rb') as f:
        hmsm = pickle.load(f)

    # Get the clustering data
    with open(indir+'clustering'+sys_nm_long+'.pkl', 'rb') as f:
        clustering = pickle.load(f)

    # convert the dtrajs
    dtrajs = clustering.dtrajs
    cc = clustering.clustercenters[:,0]

    np.save('cc_'+str(traj_frac),cc)

    # determine the sorting of the hidden states
    m_cc = []
    for i in range(nstates):
        m_cc.append( np.sum( hmsm.metastable_memberships[:,i]*cc ) / np.sum(hmsm.metastable_memberships[:,i]) )
    m_cc = np.array(m_cc)

    m_cc_stack = []
    for i in range(nstates):
        m_cc_stack.append( np.hstack( (m_cc[i],i) ) )
    
    m_cc_sorted = sorted( m_cc_stack, key=operator.itemgetter(0))
    np.save('m_cc_sorted_'+str(traj_frac),m_cc_sorted)
    # get the forward mapping
    for_map = []
    lcc = np.arange(len(m_cc))
    for i in range(nstates):
        for_map.append( np.where(m_cc_sorted==lcc[i])[0][0] )
    map_Hstate.append( np.array(for_map) )

    dtraj_full = []
    for i in range(len(dtraj_CN)):
        dtraj_full.append(clustering.assign(np.reshape(dtraj_CN[i],(dtraj_CN[i].shape[0],1))))

    # Now, to the fitting
    Tmat = hmsm.transition_matrix
    Eprob = hmsm.observation_probabilities
    pi = hmsm.eigenvectors_left(k=1).T

    # initialize the viterbi solver
    Vseq = viterbi_logop.Decoder(pi, Tmat, Eprob)
    # get the soln for each traj

    for traj in range(len(dtraj_full)):

        # get the opt seq using viterbi
        fit_tmp = np.array( Vseq.Decode(dtraj_full[traj]) )
        # transform to sorted mss and save
        #fit_full[traj_frac].append( fit_tmp )
        fit_full[traj_frac].append( map_Hstate[traj_frac][fit_tmp].astype(int) )
        print 'done with traj '+str(traj)+' of '+str(len(dtraj_full))+' for traj_frac = '+str(traj_frac)

    np.save('dtraj_Vfit_'+sys_nm+'_lag-'+str(lag)+'traj_frac-'+str(traj_frac), fit_full[traj_frac])

print 'done with the fitting, trying to put it all together'
fit_full = np.array(fit_full)
fit_final = deepcopy(fit_full[0])
fit_cert = np.zeros(N_traj_sets+1)
indet_segs = []
fr_skip = 0
for traj in range(len(fit_full[0])):
    for frame in range(fit_full[0][0].shape[0]):

        # skip over frames that were already set
        if ( fr_skip > 0 ):
            fr_skip -= 1
            continue

        # for each frame, count how many predictions for each state
        state_counts = np.zeros(len(m_cc)).astype(int)
        for state in range(len(m_cc)):
            state_counts[state] = len(np.where( fit_full[:,traj,frame] == state )[0])
        
        max_ind = np.where( state_counts == np.max(state_counts) )[0]
        if ( len(max_ind) == 1 ): # there is a well-defined choice for this state
            fit_final[traj][frame] = deepcopy(max_ind[0])
            # keep track of the certainty of the predictions
            fit_cert[np.max(state_counts)] += 1
        elif ( len(max_ind) > 1 ): # there is a tie between predictions, place the transition in the middle of the undetermined region
            # find the next well determined frame
            flag_det = False
            fr2 = 1
            # before starting, check that the current frame isn't the last one
            if ( frame+1 == fit_full[0][0].shape[0] ):
                det_val = deepcopy(fit_final[traj][frame-1])
                # keep track of the certainty of the predictions
                fit_cert[np.max(state_counts)] += 1
                # make sure we are just updating the current frame
                fr2 = 0
                flag_det = True

            while ( not flag_det ):

                # for each frame, count how many predictions for each state
                state_counts_tmp = np.zeros(len(m_cc)).astype(int)
                for state_tmp in range(len(m_cc)):
                    state_counts_tmp[state_tmp] = len(np.where( fit_full[:,traj,frame+fr2] == state_tmp )[0])
                max_ind_tmp = np.where( state_counts_tmp == np.max(state_counts_tmp) )[0]
                if ( len(max_ind_tmp) == 1 ): # there is a well-defined choice for this state
                    det_val = deepcopy(max_ind_tmp[0])
                    # keep track of the certainty of the predictions
                    fit_cert[np.max(state_counts_tmp)] += 1
                    flag_det = True
                elif ( frame+fr2+1 == fit_full[0][0].shape[0] ): # check if this frame is the last one
                    # enforce a transition, this assumes a tie between two states only
                    grid = np.where( max_ind != fit_final[traj][frame-1] )[0][0]
                    det_val = deepcopy(max_ind[grid])
                    # keep track of the certainty of the predictions
                    fit_cert[np.max(state_counts_tmp)] += 1
                    flag_det = True
                else:
                    fr2 += 1

            # set the values for all the indeterminate frames
            for fr3 in range(fr2+1):
                if ( fr3 < np.ceil(fr2/2).astype(int) ):
                    fit_final[traj][frame+fr3] = deepcopy(fit_final[traj][frame-1])
                else:
                    fit_final[traj][frame+fr3] = deepcopy(det_val)

            # skip the frames just addressed
            fr_skip = deepcopy(fr2)
            indet_segs.append(fr_skip)
        else:
            raise ValueError('Problem choosing the optimal states')


np.save('dtraj_Vfit_opt'+sys_nm+'_lag-'+str(lag), fit_final)
np.save('frac_opt_pred'+sys_nm+'_lag-'+str(lag), fit_cert)
np.save('indet_segs'+sys_nm+'_lag-'+str(lag), indet_segs)
np.save('map_Hstate'+sys_nm+'_lag-'+str(lag), map_Hstate)



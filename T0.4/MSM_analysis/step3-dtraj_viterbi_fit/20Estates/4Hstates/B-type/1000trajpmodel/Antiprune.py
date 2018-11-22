# Since the various feature trajectories were filtered with different lag times, we need to recover a full trajectory for each
# we will simply fill in the gaps with the value at the previous timeframe

import numpy as np

N_Hstates = 4

N_BB_ss = 9
#lag = [ 16, 16, 16, 16, 16, 16, 16, 16, 16 ]
lag = [ 8, 8, 8, 8, 8, 8, 8, 8, 8 ]
for ss in range(0,N_BB_ss):
    path = 'BBss'+str(ss)+'/hmsm_lag-'+str(lag[ss])+'/'
    dtraj = np.load(path+'dtraj_Vfit_optBBss'+str(ss)+'_20Estates_'+str(N_Hstates)+'Hstates_lag-'+str(lag[ss])+'.npy')
    ntraj = dtraj.shape[0]
    nfr = dtraj.shape[1]
    dtraj_apruned = np.zeros(shape=(ntraj,nfr*lag[ss]))
    for traj in range(ntraj):
        for fr in range(nfr):
            dtraj_apruned[traj,fr*lag[ss]:(fr+1)*lag[ss]] = dtraj[traj,fr]*np.ones(lag[ss])
    np.save( path+'dtraj_Vfit_opt_apruned_BBss'+str(ss)+'_20Estates_'+str(N_Hstates)+'Hstates_lag-'+str(lag[ss])+'.npy', dtraj_apruned )

N_AB_ss = 8
#lag = [ 16, 16, 16, 16, 16, 16, 16, 16 ]
lag = [ 8, 8, 8, 8, 8, 8, 8, 8 ]
for ss in range(0,N_AB_ss):
    path = 'AB-Bss'+str(ss)+'/hmsm_lag-'+str(lag[ss])+'/'
    dtraj = np.load(path+'dtraj_Vfit_optAB-Bss'+str(ss)+'_20Estates_'+str(N_Hstates)+'Hstates_lag-'+str(lag[ss])+'.npy')
    ntraj = dtraj.shape[0]
    nfr = dtraj.shape[1]
    dtraj_apruned = np.zeros(shape=(ntraj,nfr*lag[ss]))
    for traj in range(ntraj):
        for fr in range(nfr):
            dtraj_apruned[traj,fr*lag[ss]:(fr+1)*lag[ss]] = dtraj[traj,fr]*np.ones(lag[ss])
    np.save( path+'dtraj_Vfit_opt_apruned_AB-Bss'+str(ss)+'_20Estates_'+str(N_Hstates)+'Hstates_lag-'+str(lag[ss])+'.npy', dtraj_apruned )


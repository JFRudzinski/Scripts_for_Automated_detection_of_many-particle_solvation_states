# we want to make sure that the various models made reasonably similar predictions

import numpy as np

N_Hstates = 2

N_BB_ss = 9
lag = [ 16, 16, 16, 16, 16, 16, 16, 16, 16 ]
for ss in range(N_BB_ss):
    print '**** type BB, ss '+str(ss)+' ****'
    print '\n'
    path = 'BBss'+str(ss)+'/hmsm_lag-'+str(lag[ss])+'/'
    dtraj = np.load(path+'dtraj_Vfit_optBBss'+str(ss)+'_25Estates_'+str(N_Hstates)+'Hstates_lag-'+str(lag[ss])+'.npy')
    Nfr_tot = dtraj.shape[0]*dtraj.shape[1]
    frac = np.load(path+'frac_opt_predBBss'+str(ss)+'_25Estates_'+str(N_Hstates)+'Hstates_lag-'+str(lag[ss])+'.npy')
    Nfr_frac = np.sum(frac)
    indet = np.load(path+'indet_segsBBss'+str(ss)+'_25Estates_'+str(N_Hstates)+'Hstates_lag-'+str(lag[ss])+'.npy')
    if ( len(indet) != 0 ):
        Nfr_indet = np.sum(indet)
    else:
        Nfr_indet = 0
    print 'Nfr_tot = '+str(Nfr_tot)
    print 'Nfr_tot - (Nfr_frac+Nfr_indet) = '+str(Nfr_tot - (Nfr_frac+Nfr_indet))
    print '\n'
    frac[len(frac)/2] += Nfr_indet # need to add the contributions from indeterminate frames
    print '\n'
    print 'frac[N_model_agree] = '
    print frac/np.sum(frac)
    if ( len(indet) != 0 ):
        hist,edges = np.histogram(indet,bins=np.arange(np.max(indet))+1)
    else:
        hist,edges = 0.,0.
    print '\n'
    print 'histogram of indeterminate segments => '
    print edges
    print hist
    print '\n'

N_AB_ss = 8
lag = [ 16, 16, 16, 16, 16, 16, 16, 16 ]
for ss in range(N_AB_ss):
    print '**** type AB-B, ss '+str(ss)+' ****'
    print '\n'
    path = 'AB-Bss'+str(ss)+'/hmsm_lag-'+str(lag[ss])+'/'
    dtraj = np.load(path+'dtraj_Vfit_optAB-Bss'+str(ss)+'_25Estates_'+str(N_Hstates)+'Hstates_lag-'+str(lag[ss])+'.npy')
    Nfr_tot = dtraj.shape[0]*dtraj.shape[1]
    frac = np.load(path+'frac_opt_predAB-Bss'+str(ss)+'_25Estates_'+str(N_Hstates)+'Hstates_lag-'+str(lag[ss])+'.npy')
    Nfr_frac = np.sum(frac)
    indet = np.load(path+'indet_segsAB-Bss'+str(ss)+'_25Estates_'+str(N_Hstates)+'Hstates_lag-'+str(lag[ss])+'.npy')
    if ( len(indet) != 0 ):
        Nfr_indet = np.sum(indet)
    else:
        Nfr_indet = 0
    print 'Nfr_tot = '+str(Nfr_tot)
    print 'Nfr_tot - (Nfr_frac+Nfr_indet) = '+str(Nfr_tot - (Nfr_frac+Nfr_indet))
    print '\n'
    frac[len(frac)/2] += Nfr_indet # need to add the contributions from indeterminate frames
    print '\n'
    print 'frac[N_model_agree] = '
    print frac/np.sum(frac)
    if ( len(indet) != 0 ):
        hist,edges = np.histogram(indet,bins=np.arange(np.max(indet))+1)
    else:
        hist,edges = 0.,0.
    print '\n'
    print 'histogram of indeterminate segments => '
    print edges
    print hist
    print '\n'

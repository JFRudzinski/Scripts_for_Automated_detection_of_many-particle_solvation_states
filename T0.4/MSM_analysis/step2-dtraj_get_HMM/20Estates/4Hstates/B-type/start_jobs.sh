Nss=7
lag=( 8 8 8 8 8 8 8 8 )
dmin=( 0.23 0.32 0.58 0.74 0.83 0.96 1.35 1.62 )
for ss in $(seq 0 ${Nss})
do
    cd AB-Bss${ss}/
    cd hmsm_lag-${lag[$ss]}/
    qsub qsub_HMM_mpi_1.sh
    cd ../../
done


Nss=8
lag=( 8 8 8 8 8 8 8 8 8 )
dmin=( 0.19 0.37 0.46 0.37 0.71 0.93 1.09 1.15 1.33 )
for ss in $(seq 0 ${Nss})
do
    cd BBss${ss}/
    cd hmsm_lag-${lag[$ss]}/
    qsub qsub_HMM_mpi_1.sh
    cd ../../
done


Temp=0.4

Nss=7
#lag=( 16 16 16 16 16 16 16 16 )
lag=( 8 8 8 8 8 8 8 8 )
for ss in $(seq 0 ${Nss})
do
    cd AB-Bss${ss}/
    cd hmsm_lag-${lag[$ss]}/
    rm 1*
    #rm tjob*
    #rm log.dat
    qsub ./qsub_HMM_mpi_1.sh
    cd ../../
done


Nss=8
#lag=( 16 16 16 16 16 16 16 16 16 )
lag=( 8 8 8 8 8 8 8 8 8 )
for ss in $(seq 0 ${Nss})
do
    cd BBss${ss}/
    cd hmsm_lag-${lag[$ss]}/
    rm 1*
    #rm tjob*
    #rm log.dat
    qsub ./qsub_HMM_mpi_1.sh
    cd ../../
done


Nss=7
lag=( 4 4 4 4 4 4 4 4 )
for ss in $(seq 0 ${Nss})
do
    mkdir AB-Bss${ss}/
    cd AB-Bss${ss}/
    mkdir hmsm_lag-${lag[$ss]}
    #cd hmsm_lag-${lag[$ss]}/
    #sed -i "/indir_CN =/c\indir_CN = '/ptmp/jrudz/Kobb_Andersen/T0.4/MSM_analysis/2016_12-Dec_13/step1-dtraj_mpi/'" ./Get_HMM_AB-Ass${ss}.py
    #sed -i "/Nign =/c\Nign = 3" ./Get_HMM_AB-Ass${ss}.py
    #sed -i "/lags =/c\lags = [${lag[$ss]}]" ./Get_HMM_AB-Ass${ss}.py
    #cp ../../../../../../../../../../T0.5/MSM_analysis/2016_12-Dec_13/step2-dtraj_get_HMM/25Estates/6Hstates/A-type/1000trajpmodel/AB-Ass${ss}/its_errors/qsub_HMM_mpi_1_DRACO.sh ./
    #rm tjob*
    #rm log.dat
    #qsub qsub_HMM_mpi_1_DRACO.sh
    cd ../
done


Nss=8
lag=( 4 4 4 4 4 4 4 4 4 )
for ss in $(seq 0 ${Nss})
do
    mkdir BBss${ss}
    cd BBss${ss}/
    mkdir hmsm_lag-${lag[$ss]}
    #cd hmsm_lag-${lag[$ss]}/
    #sed -i "/indir_CN =/c\indir_CN = '/ptmp/jrudz/Kobb_Andersen/T0.4/MSM_analysis/2016_12-Dec_13/step1-dtraj_mpi/'" ./Get_HMM_AAss${ss}.py
    #sed -i "/Nign =/c\Nign = 3" ./Get_HMM_AAss${ss}.py
    #sed -i "/lags =/c\lags = [${lag[$ss]}]" ./Get_HMM_AAss${ss}.py
    #cp ../../../../../../../../../../T0.5/MSM_analysis/2016_12-Dec_13/step2-dtraj_get_HMM/25Estates/6Hstates/A-type/1000trajpmodel/AAss${ss}/its_errors/qsub_HMM_mpi_1_DRACO.sh ./
    #rm tjob*
    #rm log.dat
    #qsub qsub_HMM_mpi_1_DRACO.sh
    cd ../
done


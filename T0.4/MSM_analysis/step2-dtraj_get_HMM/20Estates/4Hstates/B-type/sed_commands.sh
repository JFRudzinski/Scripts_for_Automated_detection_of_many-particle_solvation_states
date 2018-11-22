Nss=7
NHstate=5
lag=( 8 8 8 8 8 8 8 8 )
dmin=( 0.23 0.32 0.58 0.74 0.83 0.96 1.35 1.62 )
for ss in $(seq 0 ${Nss})
do
    cd AB-Bss${ss}/
    mv hmsm_lag-XX hmsm_lag-${lag[$ss]}
    cd hmsm_lag-${lag[$ss]}/
    #cp ../../../../../../25Estates/${NHstate}Hstates/B-type/1000trajpmodel/AB-Bss${ss}/hmsm_lag-*/Get_HMM_AB-Bss${ss}.py ./
    #sed -i "/indir_CN =/c\indir_CN = '/ptmp/jrudz/Kobb_Andersen/T0.4/MSM_analysis/2016_12-Dec_13/step1-dtraj_mpi/'" ./Get_HMM_AB-Bss${ss}.py
    #sed -i "/Nign =/c\Nign = 3" ./Get_HMM_AB-Bss${ss}.py
    sed -i "/n_Estates =/c\n_Estates = 20" ./Get_HMM_AB-Bss${ss}.py
    sed -i "/lags =/c\lags = [${lag[$ss]}]" ./Get_HMM_AB-Bss${ss}.py
    sed -i "/dmin =/c\dmin = ${dmin[$ss]}" ./Get_HMM_AB-Bss${ss}.py
    #cp ../../../../../../../../../../T0.5/MSM_analysis/2016_12-Dec_13/step2-dtraj_get_HMM/25Estates/6Hstates/B-type/1000trajpmodel/AB-Bss${ss}/its_errors/qsub_HMM_mpi_1_DRACO.sh ./
    #rm tjob*
    #rm log.dat
    #qsub qsub_HMM_mpi_1_DRACO.sh
    cd ../../
done


Nss=8
lag=( 8 8 8 8 8 8 8 8 8 )
dmin=( 0.19 0.37 0.46 0.37 0.71 0.93 1.09 1.15 1.33 )
for ss in $(seq 0 ${Nss})
do
    cd BBss${ss}/
    mv hmsm_lag-XX hmsm_lag-${lag[$ss]}
    cd hmsm_lag-${lag[$ss]}/
    #cp ../../../../../../25Estates/${NHstate}Hstates/B-type/1000trajpmodel/BBss${ss}/hmsm_lag-*/Get_HMM_BBss${ss}.py ./
    #sed -i "/indir_CN =/c\indir_CN = '/ptmp/jrudz/Kobb_Andersen/T0.4/MSM_analysis/2016_12-Dec_13/step1-dtraj_mpi/'" ./Get_HMM_BBss${ss}.py
    #sed -i "/Nign =/c\Nign = 3" ./Get_HMM_BBss${ss}.py
    sed -i "/n_Estates =/c\n_Estates = 20" ./Get_HMM_BBss${ss}.py
    sed -i "/lags =/c\lags = [${lag[$ss]}]" ./Get_HMM_BBss${ss}.py
    sed -i "/dmin =/c\dmin = ${dmin[$ss]}" ./Get_HMM_BBss${ss}.py
    #cp ../../../../../../../../../../T0.5/MSM_analysis/2016_12-Dec_13/step2-dtraj_get_HMM/25Estates/6Hstates/B-type/1000trajpmodel/BBss${ss}/its_errors/qsub_HMM_mpi_1_DRACO.sh ./
    #rm tjob*
    #rm log.dat
    #qsub qsub_HMM_mpi_1_DRACO.sh
    cd ../../
done


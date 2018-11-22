Nss=7
N_Hstates=4
lag=( 8 8 8 8 8 8 8 8 )
lagnew=( 8 8 8 8 8 8 8 8 )
#lagnew=( 4 4 4 4 4 4 4 4 )
dmin=( 0.23 0.32 0.58 0.74 0.83 0.96 1.35 1.62 )
for ss in $(seq 0 ${Nss})
do
    cd AB-Bss${ss}/
    #mkdir hmsm_lag-${lagnew[$ss]}
    cd hmsm_lag-${lagnew[$ss]}/
    cp ../../../../../5Hstates/B-type/1000trajpmodel/AB-Bss${ss}/hmsm_lag-${lag[$ss]}/Get_HMM_* ./
    cp ../../../../../5Hstates/B-type/1000trajpmodel/AB-Bss${ss}/hmsm_lag-${lag[$ss]}/qsub* ./
    #sed -i "/indir_CN =/c\indir_CN = '/ptmp/jrudz/Kobb_Andersen/T0.4/MSM_analysis/2016_12-Dec_13/step1-dtraj_mpi/'" ./Get_HMM_AB-Bss${ss}.py
    #sed -i "/Nign =/c\Nign = 3" ./Get_HMM_AB-Bss${ss}.py
    #sed -i "/n_Estates =/c\n_Estates = 20" ./Get_HMM_AB-Bss${ss}.py
    sed -i "/n_Hstates =/c\n_Hstates = ${N_Hstates}" ./Get_HMM_AB-Bss${ss}.py
    sed -i "/lags =/c\lags = [${lagnew[$ss]}]" ./Get_HMM_AB-Bss${ss}.py
    #sed -i "/Ntraj_0 =/c\Ntraj_0 = 0" ./Get_HMM_AB-Bss${ss}.py
    #sed -i "/Ntraj_f =/c\Ntraj_f = 3" ./Get_HMM_AB-Bss${ss}.py
    #qsub qsub_HMM_mpi_1_DRACO.sh
    cd ../../
done


Nss=8
lag=( 8 8 8 8 8 8 8 8 8 )
lagnew=( 8 8 8 8 8 8 8 8 8 )
#lagnew=( 4 4 4 4 4 4 4 4 4 )
dmin=( 0.19 0.37 0.46 0.37 0.71 0.93 1.09 1.15 1.33 )
for ss in $(seq 0 ${Nss})
do
    cd BBss${ss}/
    #mkdir hmsm_lag-${lagnew[$ss]}
    cd hmsm_lag-${lagnew[$ss]}/
    cp ../../../../../5Hstates/B-type/1000trajpmodel/BBss${ss}/hmsm_lag-${lag[$ss]}/Get_HMM_* ./
    cp ../../../../../5Hstates/B-type/1000trajpmodel/BBss${ss}/hmsm_lag-${lag[$ss]}/qsub* ./
    #sed -i "/indir_CN =/c\indir_CN = '/ptmp/jrudz/Kobb_Andersen/T0.4/MSM_analysis/2016_12-Dec_13/step1-dtraj_mpi/'" ./Get_HMM_BBss${ss}.py
    #sed -i "/Nign =/c\Nign = 3" ./Get_HMM_BBss${ss}.py
    #sed -i "/n_Estates =/c\n_Estates = 20" ./Get_HMM_BBss${ss}.py
    sed -i "/n_Hstates =/c\n_Hstates = ${N_Hstates}" ./Get_HMM_BBss${ss}.py
    sed -i "/lags =/c\lags = [${lagnew[$ss]}]" ./Get_HMM_BBss${ss}.py
    #sed -i "/Ntraj_0 =/c\Ntraj_0 = 0" ./Get_HMM_BBss${ss}.py
    #sed -i "/Ntraj_f =/c\Ntraj_f = 3" ./Get_HMM_BBss${ss}.py
    #qsub qsub_HMM_mpi_1_DRACO.sh
    cd ../../
done


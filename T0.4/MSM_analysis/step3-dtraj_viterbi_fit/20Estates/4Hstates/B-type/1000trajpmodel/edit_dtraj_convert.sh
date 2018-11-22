Temp=0.4
nHstates=4
part_type='B'

Nss=7
#lag=( 16 16 16 16 16 16 16 16 )
lag=( 8 8 8 8 8 8 8 8 )
for ss in $(seq 0 ${Nss})
do
    cd AB-Bss${ss}/
    cd hmsm_lag-${lag[$ss]}/
    cp ../../../../../5Hstates/B-type/1000trajpmodel/AB-Bss${ss}/hmsm_lag-8/Dtraj_convert_HMM_AB-Bss${ss}.py ./
    cp ../../../../../5Hstates/B-type/1000trajpmodel/AB-Bss${ss}/hmsm_lag-8/qsub* ./
    sed -i "/n_Estates =/c\n_Estates = 20" ./Dtraj_convert_HMM_AB-Bss${ss}.py
    sed -i "/n_Hstates =/c\n_Hstates = ${nHstates}" ./Dtraj_convert_HMM_AB-Bss${ss}.py
    sed -i "/lag =/c\lag = ${lag[$ss]}" ./Dtraj_convert_HMM_AB-Bss${ss}.py
    sed -i "/indir =/c\indir = '/data/isilon/rudzinski/cluster_tmp/LJ/Kobb_Andersen/DiffCoeff_red-dt_30ns/T${Temp}/MSM_analysis/2016_12-Dec_13/step2-dtraj_get_HMM/20Estates/${nHstates}Hstates/${part_type}-type/1000trajpmodel/AB-Bss${ss}/hmsm_lag-${lag[$ss]}/'" ./Dtraj_convert_HMM_AB-Bss${ss}.py
    sed -i "/indir_CN =/c\indir_CN = '/data/isilon/rudzinski/cluster_tmp/LJ/Kobb_Andersen/DiffCoeff_red-dt_30ns/T${Temp}/MSM_analysis/2016_12-Dec_13/step1-dtraj_mpi/'" ./Dtraj_convert_HMM_AB-Bss${ss}.py
    sed -i "/Nign =/c\Nign = 7" ./Dtraj_convert_HMM_AB-Bss${ss}.py
    sed -i "/Ntraj_0 =/c\Ntraj_0 = 0" ./Dtraj_convert_HMM_AB-Bss${ss}.py
    cp ../../viterbi_logop.py ./
    cd ../../
done


Nss=8
#lag=( 16 16 16 16 16 16 16 16 16 )
lag=( 8 8 8 8 8 8 8 8 8 )
for ss in $(seq 0 ${Nss})
do
    cd BBss${ss}/
    cd hmsm_lag-${lag[$ss]}/
    cp ../../../../../5Hstates/B-type/1000trajpmodel/BBss${ss}/hmsm_lag-8/Dtraj_convert_HMM_BBss${ss}.py ./
    cp ../../../../../5Hstates/B-type/1000trajpmodel/BBss${ss}/hmsm_lag-8/qsub* ./
    sed -i "/n_Estates =/c\n_Estates = 20" ./Dtraj_convert_HMM_BBss${ss}.py
    sed -i "/n_Hstates =/c\n_Hstates = ${nHstates}" ./Dtraj_convert_HMM_BBss${ss}.py
    sed -i "/lag =/c\lag = ${lag[$ss]}" ./Dtraj_convert_HMM_BBss${ss}.py
    sed -i "/indir =/c\indir = '/data/isilon/rudzinski/cluster_tmp/LJ/Kobb_Andersen/DiffCoeff_red-dt_30ns/T${Temp}/MSM_analysis/2016_12-Dec_13/step2-dtraj_get_HMM/20Estates/${nHstates}Hstates/${part_type}-type/1000trajpmodel/BBss${ss}/hmsm_lag-${lag[$ss]}/'" ./Dtraj_convert_HMM_BBss${ss}.py
    sed -i "/indir_CN =/c\indir_CN = '/data/isilon/rudzinski/cluster_tmp/LJ/Kobb_Andersen/DiffCoeff_red-dt_30ns/T${Temp}/MSM_analysis/2016_12-Dec_13/step1-dtraj_mpi/'" ./Dtraj_convert_HMM_BBss${ss}.py
    sed -i "/Nign =/c\Nign = 7" ./Dtraj_convert_HMM_BBss${ss}.py
    sed -i "/Ntraj_0 =/c\Ntraj_0 = 0" ./Dtraj_convert_HMM_BBss${ss}.py
    cp ../../viterbi_logop.py ./
    cd ../../
done


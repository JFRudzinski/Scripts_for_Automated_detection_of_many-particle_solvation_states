Temp=0.4

Nss=7
#lag=( 16 16 16 16 16 16 16 16 )
lag=( 8 8 8 8 8 8 8 8 8 )
for ss in $(seq 0 ${Nss})
do
    cp driver_mpi.py AB-Bss${ss}/hmsm_lag-${lag[$ss]}/
    cd AB-Bss${ss}/
    cd hmsm_lag-${lag[$ss]}/
    sed -i "/lag =/c\lag = ${lag[$ss]}" ./driver_mpi.py
    sed -i "/sys_nm =/c\sys_nm = 'AB-Bss${ss}'" ./driver_mpi.py
    sed -i "/script_path =/c\    script_path = './Dtraj_convert_HMM_'+sys_nm+'.py'" ./driver_mpi.py
    cp ../../qsub_HMM_mpi_1.sh ./
    sed -i "/mpiexec -n 8/c\/home/theorie/rudzinski/soft/anaconda/envs/PyEmma-new/bin/mpiexec -n 8 /home/theorie/rudzinski/soft/anaconda/envs/PyEmma-new/bin/python Dtraj_convert_HMM_AB-Bss${ss}.py" ./qsub_HMM_mpi_1.sh
    cd ../../
done


Nss=8
#lag=( 16 16 16 16 16 16 16 16 16 )
lag=( 8 8 8 8 8 8 8 8 8 )
for ss in $(seq 0 ${Nss})
do
    cp driver_mpi.py BBss${ss}/hmsm_lag-${lag[$ss]}/
    cd BBss${ss}/
    cd hmsm_lag-${lag[$ss]}/
    echo $ss
    sed -i "/lag =/c\lag = ${lag[$ss]}" ./driver_mpi.py
    sed -i "/sys_nm =/c\sys_nm = 'BBss${ss}'" ./driver_mpi.py
    sed -i "/script_path =/c\    script_path = './Dtraj_convert_HMM_'+sys_nm+'.py'" ./driver_mpi.py
    cp ../../qsub_HMM_mpi_1.sh ./
    sed -i "/mpiexec -n 8/c\/home/theorie/rudzinski/soft/anaconda/envs/PyEmma-new/bin/mpiexec -n 8 /home/theorie/rudzinski/soft/anaconda/envs/PyEmma-new/bin/python Dtraj_convert_HMM_BBss${ss}.py" ./qsub_HMM_mpi_1.sh
    cd ../../
done


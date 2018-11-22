Nproc=4
clustering="/data/isilon/rudzinski/soft_backup/Clustering-0.12_static/clustering"

dtraj='../../dtraj_mss_concat_6states.dat'
traj_len=30006

lag=16
mkdir data_wtd_lag-${lag}
echo "* ${lag}" &> win
${clustering} coring -s $dtraj -w win -d wtd -v --concat-nframes ${traj_len}
mv wtd_* data_wtd_lag-${lag}

lag=32
mkdir data_wtd_lag-${lag}
echo "* ${lag}" &> win
${clustering} coring -s $dtraj -w win -d wtd -v --concat-nframes ${traj_len}
mv wtd_* data_wtd_lag-${lag}

lag=48
mkdir data_wtd_lag-${lag}
echo "* ${lag}" &> win
${clustering} coring -s $dtraj -w win -d wtd -v --concat-nframes ${traj_len}
mv wtd_* data_wtd_lag-${lag}

lag=64
mkdir data_wtd_lag-${lag}
echo "* ${lag}" &> win
${clustering} coring -s $dtraj -w win -d wtd -v --concat-nframes ${traj_len}
mv wtd_* data_wtd_lag-${lag}

lag=80
mkdir data_wtd_lag-${lag}
echo "* ${lag}" &> win
${clustering} coring -s $dtraj -w win -d wtd -v --concat-nframes ${traj_len}
mv wtd_* data_wtd_lag-${lag}

lag=96
mkdir data_wtd_lag-${lag}
echo "* ${lag}" &> win
${clustering} coring -s $dtraj -w win -d wtd -v --concat-nframes ${traj_len}
mv wtd_* data_wtd_lag-${lag}

lag=14
mkdir data_wtd_lag-${lag}
echo "* ${lag}" &> win
${clustering} coring -s $dtraj -w win -d wtd -v --concat-nframes ${traj_len}
mv wtd_* data_wtd_lag-${lag}

lag=32
mkdir data_wtd_lag-${lag}
echo "* ${lag}" &> win
${clustering} coring -s $dtraj -w win -d wtd -v --concat-nframes ${traj_len}
mv wtd_* data_wtd_lag-${lag}


#${clustering} coring -s $dtraj -w win -o dtraj_mss_concat_2states_cored.dat -v --concat-nframes 30006

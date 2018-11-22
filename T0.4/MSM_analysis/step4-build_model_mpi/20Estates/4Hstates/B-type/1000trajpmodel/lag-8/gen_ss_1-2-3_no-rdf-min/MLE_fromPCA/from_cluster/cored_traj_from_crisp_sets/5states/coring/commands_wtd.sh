Nproc=4
clustering="/data/isilon/rudzinski/soft_backup/Clustering-0.12_static/clustering"

dtraj='../../dtraj_mss_concat_5states.dat'
traj_len=30006

lag=2
mkdir data_wtd_lag-${lag}
echo "* ${lag}" &> win
${clustering} coring -s $dtraj -w win -d wtd -v --concat-nframes ${traj_len}
mv wtd_* data_wtd_lag-${lag}

lag=4
mkdir data_wtd_lag-${lag}
echo "* ${lag}" &> win
${clustering} coring -s $dtraj -w win -d wtd -v --concat-nframes ${traj_len}
mv wtd_* data_wtd_lag-${lag}

lag=6
mkdir data_wtd_lag-${lag}
echo "* ${lag}" &> win
${clustering} coring -s $dtraj -w win -d wtd -v --concat-nframes ${traj_len}
mv wtd_* data_wtd_lag-${lag}

lag=8
mkdir data_wtd_lag-${lag}
echo "* ${lag}" &> win
${clustering} coring -s $dtraj -w win -d wtd -v --concat-nframes ${traj_len}
mv wtd_* data_wtd_lag-${lag}

lag=10
mkdir data_wtd_lag-${lag}
echo "* ${lag}" &> win
${clustering} coring -s $dtraj -w win -d wtd -v --concat-nframes ${traj_len}
mv wtd_* data_wtd_lag-${lag}

lag=12
mkdir data_wtd_lag-${lag}
echo "* ${lag}" &> win
${clustering} coring -s $dtraj -w win -d wtd -v --concat-nframes ${traj_len}
mv wtd_* data_wtd_lag-${lag}

lag=14
mkdir data_wtd_lag-${lag}
echo "* ${lag}" &> win
${clustering} coring -s $dtraj -w win -d wtd -v --concat-nframes ${traj_len}
mv wtd_* data_wtd_lag-${lag}

lag=16
mkdir data_wtd_lag-${lag}
echo "* ${lag}" &> win
${clustering} coring -s $dtraj -w win -d wtd -v --concat-nframes ${traj_len}
mv wtd_* data_wtd_lag-${lag}


#${clustering} coring -s $dtraj -w win -o dtraj_mss_concat_2states_cored.dat -v --concat-nframes 30006

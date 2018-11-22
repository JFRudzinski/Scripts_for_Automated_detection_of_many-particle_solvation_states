Nproc=4
clustering="/data/isilon/rudzinski/soft_backup/Clustering-0.12_static/clustering"

dtraj='../../dtraj_mss_concat_5states.dat'

#${clustering} coring -s $dtraj -w win -d wtd -v --concat-nframes 30006

${clustering} coring -s $dtraj -w win -o dtraj_mss_concat_5states_cored.dat -v --concat-nframes 30006

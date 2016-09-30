#/bin/bash

export PRECISION=1
#make header
for i in `seq 2 2 16` `seq 20 4 128` #`seq 136 8 512`
do
    for j in `seq 2 2 16` `seq 20 4 128` #`seq 136 8 512`
    do
        for k in `seq 2 2 16` `seq 20 4 128` #`seq 136 8 512`
        do
            export BLOCKDEFS="-DS_BLOCKSIZE_MR=4 -DS_BLOCKSIZE_NR=4 -DS_BLOCKSIZE_MC=$((4*$i)) -DS_BLOCKSIZE_KC=$((4*$j)) -DS_BLOCKSIZE_NC=$((32*$k)) -DPRECISION=1"
#            make runperf
            export BLOCKDEFS="-DS_BLOCKSIZE_MR=6 -DS_BLOCKSIZE_NR=6 -DS_BLOCKSIZE_MC=$((6*$i)) -DS_BLOCKSIZE_KC=$((6*$j)) -DS_BLOCKSIZE_NC=$((36*$k)) -DPRECISION=1"
#            make runperf
        done
    done
done

export PRECISION=2
make header
for i in `seq 2 2 16` `seq 20 4 128` #`seq 136 8 512`
do
    for j in `seq 2 2 16` `seq 20 4 128` #`seq 136 8 512`
    do
        for k in `seq 2 2 16` `seq 20 4 128` #`seq 136 8 512`
        do
            export BLOCKDEFS="-DD_BLOCKSIZE_MR=4 -DD_BLOCKSIZE_NR=4 -DD_BLOCKSIZE_MC=$((4*$i)) -DD_BLOCKSIZE_KC=$((4*$j)) -DD_BLOCKSIZE_NC=$((32*$k)) -DPRECISION=2"
            make runperf
        done
    done
done


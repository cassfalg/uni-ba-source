#/bin/bash

make header
#for i in 1 2 3 4 5 6 7 8
for i in 4 5 6 7 8 9
do
    export S_BLOCKSIZE_MR=$i
#    for j in 1 2 3 4 5 6 7 8
    for j in 4 5 6 7 8 9
    do
        export S_BLOCKSIZE_NR=$j
#        for k in `seq 1 32` 40 48 56 64 72 80 88 96
        for k in `seq 6 32` 40 48 56 64 72 80 88 96
        do
            export S_BLOCKSIZE_MC=$(($i*$k))
            export S_BLOCKSIZE_KC=$(($j*$k))
            export S_BLOCKSIZE_NC=$((64*$k))
#            echo $S_BLOCKSIZE_MR $S_BLOCKSIZE_NR $S_BLOCKSIZE_MC $S_BLOCKSIZE_KC $S_BLOCKSIZE_NC
            export BLOCKDEFS="-DS_BLOCKSIZE_MR=$i -DS_BLOCKSIZE_NR=$j -DS_BLOCKSIZE_MC=$(($i*$k)) -DS_BLOCKSIZE_KC=$(($j*$k)) -DS_BLOCKSIZE_NC=$((4*$i*$j*$k))"
            make runperf
        done
    done
done


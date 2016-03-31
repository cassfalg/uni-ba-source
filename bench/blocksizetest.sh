#/bin/bash

make header
#for i in 1 2 3 4 5 6 7 8
for i in 2 3 4 5 6 7 8
do
    export D_BLOCKSIZE_MR=$i
#    for j in 1 2 3 4 5 6 7 8
    for j in 2 3 4 5 6 7 8
    do
        export D_BLOCKSIZE_NR=$j
#        for k in `seq 1 32` 40 48 56 64 72 80 88 96
        for k in `seq 4 32` 40 48 56 64 72 80 88 96
        do
            export D_BLOCKSIZE_MC=$(($i*$k))
            export D_BLOCKSIZE_KC=$(($j*$k))
            export D_BLOCKSIZE_NC=$((64*$k))
            echo $D_BLOCKSIZE_MR $D_BLOCKSIZE_NR $D_BLOCKSIZE_MC $D_BLOCKSIZE_KC $D_BLOCKSIZE_NC
            export BLOCKDEFS="-DD_BLOCKSIZE_MR=$i -DD_BLOCKSIZE_NR=$j -DD_BLOCKSIZE_MC=$(($i*$k)) -DD_BLOCKSIZE_KC=$(($j*$k)) -DD_BLOCKSIZE_NC=$((64*$k))"
            make runperf
        done
    done
done


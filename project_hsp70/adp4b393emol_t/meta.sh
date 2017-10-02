#!/bin/bash


run_sopnucleo1 () {

    RNUM1=$((1 + RANDOM % 10))
    RNUM2=$((RANDOM % 10))
    RNUM3=$(($RNUM2*10 + $RNUM1))
    echo $RNUM3

    mkdir -p Struct_data Coord
    run_sopnucleo_emol ../adp4b393.pdb ../top4b3934e.top $RNUM3
}

# for i in `seq 1 10`;
# do
#     run_sopnucleo1 $i
# done

run_sopnucleo1

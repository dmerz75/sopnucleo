#!/bin/bash


run_sopnucleo1 () {

    RNUM1=$((1 + RANDOM % 10))
    RNUM2=$((RANDOM % 10))
    RNUM3=$(($RNUM2*10 + $RNUM1))
    echo $RNUM3

    mkdir -p Struct_data Coord
    run_sopnucleo ../hsp70atpmod.pdb ../Contact_map_hsp70total $RNUM3
}

# for i in `seq 1 10`;
# do
#     run_sopnucleo1 $i
# done

run_sopnucleo1

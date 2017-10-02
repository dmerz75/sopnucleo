#!/bin/bash


run_sopnucleo () {

    mkdir -p Struct_data Coord
    run_sopnucleo hsp70atp.pdb Contact_map_intra_b_70atp 1

}


run_sopnucleo

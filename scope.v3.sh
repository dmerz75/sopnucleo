#!/bin/bash

# find . -name *.cpp
# find . -name *.h
# comment src/deprecated/*

# cscope -b -q -R
# cscope -b -q -f cscope.files
# cscope -b -q -k -R -i cscope.files


preparation () {
    # echo '' > cscope.files
    find . -name "*.c*" ! -path "./exclude*" > cscope.files
    find . -name "*.h*" ! -path "./exclude*" >> cscope.files
    # find ./src -name "*.c*" >> cscope.files
    # find ./include -name "*.h*" >> cscope.files
    echo "Makefile" >> cscope.files
}

scope () {
    cscope -b -i cscope.files
}


preparation
scope

#!/bin/sh

# this is an ugly, ugly hack for compiling the material models, but I am lazy
# should put this into the makefile at some point...

if [ $# -ne 1 ]; then
    echo "Needs 1 argument."
    echo "ARG1: C source file with material function 'calculate_stress'."
else
    CFILE=$1
    SOFILE=${CFILE%.c}.so

    gcc -g -O3 -Wall -march=native -fPIC -shared -nostartfiles -o $SOFILE $CFILE
fi


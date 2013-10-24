#!/bin/sh

if [ $# -ne 2 ] && [ $# -ne 1 ]; then
    echo "Needs 1 or 2 arguments.";
    echo "ARG1: Time for sim.";
    echo "ARG2: (Optional) Job name. If left unused, name is created from current time.";
else
    if [ $# -ne 2 ]; then
        JOB_NAME=job_`date +%Y%m%d%H%M%S`
    else
        JOB_NAME=$2
    fi;
    echo $JOB_NAME
    make clean
    if make -j4; then
        mkdir -p jobs/$JOB_NAME
        octave gen_particles.m && time --verbose ./mpm_2d -o jobs/$JOB_NAME -g generated_grid.txt -p generated_particles.txt $1
        tar --exclude=jobs --exclude-vcs -cvzf jobs/$JOB_NAME.tar.gz ../$(basename `pwd`)
    fi
fi


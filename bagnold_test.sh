#!/bin/sh

MIN_SIZE=19
STEP=1
MAX_SIZE=41

FLOWSTART=0.1

MODEL_SO=g_local_mu2.so
THREADS=1
SIMTIME=500.0


EMOD=1e9
NUMOD=0.3

# global preprocessing

for opening_size in `seq $MIN_SIZE $STEP $MAX_SIZE`;
do {
    # preprocessing
    OUTDIR=/home/$USER/backup/forward/bagnold-sc-$opening_size
    mkdir -p $OUTDIR

    cat << EOF > $OUTDIR/runcfg
timestep
{
    dt-max = 1e-2
    dt-min = 1e-10
    dt = 3e-6
    automatic-dt = 0 
    allow-dt-increase = 0
    stable-dt-threshold = 4
}

solver
{
    solver-type = explicit-usl
}

material
{
    material-file = "./$MODEL_SO"
    use-builtin = 0
    properties = {$EMOD, $NUMOD}
        # properties are Young's modulus and Poisson's ratio
    integer-properties = { }
        # no integer properties by default
}

boundary-conditions
{
    boundary-conditions-file = "builtin.so"
    use-builtin = 1
    properties = {-1, $FLOWSTART}
    integer-properties = {$opening_size }
}

implicit
{
    displacement-norm-ratio = 1e-2
    residual-norm-ratio = 1e-2

    converged-displacement-norm = 1e-8

    unstable-iteration-count = 10
}

input
{
    initial-particle-file = "column_40d_horizontal.txt"
    grid-file = "grid_200d.txt"
}

output
{
    directory = "$OUTDIR/"
    user = ${USER:-unknown}
    sample-rate = 60.0
}
EOF

    # run simulations
    # echo $opening_size;
    # ./mpm_2d -c $OUTDIR/runcfg -t $THREADS $SIMTIME

    # postprocessing
}; done;
    OUTDIR=

parallel --progress ./mpm_2d -c /home/$USER/backup/forward/bagnold-sc-{}/runcfg -t$THREADS $SIMTIME  ::: `seq $MIN_SIZE $STEP $MAX_SIZE`;

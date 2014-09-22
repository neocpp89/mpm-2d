#!/bin/sh

MIN_SIZE=0.01
STEP=0.01
MAX_SIZE=0.10

FLOWSTART=0.1

MODEL_SO=g_local_mu2.so
THREADS=4
SIMTIME=10.0

EMOD=1e9
NUMOD=0.3

# global preprocessing
octave beverloo_genp.m

for opening_size in `seq $MIN_SIZE $STEP $MAX_SIZE`;
do {
    # preprocessing
    OUTDIR=/home/mit/backup/forward-param/bev-$opening_size
    #mkdir -p $OUTDIR

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
    properties = {$opening_size, $FLOWSTART}
    integer-properties = { }
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
    initial-particle-file = "generated_particles.txt"
    grid-file = "generated_grid.txt"
}

output
{
    directory = "$OUTDIR/"
    user = ${USER:-unknown}
    sample-rate = 60.0
}
EOF

    # run simulations
    echo $opening_size;
    ./mpm_2d -c $OUTDIR/runcfg -t $THREADS $SIMTIME

    # postprocessing
}; done;



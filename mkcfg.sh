#!/bin/sh

for S in `seq 0 0.5 3`
do
for V in `seq 2.0 0.25 5.25`
do
mkdir "$S"_"$V".output
cat > "$S"_"$V".cfg  << __EOF__
timestep
{
    dt-max = 1e-2
    dt-min = 1e-10
    dt = 3.0e-6

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
    material-file = "g_local_mu2_plane_strain.so"
    use-builtin = 0
    properties = {1e8, 0.45, 0.005, -$V}
        # properties are Young's modulus and Poisson's ratio
    integer-properties = { }
        # no integer properties by default
}

boundary-conditions
{
    boundary-conditions-file = "box.so"
    use-builtin = 0
    properties = {0.0, 0.0, 1.0, 1.0}
        # for the box, the properties are:
        #   (bottomleftx, bottomlefty, width, height)
    integer-properties = { 1, 3, 1, 0 }
        # also for the box, the integer properties are the
        # fixed displacment DOF:
        #   0 - no friction
        #   1 - x only
        #   2 - y only
        #   3 - x and y both
        # the properties are listed in order for the
        #   (left, bottom, right, top)
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
    initial-particle-file = "s$S.particles"
    grid-file = "s$S.grid"
}

output
{
    directory = "${S}_${V}.output"
    user = ${USER:-unknown}
    sample-rate = 60.0
}
__EOF__
done;
done;

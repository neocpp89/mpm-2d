#!/bin/sh

if [ $# -le 4 ]; then echo "Need at least 4 arguments. "; exit; fi;

f=$1
s=$2
m=$3

echo "Using file $f."
echo "Plotting dataset x,y,$s"
echo "Estimate for max value of $s is $m."

shift 3

for arg in $@; do ./streaming_frame_plotter_scatter.py $f "x,y,$s" $arg "0,$m"; done;


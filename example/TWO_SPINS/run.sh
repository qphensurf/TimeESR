#!/bin/bash

rundir=$(pwd)
progdir="../../src"
./clean.sh

cd $progdir
make clean
make

cd $rundir
$progdir/TimeESR.x

gnuplot POPULATIONS.plot
gnuplot SPIN.plot
gnuplot CURRENT.plot

xdg-open Populations_1-4.png &
xdg-open SPINx.png &
xdg-open SPINy.png &
xdg-open SPINz.png &
xdg-open Current.png &

python bloch_sphere.py
xdg-open bloch_sphere.png &
xdg-open bloch_movie.mp4 &

#!/bin/bash

rundir=$(pwd)
progdir="/home/qns/Research/TimESR/Release/src"
./clean.sh

cd $progdir
make clean
make

cd $rundir
$progdir/TimeESR.x

gnuplot POPULATIONS.plot
gnuplot SPIN.plot
gnuplot CURRENT.plot

xdg-open Populations_1-3.png &
xdg-open SPINx.png &
xdg-open SPINy.png &
xdg-open SPINz.png &
xdg-open Current.png &

python bloch_sphere.py

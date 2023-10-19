#!/bin/bash

rundir=$(pwd)
progdir="../../src"
./clean.sh

cd $progdir
make clean
make

cd $rundir
$progdir/TimeESR.x

names='$\downarrow$ $\uparrow$ $\emptyset$'
python plot_spinpops.py 3 POPULATIONS.dat spinpops.png $names

xdg-open spinpops.png &

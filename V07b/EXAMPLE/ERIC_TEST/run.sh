#!/bin/bash

rundir=$(pwd)
progdir="./TimESR/V07b"
./clean.sh

cp $progdir/arch/make.sys.gnu $progdir/make.sys
cd $progdir
make clean
make

cd $rundir
$progdir/TimeESR.x

gnuplot POPULATIONS.plot
gnuplot SPIN.plot

xdg-open Populations_1-3.png &
xdg-open SPINx.png &
xdg-open SPINy.png &
xdg-open SPINz.png &

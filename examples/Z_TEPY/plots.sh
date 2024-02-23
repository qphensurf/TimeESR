#!/bin/bash

names='$\downarrow\downarrow$ $\uparrow\downarrow$ $\downarrow\uparrow$ $\uparrow\uparrow$ $\emptyset\downarrow$ $\emptyset\uparrow$'
numstates=6
totaldim=8
numlines=3
lines='223.904447808894 279.8805597611175 289.204578409156'
timediff=200

echo "Plotting Current"
python plot_current.py C.dat current.png $timediff $numlines $lines

echo "Plotting Populations"
python plot_pops.py $numstates POPULATIONS.dat populations.png $names $timediff $numlines $lines

echo "Plotting Denisty Matrix Elements"
python plot_dmelements.py $numstates RHO.dat $totaldim $timediff $numlines $lines

echo "Plotting Fidelity"
python plot_fidelity.py RHO.dat $totaldim fidelity.png $timediff $numlines $lines

echo "Plotting Concurrence"
python plot_concurrence.py RHO.dat $totaldim concurrence.png $timediff $numlines $lines

echo "Plotting Negativity"
python plot_negativity.py RHO.dat $totaldim negativity.png $timediff $numlines $lines


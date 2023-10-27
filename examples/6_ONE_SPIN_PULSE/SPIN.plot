set terminal pngcairo enhanced font 'Verdana,18' #set terminal svg size 350,262 fname 'Verdana' fsize 14
set grid
set xlabel "Time (ns)"
#set xrange [0:500]

set output 'SPINx.png'
set ylabel "Spin Projection x"
#set yrange [-0.6:0.6]
#set ytics -0.6,0.2,0.6
unset key
p "SpinDynamics.dat" u 1:2 w l title "Spin 1"

set output 'SPINy.png'
set ylabel "Spin Projection y"
#set yrange [-0.6:0.6]
#set ytics -0.6,0.2,0.6
unset key
p "SpinDynamics.dat" u 1:3 w l title "Spin 1"

set output 'SPINz.png'
set ylabel "Spin Projection z"
#set yrange [-0.6:0.6]
#set ytics -0.6,0.2,0.6
unset key
p "SpinDynamics.dat" u 1:4 w l title "Spin 1"


set terminal pngcairo enhanced font 'Verdana,18' #set terminal svg size 350,262 fname 'Verdana' fsize 14
set grid
set xlabel "Time (ns)"
#set xrange [0:500]

set output 'SPINx.png'
set ylabel "Sx"
#set yrange [-0.6:0.6]
#set ytics -0.6,0.2,0.6
#unset key
p "SpinDynamics.dat" u 1:2 w l title "Spin 1", "SpinDynamics.dat" u 1:5 w l title "Spin 2"

set output 'SPINy.png'
set ylabel "Sy"
#set yrange [-0.6:0.6]
#set ytics -0.6,0.2,0.6
#unset key
p "SpinDynamics.dat" u 1:3 w l title "Spin 1", "SpinDynamics.dat" u 1:6 w l title "Spin 2"

set output 'SPINz.png'
set ylabel "Sz"
#set yrange [-0.6:0.6]
#set ytics -0.6,0.2,0.6
#unset key
p "SpinDynamics.dat" u 1:4 w l title "Spin 1", "SpinDynamics.dat" u 1:7 w l title "Spin 2"

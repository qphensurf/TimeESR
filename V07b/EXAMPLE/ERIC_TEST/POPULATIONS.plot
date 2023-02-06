set terminal pngcairo enhanced font 'Verdana,18' #set terminal svg size 350,262 fname 'Verdana' fsize 14
set grid
set xlabel "Time (ns)"

set output 'Populations_1-3.png'
set ylabel "Population"
#set yrange [-0.1:1.0]
#set ytics 0.0,0.1,1.0
#set xrange [0:500]
#unset key
p "POPULATIONS.dat" u 1:2  w l title "State 1", "POPULATIONS.dat" u 1:3  w l title "State 2", "POPULATIONS.dat" u 1:4  w l title "State 3"

set output 'Populations_4-6.png'
set ylabel "Population"
#set yrange [-0.1:1.0]
#set ytics 0.0,0.1,1.0
#set xrange [0:500]
#unset key
p "POPULATIONS.dat" u 1:5  w l title "State 4", "POPULATIONS.dat" u 1:6  w l title "State 5", "POPULATIONS.dat" u 1:7  w l title "State 6"




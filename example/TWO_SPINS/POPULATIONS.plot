set terminal pngcairo enhanced font 'Verdana,18' #set terminal svg size 350,262 fname 'Verdana' fsize 14
set grid
set xlabel "Time (ns)"

set output 'Populations_1-4.png'
set ylabel "Population"
#set yrange [-0.1:1.0]
#set ytics 0.0,0.1,1.0
#set xrange [0:500]
#unset key
p "POPULATIONS.dat" u 1:2  w l title "State 1", "POPULATIONS.dat" u 1:3  w l title "State 2", "POPULATIONS.dat" u 1:4  w l title "State 3", "POPULATIONS.dat" u 1:5  w l title "State 4"

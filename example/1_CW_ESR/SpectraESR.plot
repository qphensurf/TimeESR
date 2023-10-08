set terminal pngcairo enhanced font 'Verdana,18' #set terminal svg size 350,262 fname 'Verdana' fsize 14
set grid
set xlabel "Frequency (GHz)"

set output 'SpectraESR.png'
set ylabel "Current (pA)"
#set yrange [-0.1:1.0]
#set ytics 0.0,0.1,1.0
#set xrange [0:500]
unset key
p "SpectraESR.dat" u 1:2  w l

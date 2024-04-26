set terminal pngcairo size 1200,1000 font ",15"
set output "magn_time.png"

set multiplot layout 4,1

set xlabel "Time (ps)"
set ylabel "M/M_s (adim.)"

array labels[3]=["Ox", "Oy", "Oz"]

set xrange [0:3700]

timescale=1e12
lss=3
set key t o c horizontal 
do for [i=2:4]{
 p "output.txt" every 1 u ($1*timescale):i w lp ls i lw lss ps 0.5 title sprintf("%s",labels[i-1])}

set ylabel "B_L (T)"
p "output.txt" every 1 u ($1*timescale):11 w lp ls 5 lw lss ps 0.5 title "Laser Field || O_z"

reset

set terminal gif animate delay 20

set output 'shallow_water_serial.gif'

set xrange [1:100]
set yrange [-0.2:1.2]

set xlabel 'x [m]' font ',16'
set ylabel 'height [m]' font ',16'

set object circle at 90,0.9 size scr 0.1 fc 'yellow' fillstyle solid

do for [i=0:5000:20] {
  set title sprintf('Shallow water simulation, time = %4.1f s',0.02*i) font ",18"
  p 'shallow_serial.dat' index i with filledcurve y1=-0.2 notitle lw 2 lc 'cyan'
}

unset output

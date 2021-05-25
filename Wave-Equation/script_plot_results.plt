f(x) = x*(1-x**2)

plot [0:1] f(x) w l lw 2 lc "red" title "u(x,0) = x(1-x^2)"
replot "data_to_plot.dat" using 1:2 title "Série de Fourier" w l lw 2 lc "blue"

set terminal pngcairo
set output "Resultado.png"
set title "Solução em t = 0 s"
set grid
set xlabel "x (m)"
set ylabel "u (m)" 
replot
set terminal wxt
set output
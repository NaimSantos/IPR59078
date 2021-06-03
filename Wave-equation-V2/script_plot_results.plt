f(x) = x*(1-x**2)

set key left top # colocar a legenda na esquerda

plot [0:1] f(x) w l lw 2 lc "red" title "u(x,0) = x(1-x^2)"
replot "dados_N5.dat" using 1:2 title "Série (N= 5)" w l lw 2 lc "blue"

set terminal pngcairo
set output "Resultado.png"
set title "Solução em t = 0 s"
set grid
set xlabel "x (m)"
set ylabel "u (m)" 
replot
set terminal wxt
set output
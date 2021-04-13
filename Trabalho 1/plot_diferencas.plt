plot "dados.dat" using 1:2 title "Implícito" w l lw 3 lc "red" dt '_', "dados.dat" using 1:3 title "Crank-Nicolson" w l lw 3 lc "blue"
set terminal pngcairo
set output "Graf1.png"
set grid
#set title "Perfil de Temperatura em t = 500 s"
set xlabel "Comprimento (m)"
set ylabel "Temperatura ({\260}C)" 
replot
set terminal wxt
set output


# with lines = w l
# linewidth = lw
# linecolor = lc
# dashtype = dt
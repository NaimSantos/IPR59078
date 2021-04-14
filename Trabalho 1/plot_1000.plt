plot "dados.dat" using 1:2 title "Imp-Dif" w l lw 2 lc "red" dt '_', "dados.dat" using 1:4 title "Imp-Fic" w l lw 2 lc "red", "dados.dat" using 1:3 title "CN-Imp" w l lw 2 lc "blue", "dados.dat" using 1:5 title "CN-Fic" w l lw 2 lc "blue" dt '_'
set terminal pngcairo
set output "Graf3.png"
set grid
#set title "Perfil de Temperatura em t = 1000 s"
set xlabel "Comprimento (m)"
set ylabel "Temperatura ({\260}C)" 
replot
set terminal wxt
set output


# with lines = w l
# linewidth = lw
# linecolor = lc
# dashtype = dt
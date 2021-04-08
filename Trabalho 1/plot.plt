plot "dados.dat" using 1:2 title "Implícito (Dif. Finit.)" with lines linecolor "red", "dados.dat" using 1:3 title "C-N (Dif. Finit.)" with lines linecolor "blue", "dados.dat" using 1:4 title "Implícito (Nós Ficticios)" with lines linecolor "dark-green", "dados.dat" using 1:5 title "C-N (Nós Ficticios)" with lines linecolor "orange"
set terminal pngcairo
set output "Temp500s.png"
set grid
set title "Perfil de Temperatura em t = 500 s"
set xlabel "Comprimento (m)"
set ylabel "Temperatura (C)" 
replot
set terminal wxt
set output

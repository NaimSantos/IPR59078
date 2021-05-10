plot "Plot_Analitico.dat" using 1:2 title "Solucao Analitica" w l lw 2 lc "red",
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
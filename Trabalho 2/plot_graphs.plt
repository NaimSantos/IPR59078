##
plot "dados_parciais.dat" using 1:6 title "Analitico (500 s)" w l lw 2 lc "red"
replot	"dados_diferencas_finitas_500.dat" using 1:2 title "Dif-Fin (500 s)" w l lw 2 lc "blue" dt '_'
replot "dados_parciais.dat" using 1:7 title "Analitico (750 s)" w l lw 2 lc "dark-green"
replot	"dados_diferencas_finitas_750.dat" using 1:2 title "Dif-Fin (750 s)" w l lw 2 lc "purple" dt '_'
replot	"dados_parciais.dat" using 1:8 title "Analitico (1000 s)" w l lw 2 lc "black"
replot	"dados_diferencas_finitas_1000.dat" using 1:2 title "Dif-Fin (1000 s)" w l lw 2 lc "yellow" dt '_'
set terminal pngcairo
set output "Grafico_Total.png"
set grid
set xlabel "Comprimento (m)"
set ylabel "Temperatura ({\260}C)" 
replot
set terminal wxt
set output

##
plot "dados_L2.dat" using 1:2 title "x = L/2 (m)" w l lw 2 lc "red"
set terminal pngcairo
set output "Grafico_L2.png"
set grid
set xlabel "Tempo (s)"
set ylabel "Temperatura ({\260}C)" 
replot
set terminal wxt
set output

##
plot "dados_parciais.dat" using 1:6 title "Analitico (500 s)" w l lw 2 lc "red", "dados_diferencas_finitas_500.dat" using 1:2 title "Dif-Fin (500 s)" w l lw 2 lc "blue" dt '_'
set terminal pngcairo
set output "Grafico_500.png"
set grid
#set title "Perfis de Temperatura Obtidos"
set xlabel "Comprimento (m)"
set ylabel "Temperatura ({\260}C)" 
replot
set terminal wxt
set output

##
plot "dados_parciais.dat" using 1:7 title "Analitico (750 s)" w l lw 2 lc "red", "dados_diferencas_finitas_750.dat" using 1:2 title "Dif-Fin (750 s)" w l lw 2 lc "blue" dt '_'
set terminal pngcairo
set output "Grafico_750.png"
set grid
#set title "Perfis de Temperatura Obtidos"
set xlabel "Comprimento (m)"
set ylabel "Temperatura ({\260}C)" 
replot
set terminal wxt
set output

##
plot "dados_parciais.dat" using 1:8 title "Analitico (1000 s)" w l lw 2 lc "red", "dados_diferencas_finitas_1000.dat" using 1:2 title "Dif-Fin (1000 s)" w l lw 2 lc "blue" dt '_'
set terminal pngcairo
set output "Grafico_1000.png"
set grid
#set title "Perfis de Temperatura Obtidos"
set xlabel "Comprimento (m)"
set ylabel "Temperatura ({\260}C)" 
replot
set terminal wxt
set output



# with lines = w l
# linewidth = lw
# linecolor = lc
# dashtype = dt
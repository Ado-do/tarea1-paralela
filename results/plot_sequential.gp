# CSV input
set terminal pdfcairo
set datafile separator ","
set key autotitle columnhead

# Etiquetas y grid
#set title "Multiplicación de Matrices Secuenciales"
set xlabel "Tamaño de Matrices (n)"
set ylabel "Tiempo de Ejecución (segundos)"
set key top left
set grid

# Escala log
set logscale xy
set format x "%g"

# Graficar tiempos
#set terminal qt 0
#plot for [i=2:4] "experimentos_sequenciales.csv" using 1:i:xtic(1) with linespoints linewidth 2

#set terminal qt 0
set output "sequential_cache.pdf"
plot for [i=3:6] "experimentos_sequenciales.csv" using 1:i:xtic(1) with linespoints linewidth 2

#set terminal qt 1
set output "sequential_strassen.pdf"
plot for [i=7:10] "experimentos_sequenciales.csv" using 1:i:xtic(1) with linespoints linewidth 2

#set terminal qt 2
set output "sequential.pdf"
plot "experimentos_sequenciales.csv" using 1:2:xtic(1) with linespoints linewidth 2, \
    "" using 1:5:xtic(1) with linespoints linewidth 2, \
    "" using 1:7:xtic(1) with linespoints linewidth 2

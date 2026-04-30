# Configuración general
set terminal pdfcairo
set datafile separator ','

TARGET_N = 4096

set grid

# Tiempo de Ejecución
set output 'parallel_time.pdf'
# set title sprintf("Tiempo de Ejecución vs Procesadores (N = %d)", TARGET_N)
set xlabel "Número de Hilos / Procesadores (p)"
set ylabel "Tiempo (segundos)"
set key top right

# Columnas: 3=p, 4=tiled, 5=strassen, 6=hybrid
plot 'parallel_metrics.csv' using 3:($2==TARGET_N ? $4 : 1/0) with linespoints lw 2 pt 7 title 'Tiled', \
     ''                     using 3:($2==TARGET_N ? $5 : 1/0) with linespoints lw 2 pt 5 title 'Strassen', \
     ''                     using 3:($2==TARGET_N ? $6 : 1/0) with linespoints lw 2 pt 9 title 'Hybrid'

# Aceleración (Speedup)
set output 'parallel_speedup.pdf'
# set title sprintf("Aceleración (Speedup) vs Procesadores (N = %d)", TARGET_N)
set ylabel "Speedup"
set key top left

# Columnas: 3=p, 7=tiled_S, 9=strassen_S, 11=hybrid_S
plot 'parallel_metrics.csv' using 3:($2==TARGET_N ? $7 : 1/0) with linespoints lw 2 pt 7 title 'Tiled', \
     ''                     using 3:($2==TARGET_N ? $9 : 1/0) with linespoints lw 2 pt 5 title 'Strassen', \
     ''                     using 3:($2==TARGET_N ? $11: 1/0) with linespoints lw 2 pt 9 title 'Hybrid', \
     x with lines dashtype 2 lc rgb 'black' title 'Speedup Ideal'

# Eficiencia
set output 'parallel_eficciency.pdf'
# set title sprintf("Eficiencia vs Procesadores (N = %d)", TARGET_N)
set ylabel "Eficiencia"
set yrange [0:1.2]
set key top right

# Columnas: 3=p, 8=tiled_E, 10=strassen_E, 12=hybrid_E
plot 'parallel_metrics.csv' using 3:($2==TARGET_N ? $8 : 1/0) with linespoints lw 2 pt 7 title 'Tiled', \
     ''                     using 3:($2==TARGET_N ? $10: 1/0) with linespoints lw 2 pt 5 title 'Strassen', \
     ''                     using 3:($2==TARGET_N ? $12: 1/0) with linespoints lw 2 pt 9 title 'Hybrid', \
     1.0 with lines dashtype 2 lc rgb 'black' title 'Eficiencia Ideal'

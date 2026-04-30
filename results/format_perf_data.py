import pandas as pd

df_seq = pd.read_csv('perf_sequential.csv')
df_par = pd.read_csv('perf_parallel.csv')

df_seq['miss_rate'] = (df_seq['cache_misses'] / df_seq['cache_references']) * 100
df_par['miss_rate'] = (df_par['cache_misses'] / df_par['cache_references']) * 100

# PREPARAR DATOS PARA GNUPLOT
# Tasa de Fallas (Clásica vs Bloques)
g1 = "perf_plot1.csv"
df_g1 = df_seq[df_seq['algo'].isin(['classic', 'cache'])].pivot(index='n', columns='algo', values='miss_rate')
# df_g1.to_csv('plot1_missrate.csv', sep=',')
df_g1.to_csv(g1)

# Ciclos y Fallas (Clásica vs Bloques)
g2 = "perf_plot2.csv"
df_g2 = df_seq[df_seq['algo'].isin(['classic', 'cache'])].pivot(index='n', columns='algo', values=['cycles', 'cache_misses'])
df_g2.columns = [f"{col[0]}_{col[1]}" for col in df_g2.columns] # Aplana los nombres de columnas
# df_g2.to_csv('plot2_cycles_misses.csv', sep='\t')
df_g2.to_csv(g2)

# Ciclos e Instrucciones
N = 2048
g3 = "perf_plot3.csv"
df_g3 = df_seq[df_seq['n'] == N][['algo', 'cycles', 'instructions']]
# df_g3.to_csv('plot3_bars.dat', sep='\t', index=False)
df_g3.to_csv(g3)

# Gráfico 4: Secuencial vs Paralelo (Bloques/Cache)
g4 = "perf_plot4.csv"
df_seq_cache = df_seq[df_seq['algo'] == 'cache'][['n', 'cycles', 'cache_references']].set_index('n')
df_par_tiled = df_par[df_par['algo'] == 'tiled'][['n', 'cycles', 'cache_references']].set_index('n')
df_g4 = df_seq_cache.join(df_par_tiled, lsuffix='_seq', rsuffix='_par')
# df_g4.to_csv('plot4_seq_vs_par.csv', sep='\t')
df_g4.to_csv(g4)

print(f"Archivos csv generados para Gnuplot: {g1} {g2} {g3} {g4}")

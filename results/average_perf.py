import pandas as pd
import sys

if len(sys.argv) < 2:
    print(f"Uso: {sys.argv[0]} <perf_type.csv>")
    sys.exit(1)

perf = sys.argv[1]
df = pd.read_csv(perf)
new_df = df.groupby(["n", "algo"]).mean().reset_index()
new_df.to_csv(perf, index=False)

print(f"Archivo {perf} promediado.")

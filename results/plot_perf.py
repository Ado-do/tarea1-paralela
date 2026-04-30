import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def plot_csv():
    #sns.set_theme(style="whitegrid")

    seq = "perf_sequential.csv"
    par = "perf_parallel.csv"

    if not os.path.exists(seq) or not os.path.exists(par):
        print(f"Error: {seq} y {par} no generados.")
        return

    df_seq = pd.read_csv("perf_sequential.csv")
    df_par = pd.read_csv("perf_parallel.csv")

    # 1. Tasa de Fallos de Caché (Clásica vs Bloques)
    df_seq["miss_rate"] = (df_seq["cache_misses"] / df_seq["cache_references"]) * 100
    df_cache_comp = df_seq[df_seq["algo"].isin(["classic", "cache"])] # Filtrar algos

    plt.figure(figsize=(8, 5))
    sns.barplot(data=df_cache_comp, x="n", y="miss_rate", hue="algo")

    # plt.title("Comparación de Tasa de Fallos de Caché: Clásica vs Bloques")
    plt.ylabel("Tasa de Fallos (%)")
    plt.xlabel("Tamaño de Matriz (n)")
    plt.legend(title="Algoritmo")
    plt.tight_layout()
    plt.savefig("perf1_missrate.pdf", format="pdf", dpi=300)
    plt.close()

    # 2. Rendimiento por Algoritmo (n = 2048) (log)
    N = 2048
    df_2048 = df_seq[df_seq["n"] == N].copy()

    # Transformar datos usando melt para agrupar ciclos e instrucciones
    df_melted = df_2048.melt(id_vars="algo", value_vars=["cycles", "instructions"], var_name="métrica", value_name="cantidad")

    # Renombrar para mejor visualización en la leyenda
    # df_melted["métrica"] = df_melted["métrica"].replace({"cycles": "Ciclos de CPU", "instructions": "Instrucciones"})

    plt.figure(figsize=(9, 6))
    sns.barplot(data=df_melted, x="algo", y="cantidad", hue="métrica")

    n_val = df_2048["n"].iloc[0]
    # plt.title(f"Comparación de Ciclos e Instrucciones por Algoritmo (n = {n_val})")
    plt.ylabel("Cantidad (log)")
    plt.xlabel("Variante del Algoritmo")
    plt.yscale("log")

    plt.tight_layout()
    plt.yscale("log")
    plt.savefig("perf2_scale.pdf", format="pdf", dpi=300)
    plt.close()


    # 3. Transición Secuencial vs Paralelo
    df_seq_comp = df_seq.copy()
    df_seq_comp["modo"] = "secuencial"
    df_seq_comp["algo"] = df_seq_comp["algo"].replace({"cache": "bloques"})

    df_par_comp = df_par.copy()
    df_par_comp["modo"] = "paralelo"
    df_par_comp["algo"] = df_par_comp["algo"].replace({"tiled": "bloques"})

    # Unir ambos datasets
    df_combined = pd.concat([df_seq_comp, df_par_comp])
    df_combined = df_combined[df_combined["n"] >= 512]

    # Filtrar Strassen y Bloques, n >= 512
    df_focus = df_combined[df_combined["algo"].isin(["strassen", "bloques"])]

    # 2 subgraficos
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # 1) Strassen
    sns.barplot(data=df_focus[df_focus["algo"] == "strassen"], x="n", y="cache_misses", hue="modo", ax=axes[0])
    axes[0].set_title("Strassen")
    axes[0].set_ylabel("Fallos de Caché Absolutos")
    axes[0].set_xlabel("Tamaño de Matriz (n)")
    #axes[0].set_yscale("log")

    # 2) Bloques/Tiled
    sns.barplot(data=df_focus[df_focus["algo"] == "bloques"], x="n", y="cache_misses", hue="modo", ax=axes[1])
    axes[1].set_title("Cache-friendly")
    axes[1].set_ylabel("Fallos de Caché Absolutos")
    axes[1].set_xlabel("Tamaño de Matriz (n)")
    #axes[1].set_yscale("log")

    plt.tight_layout()
    plt.savefig("perf3_seqvspar.pdf", format="pdf", dpi=300)
    plt.close()

    print("Listo, archivos perf pdf generados.")


if __name__ == "__main__":
    plot_csv()

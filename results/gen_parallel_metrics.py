import pandas as pd
import sys


def main():
    if len(sys.argv) != 2:
        print("Uso: python generate_metrics.py <paralell.csv>")
        sys.exit(1)

    input_csv = sys.argv[1]
    output_csv = input_csv.replace(".csv", "_metrics.csv")

    df = pd.read_csv(input_csv)
    algos = df.columns[2:]

    for algo in algos:
        base_cases = df[df["p"] == 1]

        t1_mapping = dict(zip(base_cases["n"], base_cases[algo]))

        df[f"{algo}_speedup"] = df["n"].map(t1_mapping) / df[algo]

        df[f"{algo}_efficiency"] = df[f"{algo}_speedup"] / df["p"]

    df.to_csv(output_csv)
    print(f"Metricas guardadas en: {output_csv}")


if __name__ == "__main__":
    main()

#!/bin/bash

# Config
EXE="./build/sequential_main"
OUT_CSV="perf_sequential.csv"
SIZES=(128 256 512 1024 2048)
ALGOS=("classic" "cache" "strassen")
EVENTS="cycles,instructions,cache-references,cache-misses"

# Header CSV
echo "n,algo,cycles,instructions,cache_references,cache_misses" > $OUT_CSV

echo "** Comenzando profiling..."

for n in "${SIZES[@]}"; do
    for algo in "${ALGOS[@]}"; do
        echo "Profiling: n=$n, algo=$algo"
        
        RAW_DATA=$(perf stat -x, -e $EVENTS $EXE $n $algo 2>&1)

        # Extraemos usando grep y cut (formato perf -x: value,,event,...)
        CYCLES=$(echo "$RAW_DATA" | grep "cycles" | cut -d, -f1)
        INSTR=$(echo "$RAW_DATA" | grep "instructions" | cut -d, -f1)
        REFS=$(echo "$RAW_DATA" | grep "cache-references" | cut -d, -f1)
        MISSES=$(echo "$RAW_DATA" | grep "cache-misses" | cut -d, -f1)

        # Pal CSV
        echo "$n,$algo,$CYCLES,$INSTR,$REFS,$MISSES" >> $OUT_CSV
    done
done

echo "** Experimantos finalizados. Resultados guardados a $OUT_CSV"

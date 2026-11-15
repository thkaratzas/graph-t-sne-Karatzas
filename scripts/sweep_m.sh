#!/bin/bash
# Sweep the leaf capacity m for the sequential Barnes-Hut.
# Run from the repository root:  bash scripts/sweep_m.sh
# Output: results/sweep_m.csv

GRAPHS="k49_norm_10NN com-Amazon coPapersDBLP"
MS="1 2 4 8 16 32"
REPS=3
OUT=results/sweep_m.csv

echo "tag,graph,n,m,ncells,depth,build,com,walk,zerr,ferr,rep" > $OUT

for g in $GRAPHS; do
    for m in $MS; do
        for r in $(seq 1 $REPS); do
            echo "$g  m=$m  rep=$r"
            ./bin/seq data/$g $m | grep "^CSV_BH" | sed "s/$/,$r/" >> $OUT
        done
    done
done

echo "done: $OUT"
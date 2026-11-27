#!/bin/bash
# Correctness numbers of the report: errors vs reference, tree diagnostics,
# agreement of the four versions. Run from the repository root.
# Output: results/check-<machine>.txt
OUT=results/check-$(hostname -s).txt
P=4

{
for g in k49_norm_10NN com-Amazon coPapersDBLP; do
    echo "===== $g"
    ./bin/seq data/$g | grep -E "beta|err|tree|leaves|Z bh|Z direct|BH vs"
    OMP_NUM_THREADS=$P ./bin/omp data/$g > /dev/null
    CILK_NWORKERS=$P ./bin/cilk data/$g > /dev/null
    PT_NUM_THREADS=$P ./bin/pthreads data/$g > /dev/null
    for v in omp cilk pthreads; do
        ./bin/compare results/$g-$v.bin results/$g-seq.bin
    done
done
} | tee $OUT

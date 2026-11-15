#!/bin/bash
# Timing of every phase, for every version, graph and thread count.
# Run from the repository root:
#   bash scripts/measure.sh                         all graphs
#   bash scripts/measure.sh "k49_norm_10NN com-Amazon"   only these
# Output: results/times-<machine>.csv and results/machine-<machine>.txt

GRAPHS=${1:-"k49_norm_10NN com-Amazon coPapersDBLP"}
REPS=3
CHUNK=256
MACHINE=$(hostname)
OUT=results/times-$MACHINE.csv

# thread counts: 1 2 4 8 ... up to the number of cores, and the number of cores itself
NCORES=$(nproc)
THREADS=""
p=1
while [ $p -lt $NCORES ]; do
    THREADS="$THREADS $p"
    p=$((p * 2))
done
THREADS="$THREADS $NCORES"

# the machine, for the report
{
    echo "== date";  date
    echo "== lscpu"; lscpu
    echo "== memory"; free -h
    echo "== gcc";   gcc --version | head -1
    echo "== clang"; clang --version | head -1
    echo "== flags"; grep -E "^CFLAGS" Makefile
} > results/machine-$MACHINE.txt

# header only if the file is new
if [ ! -f $OUT ]; then
    echo "version,graph,threads,schedule,rep,P,attr,build,com,walk,grad,direct" > $OUT
fi

# one run: prints one CSV line from the program's output
run() {   # run <version> <graph> <threads> <schedule> <rep>
    ./bin/$1 data/$2 | awk -v pre="$1,$2,$3,$4,$5" '
        /^P:/          { P = $2 }
        /^attr:/       { A = $2 }
        /^build:/      { B = $2; C = $5; W = $8 }
        /^grad:/       { G = $2 }
        /^direct rep:/ { D = $3 }
        END { if (D == "") D = "NA";
              print pre "," P "," A "," B "," C "," W "," G "," D }' >> $OUT
}

for g in $GRAPHS; do
    for r in $(seq 1 $REPS); do
        echo "seq $g rep=$r"
        run seq $g 1 none $r
    done

    for p in $THREADS; do
        for r in $(seq 1 $REPS); do
            for s in dynamic static; do
                if [ $s = dynamic ]; then os="dynamic,$CHUNK"; else os="static"; fi

                echo "omp $g p=$p $s rep=$r"
                OMP_NUM_THREADS=$p OMP_SCHEDULE=$os run omp $g $p $s $r

                echo "pthreads $g p=$p $s rep=$r"
                PT_NUM_THREADS=$p PT_SCHEDULE=$s PT_CHUNK=$CHUNK run pthreads $g $p $s $r
            done

            echo "cilk $g p=$p rep=$r"
            CILK_NWORKERS=$p run cilk $g $p ws $r
        done
    done
done

echo "done: $OUT"
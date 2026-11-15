#!/bin/bash
# Timing of every phase, for every version, graph and thread count.
# Run from the repository root:
#   bash scripts/measure.sh                              all graphs
#   bash scripts/measure.sh "k49_norm_10NN com-Amazon"   only these
# Output: results/times-<machine>.csv and results/machine-<machine>.txt

GRAPHS=${1:-"k49_norm_10NN com-Amazon coPapersDBLP"}
REPS=3
CHUNK=256
MACHINE=$(hostname -s)
OUT=results/times-$MACHINE.csv

# thread counts: every count up to 16 cores, powers of two above that
NCORES=$(getconf _NPROCESSORS_ONLN)
if [ $NCORES -le 16 ]; then
    THREADS=$(seq 1 $NCORES)
else
    THREADS=""
    p=1
    while [ $p -lt $NCORES ]; do
        THREADS="$THREADS $p"
        p=$((p * 2))
    done
    THREADS="$THREADS $NCORES"
fi

# the machine, for the report
{
    echo "== date"
    date
    if [ "$(uname)" = Darwin ]; then
        echo "== cpu"
        sysctl -n machdep.cpu.brand_string
        sysctl hw.perflevel0.physicalcpu hw.perflevel1.physicalcpu
        sysctl hw.l1dcachesize hw.perflevel0.l2cachesize hw.perflevel1.l2cachesize
        echo "== memory (bytes)"
        sysctl -n hw.memsize
        echo "== os"
        sw_vers
    else
        echo "== lscpu"
        lscpu
        echo "== memory"
        free -h
    fi
    echo "== compilers and flags"
    make -s info
} > results/machine-$MACHINE.txt 2>&1

# header only if the file is new
if [ ! -f $OUT ]; then
    echo "version,graph,threads,schedule,rep,P,attr,build,com,walk,grad,direct" > $OUT
fi

# one run: appends one CSV line from the program's output
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
    echo "warm-up $g (not recorded)"
    ./bin/seq data/$g > /dev/null

    # sequential at the start ...
    for r in $(seq 1 $REPS); do
        echo "seq $g start rep=$r"
        run seq $g 1 start $r
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

    # ... and at the end: if these are slower, the machine got hot (throttling)
    for r in $(seq 1 $REPS); do
        echo "seq $g end rep=$r"
        run seq $g 1 end $r
    done
done

echo "done: $OUT"

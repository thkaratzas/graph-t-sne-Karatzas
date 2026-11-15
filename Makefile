CC     = gcc
CILKCC = clang
CFLAGS = -O2 -Wall -Istarter -Isrc
LIBS   = -lmatio -lm

COMMON = src/p_common.c src/io.c src/bh_common.c src/quadtree.c starter/data.c
SEQ    = src/main_seq.c src/p_seq.c src/attr_seq.c src/rep_direct.c src/bh_seq.c src/grad_seq.c
OMP    = src/main_omp.c src/p_omp.c src/attr_omp.c src/rep_direct_omp.c src/quadtree_omp.c src/bh_omp.c src/grad_omp.c
CILK   = src/main_cilk.c src/p_cilk.c src/attr_cilk.c src/rep_direct_cilk.c src/quadtree_cilk.c src/bh_cilk.c src/grad_cilk.c

all: bin/seq bin/omp bin/cilk bin/compare


bin/cilk: $(CILK) $(COMMON) src/tsne.h
	mkdir -p bin
	$(CILKCC) $(CFLAGS) -fopencilk -o $@ $(CILK) $(COMMON) $(LIBS)

bin/seq: $(SEQ) $(COMMON) src/tsne.h
	mkdir -p bin
	$(CC) $(CFLAGS) -o $@ $(SEQ) $(COMMON) $(LIBS)

bin/omp: $(OMP) $(COMMON) src/tsne.h
	mkdir -p bin
	$(CC) $(CFLAGS) -fopenmp -o $@ $(OMP) $(COMMON) $(LIBS)


bin/compare: src/compare.c
	mkdir -p bin
	$(CC) $(CFLAGS) -o $@ src/compare.c -lm

clean:
	rm -rf bin
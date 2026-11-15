CC     = gcc
CFLAGS = -O2 -Wall -Istarter -Isrc
LIBS   = -lmatio -lm

COMMON = src/p_common.c src/io.c starter/data.c
SEQ    = src/main_seq.c src/p_seq.c src/attr_seq.c src/rep_direct.c src/quadtree.c src/bh_seq.c src/grad_seq.c

bin/seq: $(SEQ) $(COMMON) src/tsne.h
	mkdir -p bin
	$(CC) $(CFLAGS) -o $@ $(SEQ) $(COMMON) $(LIBS)

clean:
	rm -rf bin
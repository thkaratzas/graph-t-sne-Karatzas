CC     = gcc
CFLAGS = -O2 -Wall -Istarter -Isrc
LIBS   = -lmatio -lm

COMMON = src/p_common.c starter/data.c

bin/seq: src/main_seq.c src/p_seq.c $(COMMON) src/tsne.h
	mkdir -p bin
	$(CC) $(CFLAGS) -o $@ src/main_seq.c src/p_seq.c $(COMMON) $(LIBS)

clean:
	rm -rf bin
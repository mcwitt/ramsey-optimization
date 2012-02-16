dsfmt_dir = $(HOME)/local/src/dSFMT-src-2.1
dsfmt_flags = -I$(dsfmt_dir) -DDSFMT_MEXP=2203 -DHAVE_SSE2

VPATH = $(dsfmt_dir)
CFLAGS = -Wall -std=gnu99 -O3 $(dsfmt_flags)
OBJECTS = dSFMT.o ramsey.o ramsey2.o

all: pt.out demon.out demon2.out

defs.h: gendefs.py
	python gendefs.py

ramsey.o ramsey2.o: defs.h

ramsey.o: ramsey.c ramsey.h

ramsey2.o: ramsey2.c ramsey2.h

pt.out demon.out demon2.out: dSFMT.o

pt.out: pt.c ramsey.o

demon.out: demon.c ramsey.o

demon2.out: demon2.c ramsey2.o
	$(CC) $^ $(CFLAGS) -DQUADRATIC -o $@

%.out: %.c
	$(CC) $^ $(CFLAGS) -o $@

clean:
	$(RM) $(OBJECTS)
	$(RM) defs.h

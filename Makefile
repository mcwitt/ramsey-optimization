dsfmt_dir = $(HOME)/local/src/dSFMT-src-2.1
dsfmt_flags = -I$(dsfmt_dir) -DDSFMT_MEXP=2203 -DHAVE_SSE2
sims = pt.out demon.out demon2.out demon2-2.out

VPATH = $(dsfmt_dir)
CFLAGS = -Wall -std=gnu99 -O3 $(dsfmt_flags)
LDFLAGS = -lm

all: $(sims) energy.out

%.out: %.o
	$(CC) $^ $(CFLAGS) $(LDFLAGS) -o $@

$(sims): dSFMT.o

defs.h: gendefs.py
	python gendefs.py

ramsey.o ramsey2.o: defs.h

pt.out: ramsey.o

demon.out demon2.out demon2-2.out: CFLAGS := $(CFLAGS) -DQUADRATIC
demon.out: ramsey.o
demon2.out: ramsey2.o
demon2-2.out: ramsey2.o

clean:
	$(RM) defs.h *.o

dsfmt_dir = $(HOME)/local/src/dSFMT-src-2.1
dsfmt_flags = -I$(dsfmt_dir) -DDSFMT_MEXP=2203 -DHAVE_SSE2
sims = pt.out sa.out demon.out demon2.out demon2-2.out genetic.out

VPATH = $(dsfmt_dir)
CFLAGS = -Wall -std=gnu99 -O3 $(dsfmt_flags) -g
LDFLAGS = -lm

all: $(sims) energy.out test.out

%.out: %.o
	$(CC) $^ $(CFLAGS) $(LDFLAGS) -o $@

$(sims): dSFMT.o

defs.h: gendefs.py
	python gendefs.py $R $S $(NV)

ramsey.o ramsey2.o: defs.h

pt.out: ramsey.o

sa.out: CFLAGS := $(CFLAGS) -DLITTLE_ENDIAN
sa.out: ramsey.o

demon2.out demon2-2.out: CFLAGS := $(CFLAGS) -DQUADRATIC
demon.out: ramsey.o
demon2.out: ramsey2.o
demon2-2.out: ramsey2.o

genetic.out: CFLAGS := $(CFLAGS) -DSGA_CHROMLEN=NED
genetic.out: ramsey.o sga.o

# TEST PROGRAMS
test.o: test.c
test.out: test.o ramsey.o dSFMT.o

gtest.out: sga.o dSFMT.o

clean:
	$(RM) defs.h *.o

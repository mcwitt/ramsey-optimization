dsfmt_dir = $(HOME)/local/src/dSFMT-src-2.1
dsfmt_flags = -I$(dsfmt_dir) -DDSFMT_MEXP=2203 -DHAVE_SSE2
sims = pt sa demon demon2 demon2-2 genetic eo

VPATH = $(dsfmt_dir)
CFLAGS = -Wall -std=gnu99 -O3 $(dsfmt_flags)
#CFLAGS = -Wall -std=gnu99 -O3 $(dsfmt_flags) -g
LDFLAGS = -lm

all: $(sims) energy

$(sims): dSFMT.o

defs.h: gendefs.py
	python gendefs.py $R $S $(NV)

ramsey.o ramsey2.o: defs.h

pt: ramsey.o

sa: ramsey.o

demon2 demon2-2: CFLAGS := $(CFLAGS) -DQUADRATIC
demon: ramsey.o
demon2: ramsey2.o
demon2-2: ramsey2.o

genetic: ramsey.o sga.o

eo: ramsey.o qselect.o

checks: ramsey.o dSFMT.o

.PHONY: run_checks
run_checks: checks
	./checks
	$(RM) checks

clean:
	$(RM) defs.h *.o

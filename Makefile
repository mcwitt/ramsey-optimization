dsfmt_dir = $(HOME)/local/src/dSFMT-src-2.1
dsfmt_flags = -I$(dsfmt_dir) -DDSFMT_MEXP=2203 -DHAVE_SSE2
need_dsfmt = pt sa demon demon2 demon2-2 genetic eo checks

VPATH = $(dsfmt_dir)
CFLAGS = -Wall -std=gnu99 -O3
#CFLAGS = -Wall -std=gnu99 -O3 -g
LDLIBS = -lm

all: pt sa demon genetic eo energy

defs.h: gendefs.py
	python gendefs.py $R $S $(NV)

$(need_dsfmt): CFLAGS += $(dsfmt_flags)
$(need_dsfmt): dSFMT.o

ramsey.o ramsey2.o: defs.h
pt sa demon genetic eo checks: ramsey.o
demon2 demon2-2: ramsey2.o
demon2 demon2-2: CFLAGS += -DQUADRATIC
genetic: sga.o
eo: qselect.o

.PHONY: run_checks
run_checks: checks
	./checks
	$(RM) checks

clean:
	$(RM) defs.h *.o

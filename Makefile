dsfmt_dir = $(HOME)/local/src/dSFMT-src-2.1
dsfmt_flags = -I$(dsfmt_dir) -DDSFMT_MEXP=2203 -DHAVE_SSE2

VPATH = $(dsfmt_dir)
#CFLAGS = -Wall -std=gnu99 -lm -O3 $(dsfmt_flags) -pg
CFLAGS = -Wall -std=gnu99 -lm -O3 $(dsfmt_flags)
OBJECTS = dSFMT.o

.PHONY: ramsey

ramsey_dbg: ramsey.c dSFMT.o
	gcc -g -DDEBUG -DNV=$(NV) -DS=$S -o $@_$S-$S-$(NV).out $^ $(CFLAGS)

ramsey: ramsey.c dSFMT.o
	gcc -g -DNV=$(NV) -DS=$S -o $@_$S-$S-$(NV).out $^ $(CFLAGS)

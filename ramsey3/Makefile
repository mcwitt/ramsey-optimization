dsfmt_dir = $(HOME)/local/src/dSFMT-src-2.1
dsfmt_flags = -I$(dsfmt_dir) -DDSFMT_MEXP=2203 -DHAVE_SSE2

VPATH = $(dsfmt_dir)
CFLAGS = -Wall -std=gnu99 -lm -O3 $(dsfmt_flags) -DNV=$(NV) -DS=$S -DNED=$(NED) -DNSG=$(NSG) -DNSGFE=$(NSGFE)
OBJECTS = dSFMT.o

.PHONY: ramsey

auto:
	python compile.py

ramsey: ramsey.c init_reps.c init_tabs.c io.c run.c sweep.c temper.c update_fields.c dSFMT.o
	$(CC) -o $@_$S-$S-$(NV).out $^ $(CFLAGS)

clean:
	$(RM) $(OBJECTS)
	$(RM) *.out

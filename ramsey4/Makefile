dsfmt_dir = $(HOME)/local/src/dSFMT-src-2.1
dsfmt_flags = -I$(dsfmt_dir) -DDSFMT_MEXP=2203 -DHAVE_SSE2

VPATH = $(dsfmt_dir)
CFLAGS = -Wall -std=gnu99 -lm -O3 $(dsfmt_flags) -DNV=$(NV) -DR=$R -DS=$S -DNED=$(NED) -DNSGR=$(NSGR) -DNSGS=$(NSGS) -DNSGFER=$(NSGFER) -DNSGFES=$(NSGFES)
OBJECTS = dSFMT.o

.PHONY: compile ramsey

compile:
	python compile.py $R $S $(NV)

ramsey: ramsey.c init_tabs.c init_reps.c io.c run.c sweep.c temper.c update_fields.c dSFMT.o
	$(CC) -o $@_$R-$S-$(NV).out $^ $(CFLAGS)

clean:
	$(RM) $(OBJECTS)
	$(RM) *.out

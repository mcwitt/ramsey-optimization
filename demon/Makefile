dsfmt_dir = $(HOME)/local/src/dSFMT-src-2.1
dsfmt_flags = -I$(dsfmt_dir) -DDSFMT_MEXP=2203 -DHAVE_SSE2

VPATH = ..:$(dsfmt_dir)
CFLAGS = -I.. -Wall -std=gnu99 -lm -O3 $(dsfmt_flags)
#CFLAGS = -I.. -Wall -std=gnu99 -lm -O3 $(dsfmt_flags) -DFULL_OUTPUT
#CFLAGS = -I.. -Wall -std=gnu99 -lm $(dsfmt_flags) -g	# FOR DEBUGGING
OBJECTS = dSFMT.o

.PHONY: demon

demon: demon.c ramsey.c dSFMT.o
	$(CC) -DNV=$(NV) -DR=$R -DS=$S -DNED=$(NED) -DNSGR=$(NSGR) -DNSGS=$(NSGS) -DNSGFER=$(NSGFER) -DNSGFES=$(NSGFES) -o $@_$R-$S-$(NV).out $^ $(CFLAGS)

clean:
	$(RM) $(OBJECTS)

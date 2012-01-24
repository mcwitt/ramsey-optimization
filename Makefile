dsfmt_dir = $(HOME)/local/src/dSFMT-src-2.1
dsfmt_flags = -I$(dsfmt_dir) -DDSFMT_MEXP=2203 -DHAVE_SSE2

VPATH = $(dsfmt_dir)
#CFLAGS = -Wall -std=gnu99 -lm -O3 $(dsfmt_flags) -pg
CFLAGS = -Wall -std=gnu99 -lm -O3 $(dsfmt_flags)
#CFLAGS = -Wall -std=gnu99 -lm $(dsfmt_flags)
OBJECTS = dSFMT.o

.PHONY: ramsey ramsey2 ramsey3 ramsey4

ramsey: ramsey.c dSFMT.o
	gcc -DNOTIME -DNV=$(NV) -DS=$S -DNED=$(NED) -DNSG=$(NSG) -DNSGFE=$(NSGFE) -o $@_$S-$S-$(NV).out $^ $(CFLAGS)

ramsey2: ramsey2.c dSFMT.o
	gcc -DNOTIME -DNV=$(NV) -DR=$R -DS=$S -DNED=$(NED) -DNSGR=$(NSGR) -DNSGS=$(NSGS) -DNSGFER=$(NSGFER) -DNSGFES=$(NSGFES) -o $@_$R-$S-$(NV).out $^ $(CFLAGS)

ramsey3: ramsey3.c dSFMT.o
	gcc -DNOTIME -DNV=$(NV) -DS=$S -DNED=$(NED) -DNSG=$(NSG) -DNSGFE=$(NSGFE) -o $@_$S-$S-$(NV).out $^ $(CFLAGS)

ramsey4: ramsey4.c dSFMT.o
	gcc -DNV=$(NV) -DR=$R -DS=$S -DNED=$(NED) -DNSGR=$(NSGR) -DNSGS=$(NSGS) -DNSGFER=$(NSGFER) -DNSGFES=$(NSGFES) -o $@_$R-$S-$(NV).out $^ $(CFLAGS)

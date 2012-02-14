dsfmt_dir = $(HOME)/local/src/dSFMT-src-2.1
dsfmt_flags = -I$(dsfmt_dir) -DDSFMT_MEXP=2203 -DHAVE_SSE2

VPATH = ..:$(dsfmt_dir)
#CFLAGS = -I.. -Wall -std=gnu99 -lm -O3 $(dsfmt_flags)
CFLAGS = -I.. -Wall -std=gnu99 -lm -O3 $(dsfmt_flags) -DFULL_OUTPUT
#CFLAGS = -I.. -Wall -std=gnu99 -lm $(dsfmt_flags) -DDEBUG		# FOR DEBUGGING
#CFLAGS = -I.. -Wall -std=gnu99 -lm $(dsfmt_flags) -DDEBUG -g	# FOR DEBUGGING W/ GDB
OBJECTS = dSFMT.o

.PHONY: demon
.PHONY: demon2
.PHONY: demon3

demon2: demon2.c ramsey2.c dSFMT.o
	$(CC) -DNV=$(NV) -DR=$R -DS=$S -DNED=$(NED) -DNSGR=$(NSGR) -DNSGS=$(NSGS) -DNSGFER=$(NSGFER) -DNSGFES=$(NSGFES) -o $@_$R-$S-$(NV).out $^ $(CFLAGS)

#demon3: demon3.c ramsey3.c dSFMT.o
#	$(CC) -DNV=$(NV) -DR=$R -DS=$S -DNED=$(NED) -DNSGR=$(NSGR) -DNSGS=$(NSGS) -DNSGFER=$(NSGFER) -DNSGFES=$(NSGFES) -o $@_$R-$S-$(NV).out $^ $(CFLAGS)

#demon: demon.c ramsey.c dSFMT.o
#	$(CC) -DNV=$(NV) -DR=$R -DS=$S -DNED=$(NED) -DNSGR=$(NSGR) -DNSGS=$(NSGS) -DNSGFER=$(NSGFER) -DNSGFES=$(NSGFES) -o $@_$R-$S-$(NV).out $^ $(CFLAGS)

clean:
	$(RM) $(OBJECTS)

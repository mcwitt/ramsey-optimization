dsfmt_dir = $(HOME)/local/src/dSFMT-src-2.1
dsfmt_flags = -I$(dsfmt_dir) -DDSFMT_MEXP=2203 -DHAVE_SSE2

VPATH = $(dsfmt_dir)
CFLAGS = -Wall -std=gnu99 -lm -O3 $(dsfmt_flags) -DNV=$(NV) -DS=$S -DNED=$(NED) -DNOTIME
OBJECTS = dSFMT.o

.PHONY: ramsey ramsey2 ramsey3 ramsey4 test

ramsey: ramsey.c dSFMT.o
	$(CC) -DNSG=$(NSG) -DNSGFE=$(NSGFE) -o $@_$S-$S-$(NV).out $^ $(CFLAGS)

ramsey2: ramsey2.c dSFMT.o
	$(CC) -DR=$R -DNSGR=$(NSGR) -DNSGS=$(NSGS) -DNSGFER=$(NSGFER) -DNSGFES=$(NSGFES) -o $@_$R-$S-$(NV).out $^ $(CFLAGS)

ramsey3: ramsey3.c dSFMT.o
	$(CC) -DNSG=$(NSG) -DNSGFE=$(NSGFE) -o $@_$S-$S-$(NV).out $^ $(CFLAGS)

ramsey4: ramsey4.c dSFMT.o
	$(CC) -DR=$R -DNED=$(NED) -DNSGR=$(NSGR) -DNSGS=$(NSGS) -DNSGFER=$(NSGFER) -DNSGFES=$(NSGFES) -o $@_$R-$S-$(NV).out $^ $(CFLAGS)

test:
	python compile.py  		5 	35
	python compile2.py		5 5	35
	python compile3.py		5 	35
	python compile4.py		5 5	35
	time ./ramsey_5-5-35.out 	temps.txt 100 101 123 > out1.tmp
	time ./ramsey2_5-5-35.out	temps.txt 100 101 123 > out2.tmp
	time ./ramsey3_5-5-35.out	temps.txt 100 123 > out3.tmp
	time ./ramsey4_5-5-35.out	temps.txt 100 123 > out4.tmp
	diff out1.tmp out2.tmp
	diff out1.tmp out3.tmp
	diff out1.tmp out4.tmp
	$(RM) ramsey*_5-5-35.out *.tmp
	

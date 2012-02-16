VPATH = ..:$(dsfmt_dir)
CFLAGS = -I.. -Wall -std=gnu99 -lm -O3 $(dsfmt_flags)
#CFLAGS = -I.. -Wall -std=gnu99 -lm $(dsfmt_flags) -g	# FOR DEBUGGING W/ GDB
OBJECTS = dSFMT.o

dsfmt_dir = $(HOME)/local/src/dSFMT-src-2.1
dsfmt_flags = -I$(dsfmt_dir) -DDSFMT_MEXP=2203 -DHAVE_SSE2

pt.out: pt.c ramsey.c dSFMT.o
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	$(RM) $(OBJECTS)

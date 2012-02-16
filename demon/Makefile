VPATH = ..:$(dsfmt_dir)
CFLAGS = -I.. -Wall -std=gnu99 -lm -O3 $(dsfmt_flags)
#CFLAGS = -I.. -Wall -std=gnu99 -lm $(dsfmt_flags) -DDEBUG		# FOR DEBUGGING
#CFLAGS = -I.. -Wall -std=gnu99 -lm $(dsfmt_flags) -DDEBUG -g	# FOR DEBUGGING W/ GDB
OBJECTS = dSFMT.o

dsfmt_dir = $(HOME)/local/src/dSFMT-src-2.1
dsfmt_flags = -I$(dsfmt_dir) -DDSFMT_MEXP=2203 -DHAVE_SSE2

demon.out: demon.c ramsey.c dSFMT.o
	$(CC) -o $@ $^ $(CFLAGS)

demon2.out: demon2.c ramsey2.c dSFMT.o
	$(CC) -DQUADRATIC -o $@ $^ $(CFLAGS)

clean:
	$(RM) $(OBJECTS)

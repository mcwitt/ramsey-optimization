CFLAGS = -Wall -std=gnu99 -lm -O3
#CFLAGS = -Wall -std=gnu99 -lm -g

%.out: %.c
	$(CC) $(CFLAGS) $^ -o $@

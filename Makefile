
CC=gcc
CFLAGS=-O3
#DEFINES=-DDEBUG
DEFINES=-DPOLAR_SELF

.c.o:
	$(CC) -c $(CFLAGS) $(DEFINES) -I. $*.c

all:	main.o cleanup.o input.o pairs.o pbc.o surface.o energy.o polar.o
	$(CC) $(CFLAGS) $(DEFINES) *.o -o 6dpes

clean:
	rm -f *.o 6dpes



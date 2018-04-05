

all: mcf 

mcf: implicit.c  mcfutil.c  pbeampp.c  pflowup.c   pstart.c   treeup.c mcf.c       output.c   pbla.c     psimplex.c  readmin.c
	gcc -o mcf *.c 
clean:
	rm -f mcf *.o


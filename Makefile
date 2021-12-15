all: nano

clean:
	/bin/rm -f *.o

PHG_MAKEFILE_INC = /soft/apps/phg/intel-2017u4/impi-2017u3/phg-0.9.6/share/phg/Makefile.inc
include ${PHG_MAKEFILE_INC}


nano.o: PNP_coefficient.h nano.c PNP_analytic.c 


runmp:
	mpirun -np 4  ./nano

BASE = $(COSMOS_HOME)
CLARDY = $(BASE)/clardy

# CFLAGS = -O
CC =   gcc  $(CFLAGS) -I$(CLARDY)/include -DOPTUTILH="\"optutils.h\""
         -DFEATSH="\"mgfeats.h\""\  -DIOUTILH="\"ioutils.h\"" -DKDCUTILH="\"kdcutil.h\"" \
			-DMASKDATH="\"maskdat.h\"" -DMGUTILH="\"mgutils.h\""

libimacs:
	${CC} -c  ioutils.c kdcutil.c mgutils.c optutils.c sincos.c
	${AR} r libimacs.a   ioutils.o kdcutil.o mgutils.o  optutils.o sincos.o
	ranlib libimacs.a
	mv libimacs.a $(BASE)/lib
	rm *.o

libimacs_x:
	${CC} -c  -fPIC ioutils.c kdcutil.c mgutils.c optutils.c sincos.c
	${AR} r libimacs_x.a   ioutils.o kdcutil.o mgutils.o  optutils.o sincos.o
	ranlib libimacs_x.a
	mv libimacs_x.a  $(BASE)/source
	rm *.o

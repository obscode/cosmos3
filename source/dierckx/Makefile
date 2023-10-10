BASE = $(COSMOS_HOME)

libspline.a:
	gfortran -O -c  *.f
	${AR} r libspline.a  *.o
	ranlib libspline.a
	mv libspline.a $(BASE)/lib
	rm *.o


INSTALLPATH=../../
INCLUDE=pil.h pil_error.h



all : libpil.a
	cp ${INCLUDE} ${INSTALLPATH}/include

libpil.a : pil.o pil_error.o pil_sys.o
	ar cr libpil.a pil.o pil_error.o pil_sys.o
	cp libpil.a ${INSTALLPATH}/lib

pil.o : pil.c ${INCLUDE}
	gcc -c pil.c $(INCLUDE)

pil_error.o : pil_error.c ${INCLUDE}
	gcc -c pil_error.c ${INCLUDE}

pil_sys.o : pil_sys.c
	gcc -c pil_sys.c ${INCLUDE}

clean:
	rm pil.o pil_error.o pil_sys.o libpil.a

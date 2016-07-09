
A1.o: A1.FOR
	gfortran -c A1.FOR

B1.o: B1.FOR
	gfortran -c B1.FOR

C1.o: C1.FOR
	gfortran -c C1.FOR

MAIN.o: MAIN.FOR
	gfortran -c MAIN.FOR

build: A1.o B1.o C1.o MAIN.o
	gfortran A1.o B1.o C1.o MAIN.o -o SHAKE91

install: build
	cp SHAKE91 /usr/local/bin/SHAKE91

clean:
	rm -f *.o
	rm -f SHAKE91

uninstall: clean
	rm -f /usr/local/bin/SHAKE91

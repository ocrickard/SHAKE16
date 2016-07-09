
build: A1.o B1.o C1.o Main.o
	gfortran A1.o B1.o C1.o Main.o -o SHAKE91

A1.o: A1.f
	gfortran -c A1.f

B1.o: B1.f
	gfortran -c B1.f

C1.o: C1.f
	gfortran -c C1.f

Main.o: Main.f
	gfortran -c Main.f

install: build
	cp SHAKE91 /usr/local/bin/SHAKE91

clean:
	rm -f *.o
	rm -f SHAKE91

uninstall: clean
	rm -f /usr/local/bin/SHAKE91

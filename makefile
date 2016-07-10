
build: output/build/A1.o output/build/B1.o output/build/C1.o output/build/Main.o
	gfortran output/build/A1.o output/build/B1.o output/build/C1.o output/build/Main.o -o SHAKE16

output/build/A1.o: A1.f
	gfortran -c A1.f -o output/build/A1.o

output/build/B1.o: B1.f
	gfortran -c B1.f -o output/build/B1.o

output/build/C1.o: C1.f
	gfortran -c C1.f -o output/build/C1.o

output/build/Main.o: Main.f
	gfortran -c Main.f -o output/build/Main.o

install: build
	cp SHAKE16 /usr/local/bin/SHAKE16

clean:
	rm -rf output/
	rm -f SHAKE16

uninstall: clean
	rm -f /usr/local/bin/SHAKE16

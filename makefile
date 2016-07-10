
build: output/build/A1.o output/build/B1.o output/build/C1.o output/build/Main.o
	gfortran output/build/A1.o output/build/B1.o output/build/C1.o output/build/Main.o -o SHAKE16

output/build/:
	mkdir -p output/build

output/build/A1.o: A1.f output/build/
	gfortran -c A1.f -o output/build/A1.o

output/build/B1.o: B1.f output/build/
	gfortran -c B1.f -o output/build/B1.o

output/build/C1.o: C1.f output/build/
	gfortran -c C1.f -o output/build/C1.o

output/build/Main.o: Main.f output/build/
	gfortran -c Main.f -o output/build/Main.o

install: build
	cp SHAKE16 /usr/local/bin/SHAKE16

test: install
	sh test.sh

clean:
	rm -rf output/
	rm -f SHAKE16

uninstall: clean
	rm -f /usr/local/bin/SHAKE16

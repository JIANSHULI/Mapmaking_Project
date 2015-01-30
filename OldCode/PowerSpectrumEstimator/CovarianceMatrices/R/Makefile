LD=g++
LDFLAGS=-L/usr/lib64 -lfftw3 -g -L/usr/lib64 -lgsl -lgslcblas -L/nfs/blank/h4213/jsdillon/Libraries04/FFTW/lib -L/nfs/blank/h4213/jsdillon/Libraries04/GSL/lib 
CXXFLAGS= -g
IFLAGS = -I/nfs/blank/h4213/jsdillon/Libraries04/GSL/include/ -I/nfs/blank/h4213/jsdillon/Libraries04/FFTW/include/

R : R.o ../../CommonClasses/Specs.o
	$(LD) -o $@ $^ $(LDFLAGS)

R.o : R.cpp  ../../CommonClasses/Specs.h 
	$(LD) $(IFLAGS) -c $^ $(LDFLAGS)

../../CommonClasses/Specs.o: ../../CommonClasses/Specs.cpp ../../CommonClasses/Specs.h

.dummy clean :
	rm -f R *.o
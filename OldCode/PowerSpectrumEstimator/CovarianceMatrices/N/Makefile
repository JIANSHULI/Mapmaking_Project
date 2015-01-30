LD=g++
LDFLAGS=-L/usr/lib64 -lfftw3 -g -L/usr/lib64 -lgslcblas -L/nfs/blank/h4213/jsdillon/Libraries/FFTW/lib -L/nfs/blank/h4213/jsdillon/Libraries/GSL/lib 
CXXFLAGS= -g
IFLAGS = -I/nfs/blank/h4213/jsdillon/Libraries/GSL/include/ -I/nfs/blank/h4213/jsdillon/Libraries/FFTW/include/

N : N.o ../../CommonClasses/Specs.o ../../CommonClasses/CVector.o 
	$(LD) -o $@ $^ $(LDFLAGS)

N.o : N.cpp  ../../CommonClasses/Specs.h ../../CommonClasses/CVector.h
	$(LD) $(IFLAGS) -c $^ $(LDFLAGS)

../../CommonClasses/CVector.o : ../../CommonClasses/CVector.cpp ../../CommonClasses/CVector.h ../../CommonClasses/Specs.h
	$(LD) $(IFLAGS) -c $^ $(LDFLAGS)

../../CommonClasses/Specs.o: ../../CommonClasses/Specs.cpp ../../CommonClasses/Specs.h

.dummy clean :
	rm -f N *.o
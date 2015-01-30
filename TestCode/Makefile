LD=g++
LDFLAGS= -L$(FFTW_LIBRARY_PATH) -lfftw3 -L$(GSL_LIBRARY_PATH) -lgsl -lgslcblas
CXXFLAGS= -g
IFLAGS = -I$(FFTW_HEADER_PATH) -I$(GSL_HEADER_PATH)

test : test.o 
	$(LD) $(IFLAGS) -o  $@ $^ $(LDFLAGS)

test.o : test.cpp 
	$(LD) -c $(IFLAGS) $^ 

.dummy clean :
	rm -f test *.o
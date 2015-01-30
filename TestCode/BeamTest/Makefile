LD=g++
LDFLAGS= -L$(FFTW_LIBRARY_PATH) -lfftw3 -L$(GSL_LIBRARY_PATH) -lgsl -lgslcblas
CXXFLAGS= -g
IFLAGS = -I$(FFTW_HEADER_PATH) -I$(GSL_HEADER_PATH)

beamTest : beamTest.o 
	$(LD) $(IFLAGS) -o  $@ $^ $(LDFLAGS)

beamTest.o : beamTest.cpp 
	$(LD) -c $(IFLAGS) $^ 

.dummy clean :
	rm -f beamTest *.o
LD=g++
LDFLAGS= -L$(FFTW_LIBRARY_PATH) -lfftw3
CXXFLAGS= -g
IFLAGS = -I$(FFTW_HEADER_PATH)

facetedMapmaker : facetedMapmaker.o 
	$(LD) $(CXXFLAGS) $(IFLAGS) -o  $@ $^ $(LDFLAGS)

facetedMapmaker.o : facetedMapmaker.cpp 
	$(LD) $(CXXFLAGS) -c $(IFLAGS) $^ 

.dummy clean :
	rm -f facetedMapmaker *.o
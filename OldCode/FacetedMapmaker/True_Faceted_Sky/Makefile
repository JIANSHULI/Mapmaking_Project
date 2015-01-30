LD=g++
LDFLAGS= -L$(FFTW_LIBRARY_PATH) -lfftw3 -L$(GSL_LIBRARY_PATH) -lgsl -lgslcblas  -L$(HEALPIX_LIBRARY_PATH) -lhealpix_cxx -lc_utils -lcxxsupport -lfftpack -lsharp -L$(CFITSIO_LIBRARY_PATH) -lcfitsio
CXXFLAGS= -g 
IFLAGS = -I$(FFTW_HEADER_PATH) -I$(GSL_HEADER_PATH) -I$(HEALPIX_HEADER_PATH) -I$(CFITSIO_HEADER_PATH) -I/opt/local/include -I/Users/jsdillon/Desktop/Joint_Mapmaking_Power_Spectrum_Pipeline/Healpix_3.11/src/cxx/cxxsupport

True_Faceted_Sky : True_Faceted_Sky.o
	$(LD) $(IFLAGS) -o  $@ $^ $(LDFLAGS) -O3

True_Faceted_Sky.o : True_Faceted_Sky.cpp
	$(LD) -c $(IFLAGS) $^

.dummy clean :
	rm -f True_Faceted_Sky *.o

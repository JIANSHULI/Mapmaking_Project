LD=g++
LDFLAGS= -L$(FFTW_LIBRARY_PATH) -lfftw3 -L$(GSL_LIBRARY_PATH) -lgsl -lgslcblas  -L$(HEALPIX_LIBRARY_PATH) -lhealpix_cxx -lc_utils -lcxxsupport -lfftpack -lsharp -L$(CFITSIO_LIBRARY_PATH) -lcfitsio
CXXFLAGS= -g 
IFLAGS = -I$(FFTW_HEADER_PATH) -I$(GSL_HEADER_PATH) -I$(HEALPIX_HEADER_PATH) -I$(CFITSIO_HEADER_PATH) -I/opt/local/include -I/Users/jsdillon/Desktop/Joint_Mapmaking_Power_Spectrum_Pipeline/Healpix_3.11/src/cxx/cxxsupport

GSMVisibilitiesAllFreqAllPol : GSMVisibilitiesAllFreqAllPol.o
	$(LD) $(CXXFLAGS) $(IFLAGS) -o  $@ $^ $(LDFLAGS)

GSMVisibilitiesAllFreqAllPol.o : GSMVisibilitiesAllFreqAllPol.cpp
	$(LD) $(CXXFLAGS) -c $(IFLAGS) $^

.dummy clean :
	rm -f GSMVisibilitiesAllFreqAllPol *.o

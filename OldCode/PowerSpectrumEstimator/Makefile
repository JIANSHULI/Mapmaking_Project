LD=g++
LDFLAGS= -L$(FFTW_LIBRARY_PATH) -lfftw3 -L$(GSL_LIBRARY_PATH) -lgsl -lgslcblas
CXXFLAGS= -g
IFLAGS = -I$(FFTW_HEADER_PATH) -I$(GSL_HEADER_PATH)

powerSpectrumEstimator : powerSpectrumEstimator.o ../CommonClasses/CVector.o ../CommonClasses/RDMatrix.o ../CommonClasses/Toeplitz.o ../CommonClasses/Specs.o
	$(LD) $(IFLAGS) -o  $@ $^ $(LDFLAGS)

powerSpectrumEstimator.o : powerSpectrumEstimator.cpp ../CommonClasses/CVector.hpp ../CommonClasses/RDMatrix.hpp ../CommonClasses/Toeplitz.hpp ../CommonClasses/Specs.hpp
	$(LD) -c $(IFLAGS) $^ 

../CommonClasses/RDMatrix.o : ../CommonClasses/RDMatrix.cpp ../CommonClasses/RDMatrix.hpp ../CommonClasses/CVector.hpp ../CommonClasses/Specs.hpp
	$(LD) -c $(IFLAGS) $^
	mv ./RDMatrix.o ../CommonClasses/RDMatrix.o 

../CommonClasses/Toeplitz.o : ../CommonClasses/Toeplitz.cpp ../CommonClasses/Toeplitz.hpp ../CommonClasses/CVector.hpp ../CommonClasses/Specs.hpp
	$(LD) -c $(IFLAGS) $^
	mv ./Toeplitz.o ../CommonClasses/Toeplitz.o 

../CommonClasses/CVector.o : ../CommonClasses/CVector.cpp ../CommonClasses/CVector.hpp ../CommonClasses/Specs.hpp
	$(LD) -c $(IFLAGS) $^
	mv ./CVector.o ../CommonClasses/CVector.o 

../CommonClasses/Specs.o : ../CommonClasses/Specs.cpp ../CommonClasses/Specs.hpp
	$(LD) -c $(IFLAGS) $^
	mv ./Specs.o ../CommonClasses/Specs.o 

.dummy clean :
	rm -f powerSpectrumEstimator *.o
	rm -f ../CommonClasses/*.o
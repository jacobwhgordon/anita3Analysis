## This Makefile compiles the interferometer and interferometer accessories

LIB_FLAGS = -L. -L$(ROOTSYS)/lib -L$(ANITA_UTIL_INSTALL_DIR)/lib -L$(HEALPIX_PATH)/lib -L${ROOTSYS}/lib
INCLUDE_FLAGS = -I. -I$(ROOTSYS)/include -I$(ANITA_UTIL_INSTALL_DIR)/include -I$(HEALPIX_PATH)/include -I${MINUIT2_PATH}inc ${MPI_CXXFLAGS}
CXX_FLAGS = -fPIC
LIBS = -lAnitaEvent -lAnitaCorrelator -lRootFftwWrapper -lgsl -lgslcblas -lm -lfftw3 -pthread -lm -ldl -lCore -lRIO -lNet \
  -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore \
  -lThread -pthread -lm -ldl -rdynamic -lfftw3 -lMathMore -lMinuit -lMinuit2 `root-config --cflags --glibs` \
  -lAnitaAnalysis -lUCorrelator \
  -lhealpix_cxx -lcxxsupport -lsharp -lfftpack -lc_utils \
  ${MPI_LIBS}
##  -lhealpix_cxx -lcxxsupport -lsharp -lfftpack -lc_utils -lcfitsio
## compiler verbose and warning flags
_V =
_W = -Wall -Wno-sign-compare
##_G = -g -O0
##_G = 
_G = -g
_S = -std=c++11


##fitFunc: fitFunc.cxx fitFunc.h
##	g++ $(_V) $(_W) $(_G) $(_S) -O0 -o $@.o $@.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS)

runInterferometer : runInterferometer.cxx interferometryUtil.h 
	g++ $(_V) $(_W) $(_G) $(_S) -o runInterferometer runInterferometer.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) -lInterferometer -fopenmp

analyzerIterator06 : analyzerIterator06.cxx
	g++ $(_V) $(_W) $(_G) $(_S) -o analyzerIterator06 analyzerIterator06.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) 

analyzerResultsIterator06 : analyzerResultsIterator06.cxx
	g++ $(_V) $(_W) $(_G) $(_S) -o analyzerResultsIterator06 analyzerResultsIterator06.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) 

runAnitaAnalyzer : runAnitaAnalyzer.cxx libInterferometer.so 
	g++ $(_V) $(_W) $(_G) $(_S) -o runAnitaAnalyzer runAnitaAnalyzer.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) -lInterferometer -fopenmp -lAnitaAnalysis

errorCorrelationAnalysis : errorCorrelationAnalysis.cxx 
	g++ $(_V) $(_W) $(_G) $(_S) -o errorCorrelationAnalysis errorCorrelationAnalysis.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) 

##runInterferometer2 : runInterferometer2.cxx interferometryUtil.h libInterferometer.so
##	g++ $(_V) $(_W) $(_G) $(_S) -o runInterferometer2 runInterferometer2.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) -lInterferometer

##runInterfSingleEvent : runInterfSingleEvent.cxx interferometryUtil.h libInterferometer.so
##	g++ $(_V) $(_W) $(_G) $(_S) -o runInterfSingleEvent runInterfSingleEvent.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) -lInterferometer 

templateBuilder : templateBuilder.cxx libPointingResult.so
	g++ $(_V) $(_W) $(_G) $(_S) -o templateBuilder templateBuilder.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) -lInterferometer -lPointingResult

buildBgProfile : buildBgProfile.cxx libInterferometer.so
	g++ $(_V) $(_W) $(_G) $(_S) -o buildBgProfile buildBgProfile.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) -lInterferometer

buildPreInterfFeatures : buildPreInterfFeatures.cxx libInterferometer.so
	g++ $(_V) $(_W) $(_G) $(_S) -o buildPreInterfFeatures buildPreInterfFeatures.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) -lInterferometer

showPowerSpectrum : showPowerSpectrum.cxx 
	g++ $(_V) $(_W) $(_G) $(_S) -o showPowerSpectrum showPowerSpectrum.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) 

pulseBrain : pulseBrain.cxx 
	g++ $(_V) $(_W) $(_G) $(_S) -o pulseBrain pulseBrain.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) 

##polyWog : polyWog.cxx 
##	g++ $(_V) $(_W) $(_G) $(_S) -o polyWog polyWog.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) 

##convertPointingResultsToRoot : convertPointingResultsToRoot.cxx 
##	g++ $(_V) $(_W) $(_G) $(_S) -o convertPointingResultsToRoot convertPointingResultsToRoot.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) 

##libEventCutVariables.so : EventCutVariables.cxx EventCutVariables.h 
##	g++ $(_V) $(_W) $(_G) $(_S) -c EventCutVariables.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) -fopenmp
##	g++ $(_V) $(_W) $(_S) -shared -o libEventCutVariables.so EventCutVariables.o $(LIB_FLAGS) $(LIBS)  -fopenmp

libInterferometer.so : Interferometer.cxx Interferometer.h PointingResult.h libbwFilter.so
	g++ $(_V) $(_W) $(_G) $(_S) -c Interferometer.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) -lbwFilter -fopenmp
	g++ $(_V) $(_W) $(_S) -shared -o libInterferometer.so Interferometer.o $(LIB_FLAGS) $(LIBS) -lbwFilter -fopenmp

libAnalysisResult.so : AnalysisResult.h 
	g++ $(_V) $(_W) $(_G) $(_S) -c AnalysisResult.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) -lbwFilter -fopenmp
	g++ $(_V) $(_W) $(_S) -shared -o libAnalysisResult.so AnalysisResult.o $(LIB_FLAGS) $(LIBS)

libbwFilter.so : bwFilter.cxx bwFilter.h
	g++ $(_V) $(_W) $(_G) $(_S) -c bwFilter.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS)
	g++ $(_V) $(_W) $(_S) -shared -o libbwFilter.so bwFilter.o $(LIB_FLAGS) $(LIBS) 

libPointingResult.so : PointingResult.cxx PointingResult.h
	g++ $(_V) $(_W) $(_G) $(_S) -c PointingResult.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS)
	g++ $(_V) $(_W) $(_S) -shared -o libPointingResult.so PointingResult.o $(LIB_FLAGS) $(LIBS) 

pointingResultsViewer : pointingResultsViewer.cxx libInterferometer.so
	g++ $(_V) $(_W) $(_G) $(_S) -o pointingResultsViewer pointingResultsViewer.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) -lInterferometer

eventLocalizer : eventLocalizer.cxx libInterferometer.so 
	g++ $(_V) $(_W) $(_G) $(_S) -o eventLocalizer eventLocalizer.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) -lInterferometer



#buildCoastlineGraph : buildCoastlineGraph.cxx libInterferometer.so 
#	g++ $(_V) $(_W) $(_G) $(_S) -o buildCoastlineGraph buildCoastlineGraph.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) -lInterferometer

#velocityProfile : velocityProfile.cxx
#	g++ $(_V) $(_W) $(_S) -o velocityProfile velocityProfile.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS)

#pulseTemplate : pulseTemplate.cxx
#	g++ $(_V) $(_W) $(_S) -o pulseTemplate pulseTemplate.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS)

#buildCoastlineMap : buildCoastlineMap.cxx
#	g++ $(_V) $(_W) $(_S) -o buildCoastlineMap buildCoastlineMap.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) 
	
%: %.cxx
	g++ $(_V) $(_W) $(_G) $(_S) -o $@ $@.cxx fitFunc.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS) 
	
#graphtest : graphtest.cxx 
#	g++ $(_V) $(_W) -o graphtest graphtest.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS)

#vector3Dtest : vector3Dtest.cxx 
#	g++ $(_V) $(_W) -o vector3Dtest vector3Dtest.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS)

#interferometryUtil : interferometryUtil.cxx 
#	g++ $(_V) $(_W) -o interferometryUtil interferometryUtil.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS)

#interferometryUtil : interferometryUtil.cxx 
#	g++ $(_V) $(_W) -o interferometryUtil interferometryUtil.cxx $(CXX_FLAGS) $(INCLUDE_FLAGS) $(LIB_FLAGS) $(LIBS)

##interferometerDict.C:
##	rootcint interferometerDict.C -c $(ANITA_UTIL_INSTALL_DIR)/include/CalibratedAnitaEvent.h

clean : 
	@rm -f interferometer.o 
	@rm -f interferometer
	@rm -f vector3Dtest.o 
	@rm -f vector3Dtest
	@rm -f velocityProfile.o 
	@rm -f velocityProfile
	@rm -f graphtest.o 
	@rm -f graphtest
	@rm -f pulseTemplate.o 
	@rm -f pulseTemplate
	@rm -f interferometryUtil.o 
	@rm -f interferometryUtil


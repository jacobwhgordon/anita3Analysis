// this program:
//   instantiates an Analyzer and supporting stuff,
//   iterates through events of a run
//     for both circular and linear polarizations
//       calculate quality cuts and save results in an array
//       

// TODO insert cut logic
//   TODO accept flags to activate cuts
//   TODO calculate cuts (copy logic from qualityCuts.cxx and analyzerResultsViewer.cxx)
//   TODO reckon on separation of quality and analysis cuts (just use two separate arrays?)
//   TODO save the cut results to the 

#include "AnitaConventions.h"
#include "Analyzer.h"
#include "AnalysisConfig.h"
#include "AnitaEventSummary.h"
#include "AnalysisWaveform.h"
#include "CalibratedAnitaEvent.h"
#include "RawAnitaHeader.h"
#include "UsefulAnitaEvent.h"
#include "Adu5Pat.h"
#include "UsefulAdu5Pat.h"
#include "FilteredAnitaEvent.h"
//#include "SineSubtract.h"
#include "FilterStrategy.h"
#include "BasicFilters.h"
#include "GeomFilter.h"
#include "AnitaGeomTool.h"
#include "UCFilters.h"
#include "TruthAnitaEvent.h"
#include "UCUtil.h"

#include "TCanvas.h"
#include "TH2.h"
#include "TApplication.h"
#include "TFile.h"
#include "TMarker.h"
#include "TChain.h"
#include "TStyle.h"
#include "TRandom.h"
#include "InterfUtil.h"
#include "analysisCuts.h"

using namespace std;
using namespace UCorrelator;

int processKeywordParm(char* kwp);

bool useSineSubFiltering = true;
bool useAbbysFiltering = false;
int filterOption = FILTER_OPTION_NONE;
int adaptiveNotchFilterOrder = 5;
float adaptiveNotchFilterWidth = 0.25;
bool adaptiveNotchFilterCausal = true;

int main(int argc, char* argv[]) {
  cout << "Hello world " << argv[0] << endl;

  int eventNumber = 0;
  int runNumber = 345;
  int pNum = 0;
  int maxEvents = 10000000;
  bool doGraphics = false;
  bool fileOutput = false;
  bool knownSource = false;
  bool useTrigNsWindow = false;
  int dbStride = 10;
  char name[128];
  //char title[128];
  //vector<int> pulserNsTimes(0, 0);
  //double ksLon = AnitaLocations::LONGITUDE_WAIS;
  //double ksLat = AnitaLocations::LATITUDE_WAIS;
  //double ksAlt = AnitaLocations::ALTITUDE_WAIS;
  char calPulser[16] = "NONE";
  bool doCircPol = true;
  bool doLinPol = false;
  char* outputDirStr = getenv("ANITA_ANALYSIS_RESULTS_DIR");
  printf("output dir is %s \n", outputDirStr);
  char outputDir[1024];
  strcpy(outputDir, outputDirStr);
  //bool useFilter = false;
  bool simulationMode = false;
  bool getBaselineFromFile = false;
  int startEntryNum = 0;
  char baselineFilename[1024] = "baselineSample.root";
  AnalysisConfig::NormalizationOption_t normalizationType = AnalysisConfig::NormalizationOverlap;
  // process command-line parms
  char cmdLine[4096] = "";
  strcat(cmdLine, argv[0]);
  strcat(cmdLine, " ");  
  for (int i=1; i<argc; ++i) {
    strcat(cmdLine, argv[i]);
    strcat(cmdLine, " ");
    if (argv[i][0]=='-') {
      if (strlen(argv[i]) > 1) {
        char flag=argv[i][1];
        //switch (argv[i][1]) {
        if      (flag=='-')                 processKeywordParm(argv[i]);        
        else if (flag=='V' || flag=='v')    verboseMode = true; 
        else if (flag=='G' || flag=='g')    doGraphics = true; 
        else if (flag=='K' || flag=='k')    knownSource = true;
        //else if (flag=='F' || flag=='f')    useFilter = true;       // cal pulses
        else if (flag=='U' || flag=='u')    simulationMode = true; 
        else if (flag=='P' || flag=='p') {
          if (strlen(argv[i])>2) {
            if (argv[i][2]=='C' || argv[i][2]=='c') {
              doCircPol = true;
              doLinPol = false;
            }
            else if (argv[i][2]=='L' || argv[i][2]=='l') {
              doCircPol = false;
              doLinPol = true;
            }
            else if (argv[i][2]=='B' || argv[i][2]=='b') {
              doCircPol = true;
              doLinPol = true;
            }
          }
        }    
        else if (flag=='T' || flag=='t') {
          useTrigNsWindow = true;       // cal pulses
          sprintf(calPulser, "%s", "WAIS");
          if (strlen(argv[i])>2) {
            char thisParm[128];
            for (unsigned int k=2; k<=strlen(argv[i]); k++) {
              thisParm[k-2] = toupper(argv[i][k]);
            }
            sprintf(calPulser, "%s", thisParm);
          }
        }
        else if (flag=='B')    {
          getBaselineFromFile = true; 
          if (strlen(argv[i])>2) {
            char thisParm[1024];
            for (unsigned int k=2; k<=strlen(argv[i]); k++) {
              thisParm[k-2] = argv[i][k];
            }
            sprintf(baselineFilename, "%s", thisParm);
          }
        }    
        else if (flag=='N' || flag=='n') { 
          if (strlen(argv[i]) > 2) {
            char thisParm[8];
            for (unsigned int k=2; k<=strlen(argv[i]); k++) {
              thisParm[k-2] = argv[i][k];
            }
            if (strcmp(thisParm, "overlap")==0) {
              normalizationType = AnalysisConfig::NormalizationOverlap;
            } else if (strcmp(thisParm, "std")==0) {
              normalizationType = AnalysisConfig::NormalizationStandard;
            } else if (strcmp(thisParm, "none")==0) {
              normalizationType = AnalysisConfig::NormalizationNone;
            }
          }     
        } else if (flag=='S') { 
          if (strlen(argv[i]) > 2) {
            char thisParm[8];
            for (unsigned int k=2; k<=strlen(argv[i]); k++) {
              thisParm[k-2] = argv[i][k];
            }
            dbStride = atoi(thisParm);
          }            
        } else if (flag=='s') { 
          if (strlen(argv[i]) > 2) {
            char thisParm[8];
            for (unsigned int k=2; k<=strlen(argv[i]); k++) {
              thisParm[k-2] = argv[i][k];
            }
            startEntryNum = atoi(thisParm);
          }            
        } else if (flag=='M' || flag=='m') {
          if (strlen(argv[i]) > 2) {
            char thisParm[8];
            for (unsigned int k=2; k<=strlen(argv[i]); k++) {
              thisParm[k-2] = argv[i][k];
            }
            maxEvents = atoi(thisParm);
          }
        } else if (flag=='O' || flag=='o') {
          fileOutput = true;
          if (strlen(argv[i]) > 2) {
            unsigned int k=2;
            for (; k<=strlen(argv[i]); k++) {
              outputDir[k-2] = argv[i][k];
            }
            outputDir[k-2] = 0;
          }
        }
      }
    } else {
      ++pNum;
      if (pNum==1) runNumber = atoi(argv[i]); 
      else if (pNum==2) {
        eventNumber = atoi(argv[i]); 
        maxEvents = 1;
      }
    }
  }

/*   // DEPRECATED
  if (strcmp(calPulser, "WAIS")==0) {
    ksLon = AnitaLocations::LONGITUDE_WAIS;
    ksLat = AnitaLocations::LATITUDE_WAIS;
    ksAlt = AnitaLocations::ALTITUDE_WAIS;
    //pulserNsTimes.push_back(pulserNsWais);
  } 
  else if (strcmp(calPulser, "LDB")==0) {
    ksLon = AnitaLocations::LONGITUDE_LDB;
    ksLat = AnitaLocations::LATITUDE_LDB;
    ksAlt = AnitaLocations::ALTITUDE_LDB;
    //pulserNsTimes.push_back(pulserNsLdb1);  
    //pulserNsTimes.push_back(pulserNsLdb2);  
    //pulserNsTimes.push_back(pulserNsLdb3); 
  } // TODO other LDB timings
*/
  TApplication* app = new TApplication("app", &argc, argv);  
  if (!doGraphics) gROOT->SetBatch(kTRUE);
  printf("command line input is:\n"); 
  printf("%s \n", cmdLine);
  cout << "run number " << runNumber << endl;
  if (eventNumber > 0) {printf(" event number %i \n", eventNumber); }
  cout << " quality cut parameters: " << endl;
  cout << "  MIN_TRIG_PHISECTORS = " << minTrigPhisectors << endl;
  cout << "  MAX_TRIG_PHISECTORS = " << maxTrigPhisectors << endl;
  cout << "  SURF_SAT_THRESHOLD = " << surfSatThreshold << endl;
  cout << "  DC_OFFSET_THRESHOLD = " << dcOffsetThreshold << endl;
  cout << "  MIN_WF_LENGTH = " << minWfLength << endl;
  cout << "  NADIR_NOISE_THRESHOLD = " << nadirNoiseThreshold << endl;
  
  cout << "  SINE_SUBTRACT_THRESHOLD = " << sineSubThreshold << endl;
  cout << "  SINE_SUBTRACT_FREQ_BANDS = " << sineSubFreqBands << endl;

  if (useTrigNsWindow) {
    printf("  CAL_PULSER = %s \n", calPulser);  
  }
  gRandom->SetSeed(0);
  
  FFTtools::loadWisdom("wisdom.dat"); 
  char* dataDirLocal;
  char* dataDirRemote;
  if (simulationMode) {
    //dataDirLocal = getenv("ANITA_DATA_SIM_DIR");
    //dataDirLocal = "/fs/scratch/PAS0174/anita/simdata";
    dataDirLocal = "/fs/scratch/PAS0174/anita/simdata_thermalOnly";
    dataDirRemote = dataDirLocal;
  } else {
    dataDirLocal = getenv("ANITA_DATA_LOCAL_DIR");
    dataDirRemote = getenv("ANITA_DATA_REMOTE_DIR");  
  }
  char fileDir[1024];
  char calEventFileDir[1024];
  if (simulationMode) {
    sprintf(fileDir, "%s/run%i/", dataDirLocal, runNumber);
    sprintf(calEventFileDir, "%s/run%i/", dataDirRemote, runNumber);
  } else {
    sprintf(fileDir, "%s/run%03i/", dataDirLocal, runNumber);
    sprintf(calEventFileDir, "%s/run%03i/", dataDirRemote, runNumber);
  }

  char filepath[1024];
  if (simulationMode) {
    sprintf(filepath,"%s%s%i.root", fileDir, "SimulatedAnitaHeadFile" ,runNumber);        
  } else {
    if (useTrigNsWindow) {
      sprintf(filepath,"%s%s%03i.root", fileDir, "timedHeadFile" ,runNumber);          
    } else {
      sprintf(filepath,"%s%s%03i.root", fileDir, HEADER_FILENAME_ID ,runNumber);    
    }
  }
  cout << "calculated header file path is " << filepath << endl;
  TChain* headTree = new TChain("headTree");
  headTree->Add(filepath);

  RawAnitaHeader* header = 0;
  headTree->SetBranchAddress("header", &header);
  headTree->BuildIndex("eventNumber");
  cout << "header tree contains " << headTree->GetEntries() << " entries" << endl;
  //int runStartTime = 0;

  if (simulationMode) {
    sprintf(filepath,"%sSimulatedAnitaEventFile%i.root", calEventFileDir, runNumber);    
  } else {
    sprintf(filepath,"%scalEventFile%03i.root", calEventFileDir, runNumber);    
  }
  cout << "calculated event file path is " << filepath << endl;
  TChain* eventTree = new TChain("eventTree");
  eventTree->Add(filepath);
  CalibratedAnitaEvent* calEvent = 0;
  UsefulAnitaEvent* event = 0;
  if (simulationMode) {
    eventTree->SetBranchAddress("event", &event);
  } else {    
    eventTree->SetBranchAddress("event", &calEvent);
  }
  eventTree->BuildIndex("eventNumber");

  if (simulationMode) {
    sprintf(filepath,"%sSimulatedAnitaGpsFile%i.root", fileDir, runNumber);    
  } else {
    sprintf(filepath,"%sgpsEvent%03i.root", fileDir, runNumber);    
  }
  cout << "calculated gps file path is " << filepath << endl;
  TChain* gpsTree = new TChain("adu5PatTree");
  gpsTree->Add(filepath);
  Adu5Pat* adu5Pat = 0;
  int gpsEventNum = 0;
  gpsTree->SetBranchAddress("pat", &adu5Pat);
  gpsTree->SetBranchAddress("eventNumber", &gpsEventNum);
  gpsTree->BuildIndex("eventNumber");

  TChain* truthTree = 0;
  TruthAnitaEvent* simTruth = new TruthAnitaEvent;
  //double inWeight = 0;
  if (simulationMode) {
    truthTree = new TChain("truthAnitaTree");
    sprintf(filepath, "%sSimulatedAnitaTruthFile%i.root", fileDir, runNumber);
    truthTree->Add(filepath);
    truthTree->SetBranchAddress("truth", &simTruth);

    //resultTree->SetBranchAddress("pNu", &pnu);
    //resultTree->SetBranchAddress("weight", &weight);
    //resultTree->SetBranchAddress("truthPhi", &truthPhi);
    //resultTree->SetBranchAddress("truthTheta", &truthTheta);

    truthTree->BuildIndex("eventNumber");
  }

  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  //for (int a=0; a<NUM_SEAVEYS; ++a) {
  //  double ax, ay, az;
  //  geom->getAntXYZ(a, ax, ay, az, AnitaPol::kHorizontal);
  //  printf("antenna %i  (%f,%f,%f) \n", a, ax, ay, az);
  //}

  // set up the analyzer
  AnalysisConfig* cfg = new AnalysisConfig();
  //cfg->combine_npad = 1;
  //cfg->combine_unfiltered = false;
  cfg->normalization_option = normalizationType;
  
  cfg->correlator_nphi = 180; 
  cfg->correlator_ntheta = 85; 
  cfg->correlator_theta_lowest = 60; 
  cfg->correlator_theta_highest = 25; 
  cfg->zoomed_nant = 0; 
  cfg->combine_nantennas = 10;
  cfg->nmaxima = 3;

  printf("normalization type is %s \n", AnalysisConfig::getNormalizationString(normalizationType));
  //cfg->td_pad_factor = timeDomainPadFactor;
  bool iActive = true;
  cout << "Analyzer interactive code is " << iActive << endl;
  Analyzer* analyzer = new Analyzer(cfg, iActive);
  cout << "Analyzer instantiated" << endl;

  // obtain baseline db power spectrum either from file or calculated inline
  vector<vector<TGraphAligned*> > noiseSamplesG(2);
  //vector<vector<AnalysisWaveform*> > noiseSamples(2);
  //FFTWComplex* noiseFreqUpdate = noiseSamples[0][0]->updateFreq();
  bool justUse0 = false;
  if (getBaselineFromFile) {       // from file 
    TFile* baselineFile = new TFile(baselineFilename);
    TTree* blSampleTree = (TTree*) baselineFile->Get("blSampleSm");
    printf("baseline sample tree has %lli entries \n", blSampleTree->GetEntries());
    for (int p=0; p<2; ++p) for (int a=0; a<NUM_SEAVEYS; ++a) {
      noiseSamplesG[p].push_back(0);
    }
    TGraph* thisGraph=0;
    int ant, pol;
    blSampleTree->SetBranchAddress("ant", &ant);
    blSampleTree->SetBranchAddress("pol", &pol);
    blSampleTree->SetBranchAddress("blWaveform", &thisGraph);   
    for (int e=0; e<blSampleTree->GetEntries(); ++e) {
      blSampleTree->GetEntry(e);
      //printf("got baseline sample entry %i pol %i ant %i\n", e, pol, ant);
      noiseSamplesG[pol][ant] = new TGraphAligned(*thisGraph);
      //noiseSamplesG[pol][ant]->Set(min(160, noiseSamplesG[pol][ant]->GetN()));
    }
    baselineFile->Close();
  } else {                          // calculate inline: broken line (dB) power spectrum
    for (int p=0; p<2; ++p) for (int a=0; a<NUM_SEAVEYS; ++a) {
      //noiseSamples[p].push_back(new AnalysisWaveform(260));
      noiseSamplesG[p].push_back(new TGraphAligned(131));
    }
    double* noiseFreqX = noiseSamplesG[0][0]->GetX();
    double* noiseFreqY = noiseSamplesG[0][0]->GetY();
    // for now just one waveform, of rayleigh noise
    double noiseDf = 65.0/6.0/1000.0;
    double firstSlope = 50;
    double firstInt = -20;
    double secondSlope = -10;
    double secondInt = -9;
    double firstBreak = 0.15;
    double secondBreak = 0.4;
    for (int k=0; k<131; ++k) {
      double freq = k*noiseDf;
      //double env = (freq < firstBreak) ? pow(10, (firstSlope*freq + firstInt)/10.0) : pow(10, (secondSlope*freq + secondInt)/10.0);
      double powdB = (freq < firstBreak) ? (firstSlope*freq + firstInt) : 
              (freq > secondBreak) ? (secondSlope*freq + secondInt) : firstSlope*firstBreak + firstInt;
      //double powRaw = pow(10, powdB/10);
      //noiseFreqUpdate[k].re = env*1000.0;
      //noiseFreqUpdate[k].im = 0;
      noiseFreqX[k] = freq;
      noiseFreqY[k] = powdB;
    }
  }
  if (justUse0) {
    for (int p=0; p<2; ++p) {
      for (int a=1; a<NUM_SEAVEYS; ++a) {
        noiseSamplesG[p][a] = noiseSamplesG[p][0];
      }
    }
  }
  
  for (int p=0; p<0; ++p) for (int a=0; a<0; ++a) {
    for (int k=0; k<noiseSamplesG[p][a]->GetN(); ++k) {
      //printf("noise sample p=%i, a=%i, k=%i, freq=%f, powerdB=%f \n", p, a, k, noiseSamplesG[p][a]->GetX()[k], noiseSamplesG[p][a]->GetY()[k]);
      printf("   new way freq=%f, powerdB=%f \n", noiseSamplesG[p][a]->GetX()[k], noiseSamplesG[p][a]->GetY()[k]);
    }
  }


  double fmins[1] = {0.2};
  double fmaxs[1] = {1.3};
  printf("  instantiating sine subtract filter \nb");
  UCorrelator::SineSubtractFilter* sineSubFilter = new UCorrelator::SineSubtractFilter(sineSubThreshold, 0, sineSubFreqBands, fmins, fmaxs);
  //AdaptiveNotchFilter* adapFilter = new AdaptiveNotchFilter(adaptiveNotchFilterOrder, adaptiveNotchFilterWidth, adaptiveNotchFilterCausal);
  
  GeometricFilter* geomFilter = new GeometricFilter(noiseSamplesG);
  printf("geom filter instantiated \n");
  geomFilter->setDbCut(2.0);
  HybridFilter* hybridFilter = new HybridFilter();
  vector<AnitaEventSummary> eventSummaryList(0);
  vector<RawAnitaHeader> headerList(0);
  vector<Adu5Pat> gpsList(0);
  vector<int> circPolList(0);
  vector<float> pNuList(0);
  vector<float> pWgtList(0);
  vector<float> truthPhiList(0);
  vector<float> truthThetaList(0);
  //vector<vector<int> > aCuts0(0);
  //vector<vector<int> > aCuts1(0);
  vector<vector<int> > qCuts0(0);
  vector<vector<int> > qCuts1(0);
  vector<vector<vector<float> > > maxPwValList(0);
  vector<vector<int > > qualityCuts(2, vector<int>(NUM_Q_CUTS));   // 1 = passed cut
  vector<TCanvas*> canvi(4);
  vector<TGraph*> cohPower(4);
  vector<TGraph*> cohWfGr(4);
  vector<TGraph*> hilbEnvGr(4);
  vector<TH2D*> corrMap(4);
  vector<TH2D*> zoomCorrMap(4);
  int qualityCutCounts[4][NUM_Q_CUTS] = {{0}};

  gStyle->SetTitleSize(0.07, "xyz");
  gStyle->SetTitleFontSize(0.07);
  gStyle->SetTitleOffset(0.7, "x");
  gStyle->SetTitleOffset(0.6, "y");
  gStyle->SetLabelSize(0.06, "xyz");

  TFile* filterOutputFile = new TFile("filterOutput.root", "RECREATE");
  // set up filter strategies, with and without hybrid
  FilterStrategy* strategy[2] = {0};
  for (int cp=0; cp<2; ++cp) {
    strategy[cp] = new FilterStrategy(filterOutputFile);
    //strategy[cp]->addOperation(new SimplePassBandFilter(0.2, 1.2)); 
    //if (useFilter) {strategy[cp]->addOperation(sineSubFilter, true);}
    if (filterOption == FILTER_OPTION_ABBY) {
      printf("Warning: Abby's filtering is unsupported in this version will be applied \n");
      //UCorrelator::applyAbbysFilterStrategy(strategy[cp]);
    } else if (filterOption == FILTER_OPTION_SINE_SUBTRACT) {
      printf("Sine-subtraction filtering will be applied \n");
      strategy[cp]->addOperation(sineSubFilter, true);
    //} else if (filterOption == FILTER_OPTION_ADAPTIVE_NOTCH) {
    //  printf("Adaptive notch filtering will be applied \n");
    //  strategy[cp]->addOperation(adapFilter, true);
    } else if (filterOption == FILTER_OPTION_GEOMETRIC) {
      printf("Geometric notch filtering will be applied \n");
      strategy[cp]->addOperation(geomFilter, true);
    } else {
      printf("no filtering will be applied \n");
    }
    strategy[cp]->addOperation(new ALFAFilter);
  }
  
  // for circ pol
  strategy[1]->addOperation(hybridFilter, false);   
  
  printf("Filtering: \n");
  char filterDescLn[1024] = "";
  printf(" LnPol \n");
  for (int op=0; op<strategy[0]->nOperations(); ++op) {
    printf("   %s \n", strategy[0]->getOperation(op)->description());
    strcat(filterDescLn, strategy[0]->getOperation(op)->description());
    strcat(filterDescLn, "\n");
  }

  char filterDescCr[1024] = "";
  printf(" CrPol \n");
  for (int op=0; op<strategy[1]->nOperations(); ++op) {
    printf("   %s \n", strategy[1]->getOperation(op)->description());
    strcat(filterDescCr, strategy[1]->getOperation(op)->description());
    strcat(filterDescCr, "\n");
  }

  //for (int k=0; k<2; ++k) {qualityCuts[k] = vector<int>(NUM_Q_CUTS);}
  int numEvents = 0;
  float avgSnr = 0;
  int snrEventCount = 0;
  vector<vector<float> > maxPwVal(2, vector<float>(2));
  AnitaEventSummary* eventSummary = new AnitaEventSummary;
  for (int e=startEntryNum; e<headTree->GetEntries() && numEvents<maxEvents; e+=dbStride) {
    vbprintf("-----------------------------------------------------------------------------------------------\n");
    // get the event and the gps
    if (maxEvents == 1 && eventNumber > 0) {
      headTree->GetEntryWithIndex(eventNumber);
    } else {
      vcout << endl << "event sequence number " << numEvents << endl;
      vcout << " event entry " << e << " of " << headTree->GetEntries() << endl;
      headTree->GetEntry(e);
    }
    
    vbprintf("event number %8i, trigger time %10i:%09i, type=%4s, event sequence %8i, header entry %8i of %8lli \n", 
        header->eventNumber, header->triggerTime, header->triggerTimeNs, header->trigTypeAsString(), numEvents, e, headTree->GetEntries());
    if (strcmp(header->trigTypeAsString(), "RF") != 0) {
      vbprintf ("  rejecting as non-RF \n");
      delete header; header=0;
      continue;
    }
    gpsTree->GetEntryWithIndex(header->eventNumber);
    if (gpsEventNum != header->eventNumber) {
      printf("  warning: gps database error - header event %i gps event %i  skipping this event \n", header->eventNumber, gpsEventNum);
      continue;
    }
    UsefulAdu5Pat* gps = new UsefulAdu5Pat(adu5Pat);
    simTruth = 0;
    if (truthTree) {
      truthTree->GetEntryWithIndex(header->eventNumber);
    }
    
    //int nsDiff = 0;
    if (useTrigNsWindow) {
      bool calPulse = false;
      if (strcmp(calPulser, "WAIS")==0) {
        calPulse = UCorrelator::isWAISHPol(gps, header, cfg);
      }
      else if (strcmp(calPulser, "LDB")==0) {
        calPulse = UCorrelator::isLDB(header, cfg);
      }
      // sidestep seg fault on this event :148-9026525
      // TODO fix the seg fault
      //calPulse &= (header->eventNumber != 9026525);
      //int sourceDelay = (int)(gps->getTriggerTimeNsFromSource(ksLat, ksLon, ksAlt));
      //for (int k=0; k<pulserNsTimes.size(); ++k) {
      //  nsDiff = sourceDelay - header->triggerTimeNs + pulserNsTimes[k];
      //  vbprintf("trigNs=%i   sourceDelay=%i   pulseTime=%i   nsDiff=%i \n", header->triggerTimeNs, sourceDelay, pulserNsTimes[k], nsDiff);
      //  calPulse |= (nsDiff > nsDiffMin && nsDiff < nsDiffMax);
      //}
      if (!calPulse) {
        vcout << " rejecting event: trigNs not in cal pulse window" << endl;
        delete gps; gps=0;
        delete header; header=0;
        continue;
      }
    }
    // add event number var, because maybe it will save time, idk... coding is hard.
    int eventNumberTemp = header->eventNumber;

    printf("processing event number %8i, trigger time %10i:%09i, type=%4s, event sequence %8i, header entry %8i of %8lli \n", 
            eventNumberTemp, header->triggerTime, header->triggerTimeNs, header->trigTypeAsString(), numEvents, e, headTree->GetEntries());


    //add skip for corrupted? or imcomplete event entries (but only for sim) -Jacob
    if ( !simulationMode && 
         (
          eventNumberTemp == 17907826 ||
          eventNumberTemp == 55868364 ||
          eventNumberTemp == 58665137 ||
          eventNumberTemp == 58666148 ||
          eventNumberTemp == 65258379 ||
          eventNumberTemp >  84667494 ||
          eventNumberTemp == 72763437 ||
          eventNumberTemp == 64298532 ||
          eventNumberTemp == 53049784 ||
          eventNumberTemp == 37771298 ||
          eventNumberTemp == 31455791 ||
          eventNumberTemp == 58175220 ||
          eventNumberTemp == 58508540   
         )
       )
    {
      printf("Skipping event number %8i to avoid segfault", eventNumberTemp);
      continue;
    }

    eventTree->GetEntryWithIndex(eventNumberTemp);
    if (!simulationMode) {
      event = new UsefulAnitaEvent(calEvent);
    }    
    //cout << "event number " << eventNumber << endl;
    int cNum = 0;
    delete analyzer; analyzer=0;
    analyzer = new Analyzer(cfg, iActive);    // careful now airstream driver  workaround for weird Analyzer state failure
    AnitaEventSummary::PointingHypothesis linPolPeaks[2];
    for (int circPol= (doLinPol ? 0 : 1); circPol<(doCircPol ? 2 : 1); ++circPol) {
      for (int p=0; p<2; ++p) for (int k=0; k<32; ++k)  {qualityCuts[p][k]=0; /*analysisCuts[p][k]=0;*/}
      //FilterStrategy* strategy = new FilterStrategy();
      //strategy->addOperation(sineSubFilter);   
      //if (circPol) strategy->addOperation(hybridFilter);   
      //FilteredAnitaEvent* filteredEvent = new FilteredAnitaEvent(event, strategy, gps, header);
      FilteredAnitaEvent* filteredEvent = new FilteredAnitaEvent(event, strategy[circPol], gps, header);
      //AnitaEventSummary* eventSummary = new AnitaEventSummary;
      //eventSummary->zeroInternals();
      analyzer->analyze(filteredEvent, eventSummary);
      if (useTrigNsWindow) {
        if (strcmp(calPulser, "WAIS")==0) {
          vbprintf("WAIS is at %f,%f (pl) \n", eventSummary->wais.phi, -eventSummary->wais.theta);
        } else if (strcmp(calPulser, "LDB")==0) {
          vbprintf("LDB is at %f,%f (pl) \n", eventSummary->ldb.phi, -eventSummary->ldb.theta);
        }
      }
      for (uint8_t p=AnitaPol::kHorizontal; p<AnitaPol::kNotAPol; ++p) {
        char polChar = AnitaPol::polAsChar((AnitaPol::AnitaPol_t)p);
        AnitaPol::AnitaPol_t pol = (AnitaPol::AnitaPol_t) p;
        UShort_t l3TrigPat = 0;
        if (circPol) switch(polChar) {case 'H' : polChar='L'; break; case 'V' : polChar='R'; break;}
        switch (polChar) {
          case '-' : 
          case 'H' : l3TrigPat = header->l3TrigPatternH; break;
          case 'V' : l3TrigPat = header->l3TrigPattern; break;
          default : l3TrigPat = header->l3TrigPatternH | header->l3TrigPattern;
        }

        vcout << " polarization " << polChar << ":" << endl;
                vector<const TGraph*> rawWf(0);
        for (short an=0; an<NUM_SEAVEYS_ANITA3; ++an) {
          const AnalysisWaveform* thisRawWf = filteredEvent->getFilteredGraph(an, pol);
          const TGraph* thisWf = thisRawWf->even();
          rawWf.push_back(thisWf);
        }
        //for (int k=0; k<eventSummary->nPeaks[p]; ++k) {
        AnitaEventSummary::WaveformInfo thisWfI = eventSummary->coherent_filtered[p][0];
        if  (circPol == 0) {
          linPolPeaks[p] = eventSummary->peak[p][0];
        }
        //AnitaEventSummary::PointingHypothesis thisPeak = eventSummary->peak[p][k];
        //vcout << "  peak " << k << " (" << thisPeak.phi << "," << -thisPeak.theta << ") value=" << thisPeak.value << endl;
        //vcout << "   coherent reconstruction: " << endl;
        //double sI, sQ, sU; // sV;
        //sI = thisWfI.I;
        //if (circPol) {
        //if (false) {
          //sQ = thisWfI.U;
          //sU = -thisWfI.V;
          //sV = -thisWfI.Q;
        //} else {
          //sQ = thisWfI.Q;
          //sU = thisWfI.U;
          //sV = thisWfI.V;
        //}
        //vcout << "    Stokes parameters: I=" << sI << "  Q=" << sQ << "  U=" << sU << "  V=" << sV << endl; 
        //float linPol = sqrt(sQ*sQ + sU*sU) / sI;
        //vcout << "    linPol proportion " << linPol << endl;
        vbprintf(" coh snr0=%f \n", thisWfI.snr);
        vbprintf(" coh snr1=%f \n", thisWfI.snr1);
        vbprintf(" coh snr2=%f \n", thisWfI.snr2);
        vbprintf(" coh snr3=%f \n", thisWfI.snr3);
        if (circPol && p==0) {
          avgSnr += thisWfI.snr;
          ++snrEventCount;
        }
        //}


        // calculate cut results
        qualityCuts[p][Q_CUT_TRIGGER_TYPE] = (strcmp(header->trigTypeAsString(), "RF")==0) ? 0 : 1; 
        vcout << "  trigger type cut: " << qualityCuts[p][Q_CUT_TRIGGER_TYPE] << endl;
        qualityCuts[p][Q_CUT_NO_TRIGGER] = (l3TrigPat == 0);
        vcout << "  no-trigger cut: " << qualityCuts[p][Q_CUT_NO_TRIGGER] << endl;
        int numPhiSectors = 0;
        //UShort_t trigPatWrk = l3TrigPat;
        vbprintf("triggering phi-sectors: ");
        for (int ps=0; ps<16; ps++) {
          if (l3TrigPat & (1<<ps)) {
            vbprintf(" %2i", ps);
            ++numPhiSectors;
          }
        }
        vbprintf("\n");
        vcout << "    " << numPhiSectors << " L3-triggering phi-sectors" << endl;
        qualityCuts[p][Q_CUT_NUM_PHISECTOR] = (numPhiSectors < minTrigPhisectors) ? 1 : 0;
        vcout << "  triggering phi-sectors cut: " << qualityCuts[p][Q_CUT_NUM_PHISECTOR] << endl;
        qualityCuts[p][Q_CUT_PAYLOAD_BLAST] = (numPhiSectors > maxTrigPhisectors) ? 1 : 0;
        vcout << "  payload blast cut: " << qualityCuts[p][Q_CUT_PAYLOAD_BLAST] << endl;

        // short waveforms cut
        vector<short> shortWfRejects;
        for (short an=0; an<NUM_SEAVEYS_ANITA3; ++an) {
          //TGraph* thisWf = rawWf[an];
          int ch = AnitaGeomTool::getChanIndexFromAntPol(an, pol);
          if (event->fNumPoints[ch] < minWfLength) {
            shortWfRejects.push_back(an);
          }
        }
        qualityCuts[p][Q_CUT_SHORT_WF] = (shortWfRejects.size() > 0) ? 1 : 0;
        vcout << "  short waveforms cut: " << qualityCuts[p][Q_CUT_SHORT_WF] << endl;

        // wf saturation cut
        set<short> antSatRejects;        
        for (short an = 0; an<NUM_SEAVEYS_ANITA3; ++an) {
          const TGraph* thisWf = rawWf[an];
          float maxPos = *max_element(thisWf->GetY(), thisWf->GetY()+thisWf->GetN());
          float maxNeg = abs(*min_element(thisWf->GetY(), thisWf->GetY()+thisWf->GetN()));
          if (max(maxPos,maxNeg) > surfSatThreshold) {
            antSatRejects.insert(an);
          }
        }
        qualityCuts[p][Q_CUT_WF_SATURATION] = (antSatRejects.size() > 3) ? 1 : 0;
        vcout << "  saturation cut: " << qualityCuts[p][Q_CUT_WF_SATURATION] << endl;

        // DC offset cut
        set<short> dcOffsetRejects;
        for (short an=0; an<NUM_SEAVEYS_ANITA3; ++an) {
          float meanV = abs(rawWf[an]->GetMean(2));
          //vcout << "      antenna " << an << ", mean value=" << meanV << endl;
          if (meanV > dcOffsetThreshold) {
            dcOffsetRejects.insert(an);
          }
        }
        qualityCuts[p][Q_CUT_WF_DC_OFFSET] = (dcOffsetRejects.size() > 0) ? 1 : 0;
        vcout << "  DC offset cut: " << qualityCuts[p][Q_CUT_WF_DC_OFFSET] << endl;

        // nadir noise cut (top/bottom ring peak ratio)
        float peakTopRing = 0;
        float peakBottomRing = 0;
        for (short an=0; an<16; ++an) {
          const TGraph* thisWf = rawWf[an];      
          float maxPos = *max_element(thisWf->GetY(), thisWf->GetY()+thisWf->GetN());
          float maxNeg = abs(*min_element(thisWf->GetY(), thisWf->GetY()+thisWf->GetN()));
          peakTopRing = max(peakTopRing, maxPos);      
          peakTopRing = max(peakTopRing, maxNeg);      
        }
        for (short an=32; an<48; ++an) {
          const TGraph* thisWf = rawWf[an];      
          float maxPos = *max_element(thisWf->GetY(), thisWf->GetY()+thisWf->GetN());
          float maxNeg = abs(*min_element(thisWf->GetY(), thisWf->GetY()+thisWf->GetN()));
          peakBottomRing = max(peakBottomRing, maxPos);      
          peakBottomRing = max(peakBottomRing, maxNeg);      
        }
        float peakRatio = peakBottomRing/peakTopRing;
        vcout << "top / bottom ring maximum peak is " << peakTopRing << " / " << peakBottomRing << ",  ratio=" << peakRatio << endl;
        qualityCuts[p][Q_CUT_NADIR_NOISE] = (peakTopRing/peakBottomRing < nadirNoiseThreshold) ? 1 : 0;
        vcout << "  nadir noise cut " << qualityCuts[p][Q_CUT_NADIR_NOISE] << endl;

        if (eventSummary->nPeaks[p] < 2) {
          printf("   Warning: event %i nPeaks[%c] = %i \n", eventSummary->eventNumber, polChar, eventSummary->nPeaks[p]);
          fflush(stdout);
        }
        for (int k=0; k<eventSummary->nPeaks[p]; ++k) {
          vbprintf("  peak %1i (%.3f,%.3f): %.3f \n", k, eventSummary->peak[p][k].phi, eventSummary->peak[p][k].theta, eventSummary->peak[p][k].value);
        }
        //AnitaEventSummary::PointingHypothesis mainPeak = eventSummary->peak[p][0];
        for (int k=0; k<2; ++k) {maxPwVal[p][k]=-1000.;}        
        // if we're on the circPol cycle, find the maximum value this circPol map on the regions around the H- and VPol peaks
        if (circPol == 1) {
          // there should be two entries in the event summary list: one H and one V; lp is this lonPol index
          for (int lp = 0; lp<2; ++lp) {
            AnitaEventSummary::PointingHypothesis thisPeak = linPolPeaks[lp];
            float loTheta = max(thisPeak.theta-5.0, -cfg->correlator_theta_lowest);
            float hiTheta = min(thisPeak.theta+5.0, cfg->correlator_theta_highest);
            float loPhi = thisPeak.phi - 5.0;
            float hiPhi = thisPeak.phi + 5.0;
            float dPhi = 360.0 / cfg->correlator_nphi;
            // careful now airstream driver; signs of the two theta limits are opposite; that's why we add not subtract
            float dTheta = (cfg->correlator_theta_lowest + cfg->correlator_theta_highest) / cfg->correlator_ntheta;   
            //vbprintf("   peak window %cPol is (%f,%f)-(%f,%f)\n", lp==0?'H':'V', loPhi, -loTheta, hiPhi, -hiTheta);
            // find the maximum values in the H and V peak windows, for this circPol
            //vbprintf(" pol code %i \n", pol);
            TH2D* thisMap = (TH2D*)(analyzer->getCorrelationMap(pol));
            for (float phi=loPhi; phi<hiPhi; phi+=dPhi) {
              float thisPhi = phi + dPhi/2.0;
              thisPhi = fmod(thisPhi + 360, 360);
              for (float theta=-loTheta; theta>-hiTheta; theta-=dTheta) {
                float thisTheta = theta + dTheta/2.0;
                int thisBin = thisMap->FindBin(thisPhi, thisTheta);
                float thisVal = thisMap->GetBinContent(thisBin);
                //vbprintf("  phi=%f, theta=%f, val=%f \n", thisPhi, thisTheta, thisVal);
                maxPwVal[p][lp] = max(thisVal, maxPwVal[p][lp]);
              }
            }
            vbprintf("     max value in peak window is %f \n", maxPwVal[p][lp]);
          }
        }
        //if (doGraphics /*&& maxEvents==1*/) {
        if (maxEvents==1) {
          cout << "drawing polarization " << polChar << endl; 
          char canvName[8]; sprintf(canvName, "%cpolCanvas", polChar);
          char canvTitle[8]; sprintf(canvTitle, "%c-Pol", polChar);
          canvi[cNum] = new TCanvas(canvName, canvTitle, 1320, 600);
          int padNum = 0;
          //canvi[cNum]->Divide(2,2);
          canvi[cNum]->Divide(2,2);

          TMarker* ksMarker = 0;
          if (strcmp(calPulser, "WAIS")==0) {
            ksMarker = new TMarker(eventSummary->wais.phi, -eventSummary->wais.theta, kPlus);
          } else if (strcmp(calPulser, "LDB")==0) {
            ksMarker = new TMarker(eventSummary->ldb.phi, -eventSummary->ldb.theta, kPlus);
          } 
          
          char corrMapTl[128]; sprintf(corrMapTl, "Correlation map: event %i, run %i, %cPol filter=%i;#phi      ;#theta     ", header->eventNumber, runNumber, polChar, filterOption);
          corrMap[cNum] = new TH2D((const TH2D&)(*(analyzer->getCorrelationMap(pol))));
          corrMap[cNum]->SetTitle(corrMapTl);
          corrMap[cNum]->SetTitleSize(0.07, "xy");
          corrMap[cNum]->SetLabelSize(0.06, "xyz");
          corrMap[cNum]->SetTitleOffset(0.7, "x");
          corrMap[cNum]->SetTitleOffset(0.5, "y");
          corrMap[cNum]->SetLabelOffset(0.01, "xy");
          corrMap[cNum]->SetAxisColor(kWhite, "xy");
          canvi[cNum]->cd(++padNum);     
          gPad->SetRightMargin(0.15);
          corrMap[cNum]->Draw("");
          //corrMap[cNum]->GetZaxis()->SetRangeUser(-0.05, 0.07);      // forced scaling for nonimpulsive events
          //corrMap[cNum]->GetZaxis()->SetRangeUser(-0.05, 0.17);      // forced scaling for cal pulse events
          corrMap[cNum]->Draw("COLZ");
          gStyle->SetOptStat("");
          if (knownSource) {ksMarker->Draw();}
          canvi[cNum]->Draw();

          zoomCorrMap[cNum] = new TH2D((const TH2D&)(*(analyzer->getZoomedCorrelationMap(pol, 0))));
          zoomCorrMap[cNum]->SetTitle(corrMapTl);
          zoomCorrMap[cNum]->SetTitleSize(0.07, "xy");
          zoomCorrMap[cNum]->SetLabelSize(0.06, "xyz");
          zoomCorrMap[cNum]->SetTitleOffset(0.7, "x");
          zoomCorrMap[cNum]->SetTitleOffset(0.5, "y");
          zoomCorrMap[cNum]->SetLabelOffset(0.01, "xy");
          zoomCorrMap[cNum]->SetAxisColor(kWhite, "xy");
          // temporarily disabled: zoomed correlation map
          canvi[cNum]->cd(++padNum);     
          gPad->SetRightMargin(0.15);
          //gPad->SetLeftMargin(0.3);
          zoomCorrMap[cNum]->Draw("");
          //zoomCorrMap[cNum]->GetZaxis()->SetRangeUser(-0.03, 0.05);    // force scaling
          zoomCorrMap[cNum]->Draw("COLZ");
          if (knownSource) {ksMarker->Draw();}
          //canvi[cNum]->Draw();
          
          /*
          char cohPowerTl[128]; sprintf(cohPowerTl, "Reconstructed power spectrum;freq (GHz)       ;power (dB, normalized)     ");
          const TGraphAligned* cohPowerGA = analyzer->getCoherentPower(pol, 0);
          cohPower[cNum] = new TGraph((TGraph)*cohPowerGA);
          cohPower[cNum]->SetTitle(cohPowerTl);
          canvi[cNum]->cd(++padNum);     
          //cohPower[cNum]->Draw("AL");
          cohPower[cNum]->GetXaxis()->SetTitleSize(0.05);
          cohPower[cNum]->GetYaxis()->SetTitleSize(0.05);
          cohPower[cNum]->GetXaxis()->SetTitleOffset(0.95);
          cohPower[cNum]->GetYaxis()->SetTitleOffset(0.9);
          cohPower[cNum]->Draw("AL");
          */
          char cohWfTl[128]; sprintf(cohWfTl, "Coherent Reconstruction;time (ns)   ;V(mV)   ");
          const TGraphAligned* cohWfGA = analyzer->getCoherent(pol, 1)->even();
          cohWfGr[cNum] = new TGraph((TGraph)*cohWfGA);
          cohWfGr[cNum]->SetTitle(cohWfTl);
          canvi[cNum]->cd(++padNum);     
          cohWfGr[cNum]->Draw("AL");
          cohWfGr[cNum]->GetXaxis()->SetTitleSize(0.05);
          cohWfGr[cNum]->GetYaxis()->SetTitleSize(0.05);
          cohWfGr[cNum]->GetXaxis()->SetTitleOffset(0.9);
          cohWfGr[cNum]->GetYaxis()->SetTitleOffset(0.95);
          cohWfGr[cNum]->GetYaxis()->SetRangeUser(-100, 100);  // force scaling cal pulses
          cohWfGr[cNum]->Draw("AL");
          
          char hilbEnvTl[128]; sprintf(hilbEnvTl, "Reconstructed Hilbert Envelope;time (ns)   ;V(mV)   ");
          const TGraphAligned* hilbEnvGA = analyzer->getCoherent(pol, 1)->hilbertEnvelope();
          hilbEnvGr[cNum] = new TGraph((TGraph)*hilbEnvGA);
          hilbEnvGr[cNum]->SetTitle(hilbEnvTl);
          canvi[cNum]->cd(++padNum);     
          hilbEnvGr[cNum]->Draw("AL");
          hilbEnvGr[cNum]->GetXaxis()->SetTitleSize(0.05);
          hilbEnvGr[cNum]->GetYaxis()->SetTitleSize(0.05);
          hilbEnvGr[cNum]->GetXaxis()->SetTitleOffset(0.9);
          hilbEnvGr[cNum]->GetYaxis()->SetTitleOffset(0.95);
          //hilbEnv[cNum]->GetYaxis()->SetRangeUser(-100, 100);  // force scaling cal pulses
          hilbEnvGr[cNum]->Draw("AL");
          
          
          
          // TODO plot the hilbert envelope
          
          char canvasFilename[1024];
          sprintf(canvasFilename, "results/plots/interfMap_%i_%i_%cPol.%s", runNumber, eventNumber, polChar, "png");
          canvi[cNum]->SaveAs(canvasFilename);
          sprintf(canvasFilename, "results/plots/interfMap_%i_%i_%cPol.%s", runNumber, eventNumber, polChar, "root");
          canvi[cNum]->SaveAs(canvasFilename);
          sprintf(canvasFilename, "results/plots/interfMap_%i_%i_%cPol.%s", runNumber, eventNumber, polChar, "eps");
          canvi[cNum]->SaveAs(canvasFilename);
          ++cNum;
        }
        
        //if (p==0 && doGraphics) {
        if (doGraphics) {
          int cn = 0;
          char pChar;
          UShort_t trigPat = 0;
          switch (p) {
           case 0 : pChar = circPol ? 'L' : 'H'; trigPat = circPol ? header->l3TrigPatternH | header->l3TrigPattern : header->l3TrigPatternH; break;
           case 1 : pChar = circPol ? 'R' : 'V'; trigPat = circPol ? header->l3TrigPatternH | header->l3TrigPattern : header->l3TrigPattern; break;
          }
          if (false) {
            sprintf(name, "wfCanvF_%cpol", pChar);          
            TCanvas* wfCanvF = new TCanvas(name, name, 1200, 800);
            wfCanvF->Divide(8,6,0,0);
            for (int ant=0; ant<NUM_SEAVEYS; ++ant) {
              wfCanvF->cd(++cn);
              int phiSec = ant%16;
              bool trigAnt = ((trigPat & (1<<(phiSec))) ? 1 :0);     
              //printf("ant %i trigPat = %i phiSec %i trig %i \n", ant, trigPat, phiSec, trigAnt);       
              if (trigAnt) {
                gPad->SetFillColor(5);
                gPad->SetFrameFillColor(50);
              }
              //gPad->SetLogy();
              const AnalysisWaveform* thisAWF = filteredEvent->getFilteredGraph(ant, (AnitaPol::AnitaPol_t)p);
              //const TGraphAligned* thisGrA = thisAWF->power();
              const TGraphAligned* thisGrA = thisAWF->even();
              TGraph* thisGr = new TGraph(*((TGraph*)thisGrA));
              thisGr->GetYaxis()->SetRangeUser(-90, 90);
              thisGr->Draw();
              char antNumStr[3];
              sprintf(antNumStr, "%i", ant);
              TText* antNumText = new TText(80, 70, antNumStr);
              antNumText->SetTextSize(0.1);
              antNumText->Draw();
              thisGr = 0;
            }
          }
          if (false) {
            TCanvas* wfCanvU = new TCanvas("wfCanvU", "wfCanvU", 1200, 800);
            wfCanvU->Divide(8,6);
            cn = 0;
            for (int ant=0; ant<NUM_SEAVEYS; ++ant) {
              wfCanvU->cd(++cn);
              gPad->SetLogy();
              const AnalysisWaveform* thisAWF = filteredEvent->getRawGraph(ant, (AnitaPol::AnitaPol_t)p);
              const TGraphAligned* thisGrA = thisAWF->power();
              TGraph* thisGr = new TGraph(*((TGraph*)thisGrA));
              thisGr->Draw();
              thisGr = 0;
              //printf("  unfiltered waveform ant %i, dt=%f, df=%f  \n", ant, thisAWF->deltaT(), thisAWF->deltaF());
            }
          }
          TCanvas* pwrSpecCanv = new TCanvas("pwrSpecCanv", "", 700, 400);
          TGraph* avgSpectrumF = new TGraph();
          filteredEvent->getAverageSpectrum(avgSpectrumF, (AnitaPol::AnitaPol_t)p);
          pwrSpecCanv->cd(0);
          char grTitle[128]; sprintf(grTitle, "Filtered(%i) spectrum avg hPol event %i;GHz;power dB ", filterOption, filteredEvent->getUsefulAnitaEvent()->eventNumber);
          avgSpectrumF->SetTitle(grTitle);
          avgSpectrumF->Draw("AL");
          pwrSpecCanv->SaveAs("results/plots/pwrSpecCanv.png");
          pwrSpecCanv->SaveAs("results/plots/pwrSpecCanv.root");
          vector<TCanvas*> baselineCanvs(0); 
          TGraph* thisBaselineGr = 0;
          TGraph* thisPwrSpecGr = 0;
          TGraph* thisPhaseSpecGr = 0;
          if (false) {
            for (int a=0; a<NUM_SEAVEYS; ++a) {
              if (a==13) {
                thisBaselineGr = 0;
                thisPwrSpecGr = 0;
                thisPhaseSpecGr = 0;              
                char canvTitle[128];
                sprintf(canvTitle,"Baseline Hpol antenna %i", a);
                char canvName[128];
                sprintf(canvName,"blCanv_%i", a);
                TCanvas* thisCanv = 0;
                thisCanv = new TCanvas(canvName, canvTitle, 400, 750);
                thisCanv->Divide(1,3);
                thisBaselineGr = new TGraph(*noiseSamplesG[p][a]);
                thisPwrSpecGr = new TGraph(*filteredEvent->getRawGraph(a, (AnitaPol::AnitaPol_t)p)->powerdB());
                sprintf(grTitle, "Power Spectrum ant%02i %cPol", a, polChar);
                thisPwrSpecGr->SetTitle(grTitle);
                thisPhaseSpecGr = new TGraph(*filteredEvent->getRawGraph(a, (AnitaPol::AnitaPol_t)p)->phase());
                sprintf(grTitle, "Phase Spectrum ant%02i %cPol", a, polChar);
                thisPhaseSpecGr->SetTitle(grTitle);
                thisCanv->cd(1);
                thisBaselineGr->Draw("AL");
                thisCanv->cd(2);
                thisPwrSpecGr->Draw("AL");
                thisCanv->cd(3);
                thisPhaseSpecGr->SetMarkerStyle(2);
                thisPhaseSpecGr->Draw("AL");
              }
            }
          }
        }
      }
      eventSummaryList.push_back(*eventSummary);
      headerList.push_back(*header);
      gpsList.push_back(*gps);
      circPolList.push_back(circPol);
      if (simulationMode) {
        char eventnumStr[1024] = {0};
        sprintf(eventnumStr, "eventNumber == %i",header->eventNumber);
        truthTree->Draw("nuMom:weight:payloadPhi:payloadTheta",eventnumStr,"");

        pNuList.push_back((float)truthTree->GetV1()[0]);
        pWgtList.push_back((float)truthTree->GetV2()[0]);
        truthPhiList.push_back((float)truthTree->GetV3()[0]);
        truthThetaList.push_back((float)truthTree->GetV4()[0]);

        //pNuList.push_back((float)simTruth->nuMom);
        //pWgtList.push_back((float)simTruth->weight);
        //cout << simTruth->weight << endl;
        //truthPhiList.push_back((float)simTruth->payloadPhi);
        //truthThetaList.push_back((float)simTruth->payloadTheta); 
      }
      vcout << "event " << eventSummary->eventNumber << endl;
      vcout << " quality cuts pol 0: ";
      for (int k=0; k<32; ++k) {vcout << qualityCuts[0][k] << " ";}
      vcout << endl;
      vcout << " quality cuts pol 1: ";
      for (int k=0; k<32; ++k) {vcout << qualityCuts[1][k] << " ";}
      vcout << endl;

      qCuts0.push_back(qualityCuts[0]);
      qCuts1.push_back(qualityCuts[1]);
      maxPwValList.push_back(maxPwVal);
      for (int j=0; j<32; ++j) {
        int co = (circPol) ? 2 : 0;
        if (qualityCuts[0][j]>0) {++qualityCutCounts[co][j]; vcout << "    incrementing qCuts[" << co+1 << "][" << j << "]" << endl;}
        if (qualityCuts[1][j]>0) {++qualityCutCounts[co+1][j]; vcout << "    incrementing qCuts[" << co+1 << "][" << j << "]" << endl;}
      }
    
      delete filteredEvent; filteredEvent=0;
      //delete eventSummary; eventSummary=0;
      //delete strategy;    
    }
    delete header; header=0;
    delete adu5Pat; adu5Pat=0;
    delete gps; gps=0;
    delete calEvent; calEvent=0;
    delete event; event=0;    
    ++numEvents;
    if (numEvents%100 == 0) fflush(stdout);
    if (doGraphics) app->Run();
  }

  filterOutputFile->Write();
  filterOutputFile->Close();
  
  //headerFile->Close();
  //calEventFile->Close();
  //gpsFile->Close();
  
  delete analyzer;
  delete cfg;
  delete geom;
  
  avgSnr /= snrEventCount;
  printf("average HPol snr for this run was %f \n", avgSnr);
  // output the results to a file
  if (fileOutput) {
    AnitaEventSummary thisEventSum;
    RawAnitaHeader thisHeader;
    Adu5Pat thisGps;
    int thisCircPol;
    float thisPNu;
    float thisWeight;
    float thisTruthTheta;
    float thisTruthPhi;
    float thisMaxPwVal[2][2];
    int thisQCuts0[32];
    int thisQCuts1[32];
    char outputFilename[1024];
    sprintf(outputFilename, "%s/analyzerResults_run%03i_%i.root", outputDir, runNumber, startEntryNum);
    cout << "output file name is " << outputFilename << endl;
    TFile* outFile = new TFile(outputFilename, "RECREATE");
    TTree* resultTree = new TTree("resultTree", "analysis results");
    resultTree->Branch("eventSummary", "AnitaEventSummary", &thisEventSum);
    resultTree->Branch("header", "RawAnitaHeader", &thisHeader);
    resultTree->Branch("gps", "Adu5Pat", &thisGps);  
    resultTree->Branch("circPol", &thisCircPol, "circPol/I");  
    resultTree->Branch("pNu", &thisPNu, "pNu/F");  
    resultTree->Branch("weight", &thisWeight, "weight/F");
    resultTree->Branch("truthPhi", &thisTruthPhi, "truthPhi/F");
    resultTree->Branch("truthTheta", &thisTruthTheta, "truthTheta/F");  
    resultTree->Branch("qCuts0", thisQCuts0, "qCuts0[32]/I");  
    resultTree->Branch("qCuts1", thisQCuts1, "qCuts1[32]/I");  
    resultTree->Branch("maxPwVal", thisMaxPwVal, "maxPwVal[2][2]/F");

    TTree* parmTree = new TTree("parmTree", "analysis parameters");
    parmTree->Branch("cmdLine", &cmdLine, "cmdLine/C[4096]");
    parmTree->Branch("minTrigPhis", &minTrigPhisectors, "minTrigPhis/I");
    parmTree->Branch("maxTrigPhis", &maxTrigPhisectors, "maxTrigPhis/I");
    parmTree->Branch("surfSatThresh", &surfSatThreshold, "surfSatThresh/F");
    parmTree->Branch("dcThresh", &dcOffsetThreshold, "dcThresh/F");
    parmTree->Branch("minWfLength", &minWfLength, "minWfLength/I");
    parmTree->Branch("nadirThresh", &nadirNoiseThreshold, "nadirThresh/F");
    parmTree->Branch("filterLnPol", filterDescLn, "filterLnPol/C");
    parmTree->Branch("filterCrPol", filterDescCr, "filterCrPol/C");
    parmTree->Branch("simMode", &simulationMode, "simMode/O");

    parmTree->Fill();
    for (int k=0; k<eventSummaryList.size(); ++k) {
      thisEventSum = eventSummaryList[k];
      thisHeader = headerList[k];
      thisGps = gpsList[k];
      thisCircPol = circPolList[k];
      for (int j=0; j<2; ++j) {
        for (int p=0; p<2; ++p) {
          thisMaxPwVal[p][j] = maxPwValList[k][p][j];
        }
      }
      thisPNu = 0;
      thisWeight = 0;
      thisTruthPhi = 0;
      thisTruthTheta = 0;
      if (simulationMode) {
        thisPNu = pNuList[k];
        thisWeight = pWgtList[k];
        thisTruthPhi = truthPhiList[k];
        thisTruthTheta = truthThetaList[k];
      }
      cout << "  " << thisEventSum.eventNumber << " " << thisWeight << " " << thisTruthPhi << " " << thisTruthTheta << endl; 
      
      vbprintf("event %i \n", thisEventSum.eventNumber);
      vbprintf("  quality cuts pol 0: ");
      
      for (int j=0; j<32; ++j) {
        thisQCuts0[j] = qCuts0[k][j]; 
        vbprintf("%i ", thisQCuts0[j]);
      }
      vbprintf("\n");
      vbprintf("  quality cuts pol 1: ");
      for (int j=0; j<32; ++j) { 
        thisQCuts1[j] = qCuts1[k][j];
        vbprintf(" %i", thisQCuts1[j]);
      }
      vbprintf("\n");
      resultTree->Fill();      
    }
    parmTree->Write();
    resultTree->Write();
    outFile->Save();
    outFile->Close();

  }
  // write out cut counts
  for (int p=0; p<4; ++p) {
    char pChar;
    switch (p) {
     case 0 : pChar = 'H'; break;
     case 1 : pChar = 'V'; break;
     case 2 : pChar = 'L'; break;
     case 3 : pChar = 'R'; break;
    }
    cout << "Polarization " << pChar << endl;
    cout << " Quality cuts: " << endl;
    cout << "   trigger type (non-RF):         " << qualityCutCounts[p][Q_CUT_TRIGGER_TYPE] << endl;
    cout << "   no triggering phisectors:      " << qualityCutCounts[p][Q_CUT_NO_TRIGGER] << endl;
    cout << "   too few triggering phisectors: " << qualityCutCounts[p][Q_CUT_NUM_PHISECTOR] << endl;
    cout << "   payload blast:                 " << qualityCutCounts[p][Q_CUT_PAYLOAD_BLAST] << endl;
    cout << "   waveform saturation:           " << qualityCutCounts[p][Q_CUT_WF_SATURATION] << endl;
    cout << "   waveform DC offset:            " << qualityCutCounts[p][Q_CUT_WF_DC_OFFSET] << endl;
    cout << "   short waveforms                " << qualityCutCounts[p][Q_CUT_SHORT_WF] << endl;
    cout << "   nadir noise                    " << qualityCutCounts[p][Q_CUT_NADIR_NOISE] << endl;
  }    
  
  FFTtools::saveWisdom("wisdom.dat");
  cout << "Goodbye world" << endl;
  fflush(stdout);
  return 0;
}

int processKeywordParm(char* kwp) {
  int result = 0;
  stringstream thisParm(kwp);
  string keyword, parm;
  getline(thisParm, keyword, '=');
  getline(thisParm, parm, '=');
  //thisParm >> keyword;
  //transform(keyword.begin(), keyword.end(), keyword.begin(), ::tolower);
  if (keyword.compare("--MIN_TRIG_PHISECTORS")==0) {minTrigPhisectors=stoi(parm);}
  // TODO INPUT_FILE or INPUT_DIR   use {strcpy(inputFilePath, parm.c_str());}
  else if (keyword.compare("--MAX_TRIG_PHISECTORS")==0) {maxTrigPhisectors=stoi(parm);}
  else if (keyword.compare("--SURF_SAT_THRESHOLD")==0) {surfSatThreshold = stof(parm);}
  else if (keyword.compare("--DC_OFFSET_THRESHOLD")==0) {dcOffsetThreshold = stof(parm);}
  else if (keyword.compare("--MIN_WF_LENGTH")==0) {minWfLength = stoi(parm);}
  else if (keyword.compare("--NADIR_NOISE_THRESHOLD")==0) {nadirNoiseThreshold = stof(parm);}
  else if (keyword.compare("--SINE_SUBTRACT_THRESHOLD")==0) {sineSubThreshold = stof(parm);}
  else if (keyword.compare("--TIME_DOMAIN_PAD_FACTOR")==0) {timeDomainPadFactor = stoi(parm);}
  //else if (keyword.compare("--USE_ABBYS_FILTERING")==0) {filterOption = FILTER_OPTION_ABBY;}
  else if (keyword.compare("--FILTER_OPTION")==0) {filterOption=stoi(parm);}
  else if (keyword.compare("--ADAPTIVE_NOTCH_FILTER_ORDER")==0) {adaptiveNotchFilterOrder=stoi(parm);}
  else if (keyword.compare("--ADAPTIVE_NOTCH_FILTER_WIDTH")==0) {adaptiveNotchFilterWidth=stof(parm);}
  else if (keyword.compare("--ADAPTIVE_NOTCH_FILTER_CAUSAL")==0) {adaptiveNotchFilterCausal=(stoi(parm)==0);}
  //else if (keyword.compare("--SINE_SUBTRACT_FREQ_BANDS")==0) {sineSubFreqBands = stoi(parm);}
  else result = -1;
  cout << "keyword parameter " << kwp << " result code " << result << endl;
  fflush(stdout);
  return result;
}


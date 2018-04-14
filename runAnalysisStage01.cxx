// TODO separate cut application and plotting from this program
//   make damned sure that ALL cut variables get stored to the database
//   always save all events to the database regardless of cut values
//   then start moving the cut application and plotting code to another program

// calculates analysis cuts allowing command line --keywords for cut parameter values
//   see processKeywordParm method in this source to see what the keywords are

// gets quality cuts from bool array in input ROOT tree

// cut flags:
//  applies q/a cuts based on user command line -options: -<flag><type><code> use as many as you want
//  this affects the plots AND the output data trees, but NOT the cut table 
//    the printed cut table is ALWAYS made using the "cut order" () regardless of the cut flags 
//    work with variables aCutOrderStage1 and qCutOrder to control the cut table
//   <flag>:
//     C  accept only events that pass cut
//     c  accept only events that fail cut
//   <type>
//     A  enforce analysis cut
//     Q  enforce quality cut
//     a  unenforce analysis cut
//     q  unenforce quality cut
//   ### deprecated <pol>
//   ### deprecated   0  first polarization (H or L depending on doCircPol)
//   ### deprecated   1  second polarization (V or R depending on doCircPol)
//   <code>  
//     see the command arg processor near the beginning of main() for codes
//     if omitted, cuts in the cut order are applied include-passing-only
//     not case-sensitive
//   e.g.:  -CA         enforce all analysis cuts as defined by aCutOrderStage1
//          -CAcpeak    means include only events that pass the correlation peak cut
//          -Casun      means unenforce the sun cut 
//          -CA -caWAIS enforce all analysis cuts in aCutOrderStage1, except for the WAIS cut
//
// plots analysis cut variables and cut lines
// plots sky maps with cuts applied

// TODO plots of cPol cut variables:
    //  peak separation vs SNR (2D w/counts on z-axis)
    //  L/R relative peak height vs SNR (2D w/counts on z-axis)
    //  peak separation vs rotated cut y-intercept
    //  L/R relative peak height vs rotated cut y-intercept
    
    //  fix calculation of peak separation DONE
    //  1D distance between peaks  DONE   
    //  1D lesser of two cPol peaks DONE NOT TESTED
    //  2D both circPol peaks DONE NOT TESTED
    //  2D lesser cPol peak and dist between peaks DONE NOT TESTED
    //  2D both linPol peaks ALREADY THERE?

#include <iostream>
#include <sstream>
#include <iomanip>
#include "TChain.h"
#include "TTree.h"
#include "TApplication.h"
#include "TH2I.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TFile.h"
#include "TStyle.h"
#include "TEllipse.h"
#include "TPolyLine.h"
#include "TMarker.h"
#include "TGaxis.h"
#include "TF1.h"
#include "TROOT.h"
#include "TDirectory.h" 
#include "TText.h"
#include "TTimeStamp.h"

#include "RawAnitaHeader.h"
#include "AnitaEventSummary.h"
#include "Adu5Pat.h"
#include "UsefulAdu5Pat.h"
//#include "WaveformInfo.h"
#include "AnitaConventions.h"
#include "BedmapReader.h"
#include "AnitaGeomTool.h"
#include "UCUtil.h"
//#include "TurfRate.h"

#include <healpix_base.h>
#include "InterfUtil.h"
#include "analysisCuts.h"

#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_RESET   "\x1b[0m"
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_STYLE_BOLD    "\033[1m"
#define ANSI_STYLE_NOBOLD  "\033[22m"

using namespace std;

int canANITASeeStripe(double lon_anita,double lat_anita,double lat_satellite,double lon_satellite,double height_anita,double height_satellite);

int processKeywordParm(char* kwp);
vector<float> surfaceEllipseParms0(float theta, float dist, float alt, float errPh, float errTh);
vector<float> surfaceEllipseParms1(float phi, float theta, float dist, float errPh, float errTh,
        UsefulAdu5Pat* gps, BedmapReader* bedmap);

int aCutOrderStage1Limit = aCutOrderStage1.size();
int polFracMethod = 0;

vector<vector<double> > conMapPoints = {MCMURDO, VOSTOK, WAIS_DIVIDE, SOUTH_POLE, HALLEY, DAVIS, HOT_SPOT_1};
int main(int argc, char* argv[]) {

  printf("Hello world %s \n", argv[0]);
  aCutOrderStage1 =    {
                          //A_CUT_PEAK_RATIO,
                          //A_CUT_CORR_PEAK,
                          //A_CUT_HILBERT_PEAK,
                          //A_CUT_LINPOL_FRAC,
                          A_CUT_SOLAR_REFL,
                          A_CUT_CONTINENT,
                          //A_CUT_DIAGONAL,
                          A_CUT_THETA,
                          A_CUT_L3_TRIG_DIR,
                          //A_CUT_DEADTIME,
                          //A_CUT_CPOL_PEAK_SEPARATION,
                          //A_CUT_CPOL_PEAK_STRENGTH,
                          //A_CUT_SOLAR_DIRECT,
                          //A_CUT_HW_ANGLE,
                          //A_CUT_GEOSTAT_DIRECT,
                          //A_CUT_GEOSTAT_REFL,
                          //A_CUT_CIRCPOL_FRAC,
                          //A_CUT_ANYPOL_FRAC,
                          //A_CUT_SNR,
                          A_CUT_WAIS,
                          A_CUT_LDB,
                          //A_CUT_HICAL,
                          //A_CUT_SOUTH
                        };
  //TCanvas* dummyCanv = new TCanvas("dummy", "", 0,0);
  
  bool doCircPol = true;  // this enables the processing of the CPol entry immediately following the linear pol
  int dbStride = 1;
  int maxEvents = 0;
  int startRunNumber = 343;
  int endRunNumber = 343;
  int singleEventNumber = 0;
  int qualityCutMask[2][NUM_Q_CUTS] = {{0}}; for (int p=0; p<2; ++p) for (int k : qCutOrder) qualityCutMask[p][k]=1;
  int preQualityCutMask[2][NUM_PRE_Q_CUTS] = {{0}}; for (int p=0; p<2; ++p) for (int k : preQCutOrder) preQualityCutMask[p][k]=1;
  int analysisCutMask[2][NUM_A_CUTS] = {{0}};
  vector<int> whichPols = {0, 1}; 
  bool listPassingEvents = false;
  //vector<double> ksLatLon = WAIS_DIVIDE;
  bool simulationMode = false;
  bool displayGraphics = false;
  bool saveOutput = true;
  bool useEventList = false;
  //int eventNumber = 0;
  int pNum =0;
  //char basisStr[4] = "H/V";
  char pChars[3]="HV";
  char pCharsC[3]="LR";
  bool useTrigNsWindow = false;
  int dataIndex = 0;
  bool use90data = false;
  
  //int ourHpBin = 3052;

  char* inputDirStr = getenv("ANITA_ANALYSIS_RESULTS_DIR");
  char inputDirName[1024];
  strcpy(inputDirName, inputDirStr);

  char outputDir[1024] = "results/plots";
  char eventListFilename[1024] = "eventList.txt";
  char runDesc[32] = "";
  char calPulser[16] = "";
  for (int i=1; i<argc; ++i) {
    //vcout << "parameter " << i << " of " << argc << ": " << argv[i] << endl;
    if (argv[i][0] == '-') {
      if (strlen(argv[i]) > 1) {
        if (argv[i][1] == '-') {processKeywordParm(argv[i]);}
        else if (argv[i][1] == 'V' || argv[i][1] == 'v') {verboseMode = true;}
        else if (argv[i][1] == 'T' || argv[i][1] == 't') {
          useTrigNsWindow = true;
          sprintf(calPulser, "%s", "WAIS");
          if (strlen(argv[i])>2) {
            char thisParm[1024];
            for (unsigned int k=2; k<=strlen(argv[i]); k++) {
              thisParm[k-2] = toupper(argv[i][k]);
            }
            sprintf(calPulser, "%s", thisParm);
          }          
        }      // cal pulses
        else if (argv[i][1] == '9')                      { use90data = true; }
        else if (argv[i][1] == 'L' || argv[i][1] == 'l') { listPassingEvents = true; }
        else if (argv[i][1] == 'U' || argv[i][1] == 'u') { simulationMode = true; }
        else if (argv[i][1] == 'G' || argv[i][1] == 'g') { displayGraphics = true; }
        else if (argv[i][1] == 'E' || argv[i][1] == 'e') {
          useEventList = true;
          if (strlen(argv[i]) > 2) {
            //char thisParm[8];
            unsigned int k=2;
            for (; k<=strlen(argv[i]); k++) {
              eventListFilename[k-2] = argv[i][k];
            }
            eventListFilename[k-2] = 0;
          }                      
        }
        else if (argv[i][1] == 'R' || argv[i][1] == 'r') {
          if (strlen(argv[i]) > 2) {
            //char thisParm[8];
            unsigned int k=2;
            for (; k<=strlen(argv[i]); k++) {
              runDesc[k-2] = argv[i][k];
            }
            runDesc[k-2] = 0;
          }            
        } else if (argv[i][1] == 'D' || argv[i][1] == 'd') {
          if (strlen(argv[i]) > 2) {
            //char thisParm[8];
            unsigned int k=2;
            for (; k<=strlen(argv[i]); k++) {
              inputDirName[k-2] = argv[i][k];
            }
            inputDirName[k-2] = 0;
          }            
        } else if (argv[i][1] == 'P' || argv[i][1] == 'p') { 
          if ((strlen(argv[i])>2) && ((argv[i][2]=='c') || (argv[i][2]=='C'))) {
            doCircPol = true; 
            //sprintf(basisStr, "L/R"); // just for titles
            //sprintf(pChars, "LR");     // just for titles
          } 
          if (strlen(argv[i])>3) {
            if (argv[i][3] == '0') {
              whichPols = {0};
            } else if (argv[i][3] == '1') {
              whichPols = {1};
            }
          }
        } else if (argv[i][1] == 'M' || argv[i][1] == 'm') {
          if (strlen(argv[i]) > 2) {
            char thisParm[8];
            for (unsigned int k=2; k<=strlen(argv[i]); k++) {
              thisParm[k-2] = argv[i][k];
            }
            maxEvents = atoi(thisParm);
          }            
        } else if (argv[i][1] == 'S' || argv[i][1] == 's') {
          if (strlen(argv[i]) > 2) {
            char thisParm[8];
            for (unsigned int k=2; k<=strlen(argv[i]); k++) {
              thisParm[k-2] = argv[i][k];
            }
            dbStride = atoi(thisParm);
          }
        } else if (argv[i][1] == 'i' || argv[i][1] == 'I') {
          if (strlen(argv[i]) > 2) {
            char thisParm[8];
            for (unsigned int k=2; k<=strlen(argv[i]); k++) {
              thisParm[k-2] = argv[i][k];
            }
            dataIndex = atoi(thisParm);
          }

        } else if (argv[i][1] == 'O' || argv[i][1]=='o') {
          saveOutput = true;
          if (strlen(argv[i]) > 2) {
            unsigned int k=2;
            for (; k<=strlen(argv[i]); k++) {
              outputDir[k-2] = argv[i][k];
            }
            outputDir[k-2] = 0;
          }                    
        } else if (argv[i][1] == 'C' || 'c') {
          int cutSense = (argv[i][1]=='c') ? -1 : 1;
          if (strlen(argv[i]) > 2) {
            //int polNum = argv[i][3] - '0';
            if (strlen(argv[i]) > 3) {
              char thisParm[8];
              for (unsigned int k=3; k<=strlen(argv[i]); k++) {
                thisParm[k-3] = toupper(argv[i][k]);
              }
              vbprintf("cut parm is %s \n", thisParm);
              if (argv[i][2] == 'a' || argv[i][2] == 'A') {
                if (argv[i][2] == 'a') {cutSense = 0;}
                if        (strcmp(thisParm, "SUN")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_SOLAR_DIRECT] = cutSense;}}
                else if   (strcmp(thisParm, "SUNREFL")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_SOLAR_REFL] = cutSense;}}
                else if   (strcmp(thisParm, "GS")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_GEOSTAT_DIRECT] = cutSense;}}
                else if   (strcmp(thisParm, "GSREFL")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_GEOSTAT_REFL] = cutSense;}}
                else if   (strcmp(thisParm, "CONT")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_CONTINENT] = cutSense;}}
                else if   (strcmp(thisParm, "THETA")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_THETA] = cutSense;}}
                else if   (strcmp(thisParm, "HWA")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_HW_ANGLE] = cutSense;}}
                else if   (strcmp(thisParm, "PRATIO")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_PEAK_RATIO] = cutSense;}}
                else if   (strcmp(thisParm, "LINPOL")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_LINPOL_FRAC] = cutSense;}}
                else if   (strcmp(thisParm, "CIRCPOL")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_CIRCPOL_FRAC] = cutSense;}}
                else if   (strcmp(thisParm, "ANYPOL")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_ANYPOL_FRAC] = cutSense;}}
                else if   (strcmp(thisParm, "CPEAK")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_CORR_PEAK] = cutSense;}}
                else if   (strcmp(thisParm, "HPEAK")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_HILBERT_PEAK] = cutSense;}}
                else if   (strcmp(thisParm, "DIAG")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_DIAGONAL] = cutSense;}}
                else if   (strcmp(thisParm, "WAIS")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_WAIS] = cutSense;}}
                else if   (strcmp(thisParm, "LDB")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_LDB] = cutSense;}}
                else if   (strcmp(thisParm, "HICAL")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_HICAL] = cutSense;}}
                else if   (strcmp(thisParm, "SOUTH")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_SOUTH] = cutSense;}}
                else if   (strcmp(thisParm, "SNR")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_SNR] = cutSense;}}
                else if   (strcmp(thisParm, "WEDGE")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_WEDGE] = cutSense;}}
                else if   (strcmp(thisParm, "RECT")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_RECT] = cutSense;}}
                else if   (strcmp(thisParm, "CONTCIRC")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_CONT_CIRC] = cutSense;}}
                else if   (strcmp(thisParm, "SKYCIRC")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_SKY_CIRC] = cutSense;}}
                else if   (strcmp(thisParm, "L3TRIGDIR")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_L3_TRIG_DIR] = cutSense;}}
                else if   (strcmp(thisParm, "ENU")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_SIM_E_NU] = cutSense;}}
                else if   (strcmp(thisParm, "CPEAKSEP")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_CPOL_PEAK_SEPARATION] = cutSense;}}
                else if   (strcmp(thisParm, "CPEAKLESSER")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_CPOL_LESSER_PEAK] = cutSense;}}
                else if   (strcmp(thisParm, "DEADTIME")==0) {for (int p=0; p<2; ++p) {analysisCutMask[p][A_CUT_DEADTIME] = cutSense;}}
              } else if (argv[i][2] == 'q' ||argv[i][2] == 'Q') {
                if (argv[i][2] == 'q') {cutSense = 0;}
                if        (strcmp(thisParm, "TRIGTYPE")==0) {for (int p=0; p<2; ++p) {qualityCutMask[p][Q_CUT_TRIGGER_TYPE] = cutSense;}}
                else if   (strcmp(thisParm, "NPHI")==0) {for (int p=0; p<2; ++p) {qualityCutMask[p][Q_CUT_NUM_PHISECTOR] = cutSense;}}
                else if   (strcmp(thisParm, "PLB")==0) {for (int p=0; p<2; ++p) {qualityCutMask[p][Q_CUT_PAYLOAD_BLAST] = cutSense;}}
                else if   (strcmp(thisParm, "SAT")==0) {for (int p=0; p<2; ++p) {qualityCutMask[p][Q_CUT_WF_SATURATION] = cutSense;}}
                else if   (strcmp(thisParm, "DCO")==0) {for (int p=0; p<2; ++p) {qualityCutMask[p][Q_CUT_WF_DC_OFFSET] = cutSense;}}
                else if   (strcmp(thisParm, "SHORTWF")==0) {for (int p=0; p<2; ++p) {qualityCutMask[p][Q_CUT_SHORT_WF] = cutSense;}}
                else if   (strcmp(thisParm, "NN")==0) {for (int p=0; p<2; ++p) {qualityCutMask[p][Q_CUT_NADIR_NOISE] = cutSense;}}
                else if   (strcmp(thisParm, "NOTRIG")==0) {for (int p=0; p<2; ++p) {qualityCutMask[p][Q_CUT_NO_TRIGGER] = cutSense;}}
              }               // TODO else action here?
            } else {
              if (argv[i][2]=='A' || argv[i][2]=='a') {
                if (argv[i][2] == 'a') {cutSense = 0;}
                for (int k : aCutOrderStage1) {analysisCutMask[0][k] =  cutSense; analysisCutMask[1][k] =  cutSense;}
              }
              else if (argv[i][2]=='Q' || argv[i][2]=='q') {
                if (argv[i][2] == 'q') {cutSense = 0;}
                for (int k : qCutOrder) {qualityCutMask[0][k] = cutSense; qualityCutMask[1][k] = cutSense;} 
              }
            }
          } else {            
            // set all cuts on
            for (int k : qCutOrder) {
              qualityCutMask[0][k] = cutSense;
              qualityCutMask[1][k] = cutSense;
            }
            for (int k : aCutOrderStage1) {
              analysisCutMask[0][k] = cutSense;
              analysisCutMask[1][k] = cutSense;
            }
          }
        }
      }
    } else {
      ++pNum;
      if (pNum==1) {startRunNumber = atoi(argv[i]); endRunNumber = startRunNumber;}
      else if (pNum==2) {endRunNumber = atoi(argv[i]);}
      else if (pNum==3) {
        singleEventNumber = atoi(argv[i]);
        maxEvents = 1;
      }
    }
    //vcout << "i=" << i << endl;
  }
  
  TApplication* app = 0;
  if (displayGraphics) {app = new TApplication("app", &argc, argv);}
  
  char logFilename[1024]; sprintf(logFilename, "%s/analyzerResultsIterator06_%i_%i.log", outputDir, startRunNumber, endRunNumber);
  printf("log filename is %s \n", logFilename);
  logFile = fopen(logFilename, "w");
  //fprintf(logFile, "test");
  for (int i=0; i<argc; ++i) {
    lprintf("%s ", argv[i]);
  }
  lprintf("\n");  
  lprintf(" compiled %s %s \n ", __DATE__, __TIME__);
  lprintf(" run %s UTC\n", (new TTimeStamp())->AsString("s"));
  lprintf("input dir is %s \n", inputDirName);
  lprintf("output dir is %s \n", outputDir);

  aCutOrderStage1Limit = aCutOrderStage1.size();  // in case default cut order is overridden

  // set up palette for localization pdf markers
  int pdfPalette[100];
  // blue-to-yellow
  //double pdfR[3] = {0, 1, 1};
  //double pdfG[3] = {0, 1, 1};
  //double pdfB[3] = {1, 0, 0};
  // blackl-to-white
  double pdfR[3] = {1, 0, 0};
  double pdfG[3] = {1, 0, 0};
  double pdfB[3] = {1, 0, 0};
  double pdfStops[3] = {0, 1, 1};
  int pdfP = TColor::CreateGradientColorTable(3, pdfStops, pdfR, pdfG, pdfB, 100);
  for (int k=0; k<100; ++k) {pdfPalette[k] = pdfP + k;}

  //if (dbStride%2==0) ++dbStride;  
  int numPols = whichPols.size();
  float ksLon = AnitaLocations::LONGITUDE_WAIS;
  float ksLat = AnitaLocations::LATITUDE_WAIS;
  float ksAlt = AnitaLocations::ALTITUDE_WAIS;
  int pulserNsTime = 0;
  if (strcmp(calPulser, "WAIS")==0) {
    //ksLon = AnitaLocations::LONGITUDE_WAIS;
    //ksLat = AnitaLocations::LATITUDE_WAIS;
    //ksAlt = AnitaLocations::ALTITUDE_WAIS;
    //pulserNsTime = pulserNsWais;
  } 
  else if (strcmp(calPulser, "LDB")==0) {
    ksLon = AnitaLocations::LONGITUDE_LDB;
    ksLat = AnitaLocations::LATITUDE_LDB;
    ksAlt = AnitaLocations::ALTITUDE_LDB;
    //pulserNsTime = pulserNsLdb;  //LDB seavey hPol 10kV
  }

  if (!displayGraphics) gROOT->SetBatch(kTRUE);
  vbprintf("verbose mode \n");
  lprintf(" doCircPol flag (0) is %i \n", doCircPol);
  lprintf("Run parameters\n");
  if (doCircPol) {lprintf(" R/L polarization basis \n");}
  else {lprintf(" H/V polarization basis \n");}
  lprintf(" Run numbers                          %3i - %3i\n", startRunNumber, endRunNumber);
  lprintf(" Sun event cut distance (phi)         %.2f\n", sunCutDistPhi);
  lprintf(" Sun event cut distance (theta)       %.2f\n", sunCutDistTheta);
  lprintf(" Sun reflection cut distance (phi)    %.2f\n", sunReflCutDistPhi);
  lprintf(" Sun reflection cut distance (theta)  %.2f\n", sunReflCutDistTheta);
  lprintf(" Minimum triggering phi-sectors       %i\n", minTrigPhisectors);
  lprintf(" Maximum triggering phi-sectors       %i\n", maxTrigPhisectors);
  lprintf(" Surf saturation threshold            %.0f\n", surfSatThreshold);
  lprintf(" DC offset threshold                  %.2f\n", dcOffsetThreshold);
  lprintf(" Minimum waveform sample length       %i\n", minWfLength);
  lprintf(" Nadir noise threshold                %.3f\n", nadirNoiseThreshold);
  lprintf(" Minimum peak theta                   %.1f\n", minPeakTheta);
  lprintf(" Maximum peak-to-trigger angle        %.1f\n", maxHwAngle);
  lprintf(" Phi-sector trigger proximity         %i\n",   trigPhisecProx);
  lprintf(" Minimum snr                          %.1f\n", minSnr);
  lprintf(" Minimum linear polarization fraction %.3f\n", minLinPolFrac);
  lprintf(" Minimum total polarization fraction  %.3f\n", minAnyPolFrac);
  lprintf(" Maximum circ polarization fraction   %.3f\n", maxCircPolFrac);
  lprintf(" Minimum correlation peak value       %.5f\n", corrPeakThresh);
  lprintf(" Minimum hilbert peak value           %.1f\n", hilbPeakThresh);
  lprintf(" Rotated cut corr peak interecept     %.3f\n", corrPeakDiagThresh);
  lprintf(" Rotated cut hilbert peak intercept   %.1f\n", hilbPeakDiagThresh);
  lprintf(" 2nd/1st peak ratio                   %.3f\n", peakRatioThreshold);
  lprintf(" rectangle cut : %6.3f, %6.3f, %6.3f, %6.3f \n", skyRectCut[0], skyRectCut[1], skyRectCut[2], skyRectCut[3]);
  lprintf(" simulated neutrino energy cut : %6.3f, %6.3f\n", simENuCut[0], simENuCut[1]);
  for (int whichPol=0; whichPol<numPols; ++whichPol) {
    int p = whichPols[whichPol];
    lprintf("  polarization %c \n", pChars[p]);
    //TODO put this to a loop  
    for (int k=0; k<NUM_PRE_Q_CUTS; ++k) {if (preQualityCutMask[p][k] != 0) lprintf("   pre quality cut selected: %s, %i \n", PRE_Q_CUT_DESC[k], preQualityCutMask[p][k]);}
    for (int k=0; k<NUM_Q_CUTS; ++k) {if (qualityCutMask[p][k] != 0) lprintf("   quality cut selected: %s, %i \n", Q_CUT_DESC[k], qualityCutMask[p][k]);}
    for (int k=0; k<NUM_A_CUTS; ++k) {if (analysisCutMask[p][k] != 0) lprintf("   analysis cut selected: %s, %i \n", A_CUT_DESC[k], analysisCutMask[p][k]);}
  }
  lprintf("Analysis cut order: \n");
  //for (int k=0; k<aCutOrderStage1.size(); ++k)  {
  for (int k=0; k<aCutOrderStage1Limit; ++k)  {
    lprintf(" %02i: %s ", aCutOrderStage1[k], A_CUT_DESC[aCutOrderStage1[k]]);
    if      (aCutOrderStage1[k] == A_CUT_CORR_PEAK) { lprintf("cut value = %f ", corrPeakThresh);}
    else if (aCutOrderStage1[k] == A_CUT_HILBERT_PEAK) { lprintf("cut value = %f ", hilbPeakThresh);}
    else if (aCutOrderStage1[k] == A_CUT_DIAGONAL) { lprintf("corrInt = %f, hilbInt = %f ", corrPeakDiagThresh, hilbPeakDiagThresh);}
    else if (aCutOrderStage1[k] == A_CUT_PEAK_RATIO) { lprintf("cut value = %f ", peakRatioThreshold);}
    else if (aCutOrderStage1[k] == A_CUT_LINPOL_FRAC) {lprintf("cut value = %f ", minLinPolFrac);}
    else if (aCutOrderStage1[k] == A_CUT_THETA) {lprintf("min=%f, max=%f ", minPeakTheta, maxPeakTheta);}
    else if (aCutOrderStage1[k] == A_CUT_L3_TRIG_DIR) {lprintf("cut value = %i ", trigPhisecProx);}
    lprintf("\n");
  }
  lprintf("polarization fraction calculation method is %i \n", polFracMethod);
  //return 0;
  BedmapReader* bedmap = BedmapReader::Instance(false);
  double ksEa, ksNo;
  bedmap->LonLattoEaNo(ksLon, ksLat, ksEa, ksNo);

  lprintf("10 current dir is %s \n", gDirectory->GetName());
  TFile* coastGraphFile = new TFile("antarcticCoastGraph.root");
  TGraph* coastLineGr = 0;
  coastLineGr = (TGraph*) coastGraphFile->Get("Graph");
  lprintf("Antarctic coastline graph contains %i points \n", coastLineGr->GetN());
  int coastStride = 40;
  for (int k=1; k<coastLineGr->GetN()-1; ++k) {for (int j=0; j<coastStride; ++j) {coastLineGr->RemovePoint(k);}}
  coastLineGr->SetPoint(coastLineGr->GetN(), coastLineGr->GetX()[0], coastLineGr->GetY()[0]);
  for (int k=0; k<coastLineGr->GetN(); ++k) {coastLineGr->SetPoint(k, coastLineGr->GetX()[k]/1000, coastLineGr->GetY()[k]/1000);}
  lprintf("Antarctic coastline graph contains %i points \n", coastLineGr->GetN());
  coastGraphFile->Close();
  fflush(stdout);
  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  /*
  TChain* filterPowerFracTree = 0;
  float filterPowerRatio;
  //TFile* filterPowerFracFile = 0;
  
  //filterPowerFracFile = new TFile("results/filterPowerFrac_WAIS_1000.root");
  filterPowerFracTree = new TChain("filterPowerLossTree");
  filterPowerFracTree->Add("results/filterPowerFrac_WAIS_1000.root");
  filterPowerFracTree->BuildIndex("eventNumber", "pol");
  filterPowerFracTree->SetBranchAddress("powerRatio", &filterPowerRatio);
  */
  Healpix_Ordering_Scheme scheme = Healpix_Ordering_Scheme::RING;
  T_Healpix_Base<int>* healpix = new T_Healpix_Base<int>(healPixN, scheme);
  lprintf("healpix map has %i pixels \n", healpix->Npix());

  RawAnitaHeader* header = 0;
  AnitaEventSummary* eventSummary = 0;
  Adu5Pat* gpsRaw = 0;
  //TurfRate* turfRate = 0;
  int circPol = 0;
  float pnu = 0;
  float weight = 0;
  float truthPhi = 0;
  float truthTheta = 0; 
  int qCuts[2][NUM_Q_CUTS] = {{0}};
  int preQCuts[2][NUM_PRE_Q_CUTS] = {{0}};
  RawAnitaHeader* cPolHeader = 0;
  AnitaEventSummary* cPolEventSummary = 0;
  float maxPwVal[2][2];
  
  float errPhiBySnr[2][20] = {{0}};
  float errThetaBySnr[2][20] = {{0}};
  float eventCountsSnr[2][20] = {{0}};
  float acceptsSnr[2][20] = {{0}};
  
  gStyle->SetLabelSize(0.05, "xyz");
  //gStyle->SetPalette(kRainBow, 0);
  
  TChain* resultTree = new TChain("resultTree");
  //TChain* turfRateTree = new TChain("turfRateTree");
  TChain* parmTree = new TChain("parmTree");
  parmTree->SetBranchAddress("minTrigPhis", &minTrigPhisectors);
  parmTree->SetBranchAddress("maxTrigPhis", &maxTrigPhisectors);
  parmTree->SetBranchAddress("surfSatThresh", &surfSatThreshold);
  parmTree->SetBranchAddress("dcThresh", &dcOffsetThreshold);
  parmTree->SetBranchAddress("minWfLength", &minWfLength);
  parmTree->SetBranchAddress("nadirThresh", &nadirNoiseThreshold);


  // parameters output to the next stage in the root file.
  int oRunNumber;
  int oEventNumber;
  char oPol;
  //float oEllParms[4];
  float oHPeak;
  float oCPeak;
  float oCohSnr;
  float oCohSnr1;
  float oCohSnr2;
  float oCohSnr3;
  float oMaxSnr;
  int oHpBin;
  float oHpWeight;
  float oPfFrac_0;
  float oPfFrac_1;
  float oLinPolFrac;
  float oCalPulseDist;
  float oPeakRatio;
  float oLdCutVal;
  float oCPolDist;
  float oCPolLesserPeak;
  float oCPolPeak[2];
  float oEventWeight;
  //float oMaxPwVal[2][2];
  float oMinCircPeakVal;
  float oDeadTimeFrac;
  //int oQCuts[NUM_Q_CUTS];
  int oACuts[NUM_A_CUTS];

  // pos info params
  //float oEllParms[4];
  float oTheta;
  float oPhi;
  float oEa;
  float oNo;
  float oLat;
  float oLon;
  float oAlt;
  Adu5Pat* oGpsRaw = 0;
  UsefulAdu5Pat* oGps = 0;

  // me being smart for oindree
  float oLTheta;
  float oRTheta;
  float oLPhi;
  float oRPhi;
  float oLSnr2;
  float oRSnr2;
  //float oLCPeak;
  //float oRCPeak;

  float oTruthPhi;
  float oTruthTheta;

  char outputFilename[1024];
  if (use90data) {
    sprintf(outputFilename, "%s/analysisOutput_%i_%i_%i.root", outputDir, startRunNumber, endRunNumber, dataIndex);
  } else {
    sprintf(outputFilename, "%s/analysisOutput_%i_%i.root", outputDir, startRunNumber, endRunNumber);
  }
   
  TFile* outputFile = new TFile(outputFilename, "RECREATE");
  outputFile->cd();
  lprintf ("20 current directory is %s \n", gDirectory->GetName());


  //char oQCutsDef[32]; sprintf(oQCutsDef, "qCuts[%i]/I", NUM_Q_CUTS);
  char oACutsDef[32]; sprintf(oACutsDef, "aCuts[%i]/I", NUM_A_CUTS);
  TTree* outputTree0 = new TTree("analysisOutput0", "analysis0");
  TTree* outputTree1 = new TTree("analysisOutput1", "analysis1");
  outputTree0->Branch("runNumber", &oRunNumber, "runNumber/I");
  outputTree0->Branch("eventNumber", &oEventNumber, "eventNumber/I");
  outputTree0->Branch("pol", &oPol, "pol/B");                                         
  outputTree0->Branch("eventWeight", &oEventWeight, "eventWeight/F");
  outputTree0->Branch("hPeak", &oHPeak, "hPeak/F");                                       // coherent reconstruction hilbert peak value
  outputTree0->Branch("cPeak", &oCPeak, "cPeak/F");                                       // correlation map peak value
  outputTree0->Branch("cohSnr", &oCohSnr, "cohSnr/F");                                    // coherent recon SNR methid 0
  outputTree0->Branch("cohSnr1", &oCohSnr1, "cohSnr1/F");                                 // coherent recon SNR methid 1
  outputTree0->Branch("cohSnr2", &oCohSnr2, "cohSnr2/F");                                 // coherent recon SNR methid 2
  outputTree0->Branch("cohSnr3", &oCohSnr3, "cohSnr3/F");                                 // coherent recon SNR methid 3
  outputTree0->Branch("maxWfSnr", &oMaxSnr, "maxWfSnr/F");                                // SNR of maximum waveform (method 0)
  //outputTree0->Branch("ellParms", &oEllParms, "ellParms[4]/F");                           // error ellipse parameters
  outputTree0->Branch("pfFrac_0", &oPfFrac_0, "pfFrac_0/F");                              // filtered power fraction 
  outputTree0->Branch("pfFrac_1", &oPfFrac_1, "pfFrac_1/F");                              // filtered power fraction 
  outputTree0->Branch("linPolFrac", &oLinPolFrac, "linPolFrac/F");                        // linear polarization fraction
  outputTree0->Branch("deadTimeFrac", &oDeadTimeFrac, "deadTimeFrac/F");                  // deadtime fraction 
  outputTree0->Branch("peakRatio", &oPeakRatio, "peakRatio/F");                           // ratio of 2nd to 1st peak
  outputTree0->Branch("ldCutVal", &oLdCutVal, "ldCutVal/F");                              // linear discriminant (SNR and Corr peak)
  outputTree0->Branch("calPulseDist", &oCalPulseDist, "calPulseDist/F");                  // distance to cal pulser?
  outputTree0->Branch("cPolDist", &oCPolDist, "cPolDist/F");                              // separation between CPol (L/R) peaks
  outputTree0->Branch("cPolLesserPeak", &oCPolLesserPeak, "cPolLesserPeak/F");            // lesser of Lpol & Rpol peaks (normalized)
  outputTree0->Branch("cPolPeak", &oCPolPeak, "cPolPeak[2]/F");                           // circpol peak values
  //outputTree0->Branch("maxPwVal", &oMaxPwVal, "maxPwVal[2][2]/F");
  outputTree0->Branch("minCircPeakVal", &oMinCircPeakVal, "minCircPeakVal/F");            // lesser circ pol val in linPol peak window
  //outputTree0->Branch("qCuts", &oQCuts, oQCutsDef);                                       // quality cut pass/fail array
  outputTree0->Branch("aCuts", &oACuts, oACutsDef);                                       // analysis cut pass/fail array


  //outputTree0->Branch("ellParms", &oEllParms, "ellParms[4]/F");                           // error ellipse parameters
  outputTree0->Branch("theta", &oTheta, "theta/F");                                       // event theta
  outputTree0->Branch("phi", &oPhi, "phi/F");                                             // event phi
  outputTree0->Branch("ea", &oEa, "ea/F");                                                // event easting
  outputTree0->Branch("no", &oNo, "no/F");                                                // event northing
  outputTree0->Branch("lat", &oLat, "lat/F");                                             // event latitude
  outputTree0->Branch("lon", &oLon, "lon/F");                                             // event longitude
  outputTree0->Branch("alt", &oAlt, "alt/F");                                             // event altitude
  outputTree0->Branch("gps", &oGps);                                                      // payload gps info 
  outputTree0->Branch("gpsRaw", &oGpsRaw);                                                // raw payload gps info
  
  outputTree0->Branch("lTheta", &oLTheta, "lTheta/F");
  outputTree0->Branch("lPhi", &oLPhi, "lPhi/F");
  outputTree0->Branch("rTheta", &oRTheta, "rTheta/F");
  outputTree0->Branch("rPhi", &oRPhi, "rPhi/F");
  outputTree0->Branch("lSnr2", &oLSnr2, "lSnr2/F");
  outputTree0->Branch("rSnr2", &oRSnr2, "rSnr2/F");

  
  outputTree1->Branch("runNumber", &oRunNumber, "runNumber/I");
  outputTree1->Branch("eventNumber", &oEventNumber, "eventNumber/I");
  outputTree1->Branch("pol", &oPol, "pol/B");
  outputTree1->Branch("hpBin", &oHpBin, "hpBin/I");
  outputTree1->Branch("hpWeight", &oHpWeight, "hpWeight/F");
  
  // truth sim neutrino direction
  outputTree0->Branch("truthPhi", &oTruthPhi, "truthPhi/F");
  outputTree0->Branch("truthTheta", &oTruthTheta, "truthTheta/F");

  double newLat[2], newLon[2], newAlt[2], newAdj[2];

  int   minTrigPhisectors_0 = 0;
  int   maxTrigPhisectors_0 = 0;
  float surfSatThreshold_0 = 0;
  float dcOffsetThreshold_0 = 0;
  int   minWfLength_0 = 0;
  float nadirNoiseThreshold_0 = 0;
  
  char inputFilename[1024];
  bool parmMismatch = false;

  if (use90data) {
    for (int runNumber = startRunNumber; runNumber <= endRunNumber; ++runNumber) {
      if (true) {
        sprintf(inputFilename, "%s/analyzerResults_run%03i_%i.root", inputDirName, runNumber, dataIndex);
        lprintf("adding input file %s \n", inputFilename);
        resultTree->Add(inputFilename);   
        parmTree->Add(inputFilename);
        parmTree->GetEntry(0);
        if (runNumber == startRunNumber) {
          // retain parameter values from first file
          minTrigPhisectors_0 = minTrigPhisectors;
          maxTrigPhisectors_0 = maxTrigPhisectors;
          surfSatThreshold_0 = surfSatThreshold;
          dcOffsetThreshold_0 = dcOffsetThreshold;
          minWfLength_0 = minWfLength;
          nadirNoiseThreshold_0 = nadirNoiseThreshold;
        } else {
          // compare parameter values to original parms and warn if different
          if (       (minTrigPhisectors_0 != minTrigPhisectors)
                  || (maxTrigPhisectors_0 != maxTrigPhisectors)
                  || (surfSatThreshold_0 != surfSatThreshold)
                  || (dcOffsetThreshold_0 != dcOffsetThreshold)
                  || (minWfLength_0 != minWfLength)
                  || (nadirNoiseThreshold_0 != nadirNoiseThreshold)
                  )  {
            parmMismatch = true;
            lprintf("   mismatch of analyzer run parameters at run %i \n", runNumber);
          }
        }
        //sprintf(inputFilename, "%s/turfRateFile%03i.root", inputDirName, runNumber);
        //turfRateTree->Add(inputFilename);
        fflush(stdout);
      }
      // some (alot) runs were split into smaller parts, deal with them seperatly.
      // -note, alot of these dont have all the runs here, dont worry about missing files too much.
      if (   runNumber == 184 || runNumber == 191 || runNumber == 192 
          || runNumber == 193 || runNumber == 201 || runNumber == 227
          || runNumber == 230 || runNumber == 232 || runNumber == 233
          || runNumber == 237 || runNumber == 237 || runNumber == 238
          || runNumber == 245 || runNumber == 246 || runNumber == 255
          || runNumber == 267 || runNumber == 269 || runNumber == 281
          || runNumber == 287 || runNumber == 292 || runNumber == 293
          || runNumber == 296 || runNumber == 302 || runNumber == 304
          || runNumber == 316 || runNumber == 317 || runNumber == 328
          || runNumber == 330 || runNumber == 341 || runNumber == 361
          || runNumber == 364 || runNumber == 365 || runNumber == 366
          || runNumber == 369 || runNumber == 371 || runNumber == 375
          || runNumber == 382 || runNumber == 387 || runNumber == 388
          || runNumber == 389 || runNumber == 395 || runNumber == 402
          || runNumber == 409 || runNumber == 431 || runNumber == 177
          || runNumber == 179 || runNumber == 180 || runNumber == 181
          || runNumber == 182 || runNumber == 183 || runNumber == 190
          || runNumber == 197 || runNumber == 198 || runNumber == 199
          || runNumber == 204 || runNumber == 205 || runNumber == 208
          || runNumber == 221 || runNumber == 225 || runNumber == 229
          || runNumber == 240 || runNumber == 283 || runNumber == 286
          || runNumber == 419
         ) 
      {
        for (int j = 0; j<10; j++) {
          if ( j ==0 ) {
            sprintf(inputFilename, "%s/run%03i/analyzerResults_run%03i_%i.root", inputDirName, runNumber, runNumber, dataIndex);
          } else {
            sprintf(inputFilename, "%s/run%03i/analyzerResults_run%03i_%i%i.root", inputDirName, runNumber, runNumber, j, dataIndex);
          }
          lprintf("adding input file %s \n", inputFilename);
          resultTree->Add(inputFilename);
          parmTree->Add(inputFilename);
          parmTree->GetEntry(0);
          // compare parameter values to original parms and warn if different
          if (       (minTrigPhisectors_0 != minTrigPhisectors)
                  || (maxTrigPhisectors_0 != maxTrigPhisectors)
                  || (surfSatThreshold_0 != surfSatThreshold)
                  || (dcOffsetThreshold_0 != dcOffsetThreshold)
                  || (minWfLength_0 != minWfLength)
                  || (nadirNoiseThreshold_0 != nadirNoiseThreshold)
                  )  {
            parmMismatch = true;
            lprintf("   mismatch of analyzer run parameters at run %i \n", runNumber);
          }
          //sprintf(inputFilename, "%s/turfRateFile%03i.root", inputDirName, runNumber);
          //turfRateTree->Add(inputFilename);
          fflush(stdout);
        }
      }
    }
  } else {
    for (int runNumber = startRunNumber; runNumber <= endRunNumber; ++runNumber) {
      if (true) {
        sprintf(inputFilename, "%s/analyzerResults_run%03i.root", inputDirName, runNumber);
        lprintf("adding input file %s \n", inputFilename);
        resultTree->Add(inputFilename);
        parmTree->Add(inputFilename);
        parmTree->GetEntry(0);
        if (runNumber == startRunNumber) {
          // retain parameter values from first file
          minTrigPhisectors_0 = minTrigPhisectors;
          maxTrigPhisectors_0 = maxTrigPhisectors;
          surfSatThreshold_0 = surfSatThreshold;
          dcOffsetThreshold_0 = dcOffsetThreshold;
          minWfLength_0 = minWfLength;
          nadirNoiseThreshold_0 = nadirNoiseThreshold;
        } else {
          // compare parameter values to original parms and warn if different
          if (       (minTrigPhisectors_0 != minTrigPhisectors)
                  || (maxTrigPhisectors_0 != maxTrigPhisectors)
                  || (surfSatThreshold_0 != surfSatThreshold)
                  || (dcOffsetThreshold_0 != dcOffsetThreshold)
                  || (minWfLength_0 != minWfLength)
                  || (nadirNoiseThreshold_0 != nadirNoiseThreshold)
                  )  {
            parmMismatch = true;
            lprintf("   mismatch of analyzer run parameters at run %i \n", runNumber);
          }
        }
        //sprintf(inputFilename, "%s/turfRateFile%03i.root", inputDirName, runNumber);
        //turfRateTree->Add(inputFilename);
        fflush(stdout);
      }
    }
  }





  lprintf ("30 current directory is %s \n", gDirectory->GetName());

  //resultTree->BuildIndex("eventNumber");
  if (parmMismatch) {lprintf("warning: mismatched run parameters found in input files \n");}
  else {lprintf("run parameters match okay across input files\n");}

  lprintf("result tree contains %lli entries, database stride is %i \n", resultTree->GetEntries(), dbStride);
  resultTree->SetBranchAddress("header", &header);
  resultTree->SetBranchAddress("gps", &gpsRaw);
  resultTree->SetBranchAddress("eventSummary", &eventSummary);
  resultTree->SetBranchAddress("circPol", &circPol);
  resultTree->SetBranchAddress("qCuts0", &qCuts[0]);
  resultTree->SetBranchAddress("qCuts1", &qCuts[1]);
  resultTree->SetBranchAddress("maxPwVal", &maxPwVal);
  if (simulationMode) {
    resultTree->SetBranchAddress("pNu", &pnu);
    resultTree->SetBranchAddress("weight", &weight);
    resultTree->SetBranchAddress("truthPhi", &truthPhi);
    resultTree->SetBranchAddress("truthTheta", &truthTheta);
  }
  //turfRateTree->SetBranchAddress("turf", &turfRate);
  //lprintf("building turf rate file index \n");
  //turfRateTree->BuildIndex("payloadTime", "fUniqueID");
  //lprintf("finished building turf rate file index \n");
  
  resultTree->GetEntry(0);
  printf(" test read event %i   weight %f \n", eventSummary->eventNumber, weight);
  int runStartTime = header->triggerTime;
  resultTree->GetEntry(resultTree->GetEntries()-1);
  int runEndTime = header->triggerTime;

  resultTree->GetEntry(2);
  printf(" test read event %i   weight %f \n", eventSummary->eventNumber, weight);

  // instantiate the histograms

  TH1I* hwAngleHist[2] = {0};
  TH1I* corrPeakHist[2] = {0};
  TH1I* hilbPeakHist[2] = {0};
  TH1I* cm3Hist[2] = {0};
  TH1I* peakRatioHist[2] = {0};
  TH1I* anyPolHist[2] = {0};
  TH1I* snrHist0[2] = {0};
  TH1I* snrHist1[2] = {0};
  TH1I* snrHist2[2] = {0};
  TH1I* snrHist3[2] = {0};
  TH1I* timeHist[2] = {0};
  TH1I* timeNsHist[2] = {0};
  TH1I* ptHistPhi[2] = {0};
  TH1I* ptHistTheta[2] = {0};
  TH1F* ptHistPhiSnr[2] = {0};
  TH1F* ptHistThetaSnr[2] = {0};
  TH1F* acceptHistSnr[2] = {0};
  TH1I* linPolHist[2] = {0};
  
  TH2I* skyMapHist[2] = {0};
  TH2I* skyMapPlHist[2] = {0};
  TH2I* skyMapCircHist[2] = {0};
  TH2I* skyMapPlCircHist[2] = {0};
  TH2I* eventHist[2] = {0};
  TH2I* corrHilbHist[2] = {0};
  TH2I* corrSnrHist[2] = {0};
  TH2I* pointingHist[2] = {0};
  TH2I* pointingHistF[2] = {0};
  TH2I* linCircPolHist[2] = {0};
  TH2I* errSnrHist[2] = {0};
  TH2I* sunHist[2] = {0};
  TH2I* sunReflHist[2] = {0}; 
  TGraph* ptTimeHistPhi[2] = {0};
  TH2I* calPulseDistHist[2] = {0};
  TH2F* hpCountsHist[2] = {0};
  TH2F* distSnrAll[2] = {0};
  TH2F* distSnrEff[2] = {0};
  TH2F* corrDistHist[2] = {0};
  TH2F* hilbDistHist[2] = {0};
  TH2F* linPolDistHist[2] = {0};
  TGraph* headingHist = {0};
  
  TH1F* peakSepCircHist = 0;           //  separation (angle) between Lpol and RPol peaks; satellite identification
  TH1F* lesserPeakCircHist = 0;        //  lesser of LPol, RPol peaks
  TH2F* circPeakHist = 0;              //  LPol, RPol correlation peak
  TH2F* cPolPeakDistHist = 0;          //  lesser of RPol, LPol peak  vs  RPol,LPol peak separation
  
  TH2F* totSnrPFFrac[2] = {0};         //  total event count by filter power fraction, SNR
  TH2F* effSnrPFFrac[2] = {0};         //  efficiency by filter power fraction, SNR
  TH2F* distPFFrac_0[2] = {0}; 
  TH2F* distPFFrac_1[2] = {0}; 

  TH2F* linPolPlPhiHist[2] = {0};
  TH2F* linPolPlNPhiHist[2] = {0};

  TH2F* corrPeakPolHist = 0;
  TH1F* corrPeakLesserPolHist = 0;

  TH2F* lesserCircPeakSNRHist[2] = {0};
  TH2F* circPeakSepSNRHist[2] = {0};
  TH2F* lesserCircPeakLDHist[2] = {0};
  TH2F* circPeakSepLDHist[2] = {0};
  TH1F* minCircValHist[2] = {0};


  TH1F* stokesQHist = 0;
  TH1F* stokesUHist = 0;
  TH1F* stokesVHist = 0;

  TH1F* peakTimeHist[2] = {0};

  //stokesIHist = new TH1F("stokesI", "stokesI", 25, 0, 1);
  stokesQHist = new TH1F("stokesQ", "stokes:Q", 50, -1, 1);
  stokesUHist = new TH1F("stokesU", "stokes:U", 50, -1, 1);
  stokesVHist = new TH1F("stokesV", "stokes:V", 50, -1, 1);
    
    // instantiate the histograms
  //bool firstPol = true;
  char title[128];
  char name[128];
  int histNum = 0;
  for (int p=0; p<2; ++p) {
    sprintf(title, "Event reconstructions (pl-north): run %i-%i, %cPol %s", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "skyMapHist%c", pChars[p]);
    vbprintf("hist number  %i %s \n", histNum++, name);
    skyMapHist[p] = new TH2I(name, title, 180, 0, 360, 100, -60, 40);    

    sprintf(title, "Event reconstructions (pl): run %i-%i, %cPol %s", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "skyMapPlHist%c", pChars[p]);
    vbprintf("hist number  %i %s \n", histNum++, name);
    skyMapPlHist[p] = new TH2I(name, title, 180, 0, 360, 100, -60, 40);    

    sprintf(title, "Event reconstructions (pl-north): run %i-%i, %cPol %s", 
            startRunNumber, endRunNumber, pCharsC[p], runDesc);
    sprintf(name, "skyMapHist%c", pCharsC[p]);
    vbprintf("hist number  %i %s \n", histNum++, name);
    skyMapCircHist[p] = new TH2I(name, title, 180, 0, 360, 100, -60, 40);    

    sprintf(title, "Event reconstructions (pl): run %i-%i, %cPol %s", 
            startRunNumber, endRunNumber, pCharsC[p], runDesc);
    sprintf(name, "skyMapPlCircHist%c", pCharsC[p]);
    vbprintf("hist number  %i %s \n", histNum++, name);
    skyMapPlCircHist[p] = new TH2I(name, title, 180, 0, 360, 100, -60, 40);    
    
    sprintf(title, "Peak-to-trigger angle: run %i-%i, %cPol %s", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "hwAngleHist%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    hwAngleHist[p] = new TH1I(name, title, 180, 0, 180);

    sprintf(title, "Correlation peak value: run %i-%i, %cPol %s", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "corrPeakHist%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    corrPeakHist[p] = new TH1I(name, title, 100, 0, 0.6);  

    sprintf(title, "Hilbert peak value: run %i-%i, %cPol %s", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "hilbPeakHist%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    hilbPeakHist[p] = new TH1I(name, title, 100, 0, 200);  

    sprintf(title, "Event localizations, run %i-%i, %cPol %s", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "eventHist%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    eventHist[p] = new TH2I(name, title, 210, -3000, 3000, 210, -3000, 3000);

    sprintf(title, "Correlation peak vs Hilbert peak: run %i-%i, %cPol %s;hilb;corr", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "corrHilbHist%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    corrHilbHist[p] = new TH2I(name, title, 100, 0, 150, 100, 0, 0.4);  

    sprintf(title, "Correlation peak vs SNR: run %i-%i, %cPol %s;SNR;corr", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "corrSnrHist%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    corrSnrHist[p] = new TH2I(name, title, 100, 0, 15, 100, 0, 0.4);  

    sprintf(title, "3rd central moment: run %i-%i, %cPol %s", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "cm3Hist%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    cm3Hist[p] = new TH1I(name, title, 100, -5000, 5000);  

    sprintf(title, "2nd/1st peak ratio: run %i-%i, %cPol %s", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "peakRatioHist%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    peakRatioHist[p] = new TH1I(name, title, 120, 0, 1.2);  

    sprintf(title, "Calibration pulse pointing: run %i-%i, %cPol %s;Az;El", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "pointingHist%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    pointingHist[p] = new TH2I(name, title, 180, -180, 180, 100, -60, 40);  
    sprintf(name, "pointingHistF%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    pointingHistF[p] = new TH2I(name, title, 180, -18, 18, 100, -5, 5);  
    
    sprintf(title, "Total polarization: run %i-%i, %cPol %s", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "anyPolHist%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    anyPolHist[p] = new TH1I(name, title, 100, 0, 1);

    sprintf(title, "Coherent reconstruction SNR (method 0): run %i-%i, %cPol %s", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "snrHist0%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    snrHist0[p] = new TH1I(name, title, 175, 0, 35);

    sprintf(title, "Coherent reconstruction SNR (method 1): run %i-%i, %cPol %s", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "snrHist1%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    snrHist1[p] = new TH1I(name, title, 175, 0, 35);

    sprintf(title, "Coherent reconstruction SNR (method 2): run %i-%i, %cPol %s", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "snrHist2%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    snrHist2[p] = new TH1I(name, title, 175, 0, 35);

    sprintf(title, "Coherent reconstruction SNR (method 3): run %i-%i, %cPol %s", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "snrHist3%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    snrHist3[p] = new TH1I(name, title, 175, 0, 35);

    sprintf(title, "Time-resolved histogram: run %i-%i, %cPol %s", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "timeHist%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    timeHist[p] = new TH1I(name, title, 100, 0, runEndTime-runStartTime);

    sprintf(title, "Payload distance from calPulse: run %i-%i, %cPol %s", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "calPulseDistHist%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    calPulseDistHist[p] = new TH2I(name, title, 100, 0, runEndTime-runStartTime, 60, 300, 600);

    sprintf(title, "Time-resolved pointing: run %i-%i, %cPol %s", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "ptTimeHistPhi%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    ptTimeHistPhi[p] = new TGraph(0);

    sprintf(title, "ns trigger time: run %i-%i, %cPol %s", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "timeNsHist%1i", p);
    if (useTrigNsWindow) {
    vbprintf("hist number  %i %s \n", histNum++, name);
      timeNsHist[p] = new TH1I(name, title, 150, -2000, 1000);
    } else {
    vbprintf("hist number  %i %s \n", histNum++, name);
      timeNsHist[p] = new TH1I(name, title, 100, 0, 1E9);
    }

    sprintf(title, "Linear/Circular polarization: run %i-%i, %cPol %s;lin;circ", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "linCircPolHist%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    linCircPolHist[p] = new TH2I(name, title, 100, 0, 1, 100, 0, 1);

    sprintf(title, "Linear polarization fraction: run %i-%i, %cPol %s;lin;circ", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "linPolHist%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    linPolHist[p] = new TH1I(name, title, 100, 0, 1);

    if (!headingHist) {
      sprintf(title, "Payload Heading vs Time: run %i-%i, %cPol %s;time;heading", 
              startRunNumber, endRunNumber, pChars[p], runDesc);
      sprintf(name, "headingHist%1i", p);
      //headingHist = new TH2F(name, title, 200, 0, runEndTime-runStartTime, 180, 0, 360);
      vbprintf("hist number  %i %s \n", histNum++, name);
      headingHist = new TGraph(0);
    }
    sprintf(title, "Pointing error vs SNR: run %i-%i, %cPol %s;snr;err", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "errSnrHist%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    errSnrHist[p] = new TH2I(name, title, 100, 0, 15, 100, 0, 2);

    sprintf(title, "Localizations, sun coordinates: run %i-%i, %cPol %s;#phi;#theta", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "sunHist%1i", p);
    //sunHist[p] = new TH2I(sunHistNm, sunHistTl, 100, -18, 18, 100, -5, 5);
    vbprintf("hist number  %i %s \n", histNum++, name);
    sunHist[p] = new TH2I(name, title, 360, -180, 180, 100, -60, 40);

    sprintf(title, "Localizations, sun refl coordinates: run %i-%i, %cPol %s;#phi;#theta", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "sunReflHist%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    sunReflHist[p] = new TH2I(name, title, 100, -18, 18, 100, -5, 5);

    sprintf(title, "Healpix bin counts: run %i-%i, %cPol %s;ea;no", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "hpCountsList%1i", p);
    //hpCountsHist[p] = new TGraph2D();
    //hpCountsHist[p]->SetTitle(hpCountsHistTl);
    vbprintf("hist number  %i %s \n", histNum++, name);
    hpCountsHist[p] = new TH2F(name, title, 600, -3000, 3000, 600, -3000, 3000);

    sprintf(title, "Pointing error - #phi: run %i-%i, %cPol %s;#delta#phi", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "ptHistPhi%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    ptHistPhi[p] = new TH1I(name, title, 100, -18, 18);

    sprintf(title, "Pointing error - #theta: run %i-%i, %cPol %s;#delta#theta", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
     sprintf(name, "ptHistTheta%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    ptHistTheta[p] = new TH1I(name, title, 100, -5, 5);

    sprintf(title, "Pointing error by SNR - #phi: run %i-%i, %cPol %s;snr;#delta#phi", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "ptHistPhiSnr%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    ptHistPhiSnr[p] = new TH1F(name, title, 20, 0, 20);

    sprintf(title, "Pointing error by SNR - #theta: run %i-%i, %cPol %s;snr;#delta#theta", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "ptHistThetaSnr%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    ptHistThetaSnr[p] = new TH1F(name, title, 20, 0, 20);

    sprintf(title, "Efficiency by SNR: run %i-%i, %cPol %s;snr;fraction passing", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "acceptHistSnr%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    acceptHistSnr[p] = new TH1F(name, title, 14, -0.5, 13.5);
    
    sprintf(title, "Efficiency by SNR, Cal Pulser distance: run %i-%i, %cPol %s;snr; dist", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "distSnrEff%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    distSnrEff[p] = new TH2F(name, title, 14, -0.5, 13.5, 60, 300, 600);
    sprintf(name, "distSnrAll%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    distSnrAll[p] = new TH2F(name, "", 14, -0.5, 13.5, 60, 300, 600);

    sprintf(title, "Correlation peak vs Cal Pulser distance: run %i-%i, %cPol %s;dist;corr peak", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "corrDistHist%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    corrDistHist[p] = new TH2F(name, title, 60, 300, 600, 60, 0, 0.5);
  
    sprintf(title, "Hilbert peak vs Cal Pulser distance: run %i-%i, %cPol %s;dist;hilb peak", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "hilbDistHist%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    hilbDistHist[p] = new TH2F(name, title, 60, 300, 600, 60, 0, 150);
  
    sprintf(title, "Pol fraction vs Cal Pulser distance: run %i-%i, %cPol %s; Cal Pulser dist;lin pol frac", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "linPolDistHist%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    linPolDistHist[p] = new TH2F(name, title, 60, 300, 600, 50, 0, 1);
    
    sprintf(name, "totSnrPFFrac%1i", p);      
    sprintf(title, "Total Event Count by Filter power fraction, SNR: run %i-%i, %cPol %s;SNR; power frac",
          startRunNumber, endRunNumber, pChars[p], runDesc);
    vbprintf("hist number  %i %s \n", histNum++, name);
    totSnrPFFrac[p] = new TH2F(name, title, 15, 0.5, 15.5, 12, -0.1, 1.1);
    
    sprintf(title, "Efficiency vs Filtered power fraction, SNR: run %i-%i, %cPol %s;SNR; (filtered/raw) power; efficiency",
          startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "effSnrPFFrac%1i", p);      
    vbprintf("hist number  %i %s \n", histNum++, name);
    effSnrPFFrac[p] = new TH2F(name, title, 15, 0.5, 15.5, 12, -0.1, 1.1);
    
    sprintf(title, "(Filtered power / raw power) vs Cal Pulser distance (method 0): run %i-%i, %cPol %s;dist; (filtered/raw) power", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "distPFFrac_0_%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    distPFFrac_0[p] = new TH2F(name, title, 60, 300, 600, 50, -0.1, 1.1);
    
    sprintf(title, "(Filtered power / raw power) vs Cal Pulser distance (method 1): run %i-%i, %cPol %s; dist; (filtered/raw) power", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "distPFFrac_1_%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    distPFFrac_1[p] = new TH2F(name, title, 60, 300, 600, 50, -0.1, 1.1);

    sprintf(title, "LinPol frac vs. PL Az: run %i-%i, %cPol %s; Az; LinPol", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "linPolPlPhiHist%1i", p);    
    vbprintf("hist number  %i %s \n", histNum++, name);
    linPolPlPhiHist[p] = new TH2F(name, title, 180, 0, 360, 55, 0, 1.1);

    sprintf(title, "LinPol frac vs. PL-north Az: run %i-%i, %cPol %s; Az; LinPol", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "linPolPlNPhiHist%1i", p);    
    vbprintf("hist number  %i %s \n", histNum++, name);
    linPolPlNPhiHist[p] = new TH2F(name, title, 180, 0, 360, 55, 0, 1.1);
    
    sprintf(title, "Lesser circPol peak vs SNR: run %i-%i, %cPol %s; SNR; lesserPeak", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "lesserCircPeakSNRHist%1i", p);    
    vbprintf("hist number  %i %s \n", histNum++, name);
    lesserCircPeakSNRHist[p] = new TH2F(name, title, 100, 0, 20, 100, 0, 0.6);
    
    sprintf(title, "circPol peak separation vs SNR: run %i-%i, %cPol %s; SNR; peakSep", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "circPeakSepSNRHist%1i", p);    
    vbprintf("hist number  %i %s \n", histNum++, name);
    circPeakSepSNRHist[p] = new TH2F(name, title, 100, 0, 20, 180, 0, 180);
        
    sprintf(title, "Lesser circPol peak vs SNR_hilb: run %i-%i, %cPol %s; SNR_hilb; lesserPeak", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "lesserCircPeakLDHist%1i", p);    
    vbprintf("hist number  %i %s \n", histNum++, name);
    lesserCircPeakLDHist[p] = new TH2F(name, title, 90, 0, 30, 100, 0, 0.6);
    
    sprintf(title, "circPol peak separation vs SNR_hilb: run %i-%i, %cPol %s; SNR_hilb; peakSep", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "circPeakSepLDHist%1i", p);    
    vbprintf("hist number  %i %s \n", histNum++, name);
    circPeakSepLDHist[p] = new TH2F(name, title, 90, 0, 30, 180, 0, 180);    
    
    sprintf(title, "peak time: run %i-%i, %cPol %s; peak time (ns)", 
            startRunNumber, endRunNumber, pChars[p], runDesc);
    sprintf(name, "peakTimeHist%1i", p);    
    vbprintf("hist number  %i %s \n", histNum++, name);
    peakTimeHist[p] = new TH1F(name, title, 100, 0, 100);    
    
    sprintf(title, "Lesser CPol correlation at %cPol peak location : run %i-%i, %s; lesser cpol value ", 
            pChars[p], startRunNumber, endRunNumber, runDesc);
    sprintf(name, "minCircValHist%1i", p);
    vbprintf("hist number  %i %s \n", histNum++, name);
    minCircValHist[p] = new TH1F(name, title, 100, 0, 2.0);

    //firstPol = false;
  }
  
  sprintf(title, "Correlation peak by pol: run %i-%i, %s; %cpol; %cpol", 
          startRunNumber, endRunNumber, runDesc, pChars[0], pChars[1]);
  sprintf(name, "corrPeakPolHist");
  vbprintf("hist number  %i %s \n", histNum++, name);
  corrPeakPolHist = new TH2F(name, title, 100, 0, 0.6, 100, 0, 0.6);

  sprintf(title, "Lesser correlation peak by pol: run %i-%i, %s; %cpol; %cpol", 
          startRunNumber, endRunNumber, runDesc, pChars[0], pChars[1]);
  sprintf(name, "corrPeakLesserPolHist");
  vbprintf("hist number  %i %s \n", histNum++, name);
  corrPeakLesserPolHist = new TH1F(name, title, 100, 0, 0.6);

  sprintf(title, "CPol peak separation: run %i-%i, %s", 
          startRunNumber, endRunNumber, runDesc);
  sprintf(name, "peakSepCircHist");
  vbprintf("hist number  %i %s \n", histNum++, name);
  peakSepCircHist = new TH1F(name, title, 180, 0, 180);    

  sprintf(title, "Lesser CPol peak : run %i-%i, %s", 
          startRunNumber, endRunNumber, runDesc);
  sprintf(name, "lesserPeakCircHist");
  vbprintf("hist number  %i %s \n", histNum++, name);
  lesserPeakCircHist = new TH1F(name, title, 100, 0, 0.6);

  sprintf(title, "CPol peak value: run %i-%i, %s; %cpol; %cpol", 
          startRunNumber, endRunNumber, runDesc, pCharsC[0], pCharsC[1]);
  sprintf(name, "circPeakHist");
  vbprintf("hist number  %i %s \n", histNum++, name);
  circPeakHist = new TH2F(name, title, 100, 0, 0.3, 100, 0, 0.3);

  sprintf(title, "Lesser CPol peak value: run %i-%i, %s; lesser cpol corr peak; cpol peak separation ", 
          startRunNumber, endRunNumber, runDesc);
  sprintf(name, "cPolPeakDistHist");
  vbprintf("hist number  %i %s \n", histNum++, name);
  cPolPeakDistHist = new TH2F(name, title, 100, 0, 0.6, 180, 0, 180);

  //lprintf (" 40 current directory is %s \n", gDirectory->GetName());

  // So we want to set it up so some of the quality cuts (and other later cuts) happen for 
  //  both pols at the same time.  to do this we are going to LEAVE THE CURRENT Q CUTS ALONE! 
  //  and add in a new set of these variables for the polarization independent cuts.  which we will call 
  //  the preQCuts.  these cuts are:
  //   1) RF Trigger Cut (trigger type cut)
  //   2) Payload Blast cut **Failes if EITHER POL Fails** 
  //   3) Calibration pulser cut (from a1)
  //   4) Satellite stripe cut (from a2)
  //   5) cpol seperation cut (from a2/af)
  //   6) Waveform Saturation Cut (from qcuts)
  //   7) DC Offset Cut (from qcuts)
  //   8) Incomplete waveform cut (from qcuts)
  //   /) cpol strength cut (from a2/af)  (Weirdly, this does depend on linPols, so it doesnt happen here!
  //   6) L3 Trigger Cut (as an OR!)
  //
  // Defs for these have been added to analysisCuts.h
  //
  // AFTER these cuts, the starting point for the qCuts will be:
  //  (in H) the events that pass these cuts AND the L3 Trigger cut (in H)
  //  (in V) the events that pass these cuts AND the L3 Trigger cut (in V)



  /*
  int numEvents = 0;
  int qualityCutCounts[2][NUM_Q_CUTS] = {{0}};
  int qualityCutUniqueCounts[2][NUM_Q_CUTS] = {{0}};
  int qualityCutOrderedCounts[2][NUM_Q_CUTS] = {{0}};
  int analysisCutCounts[2][NUM_A_CUTS] = {{0}};
  int analysisCutUniqueCounts[2][NUM_A_CUTS] = {{0}};
  int analysisCutOrderedCounts[2][NUM_A_CUTS] = {{0}};
  int eventCounts[2] = {0};
  int qEventCounts[2] = {0};
  int aEventCounts[2] = {0};
  int qualityTotalOrderedCount[2] = {0};
  int analysisTotalOrderedCount[2] = {0};
  int peakFailCounts[2] = {0};

  int qualityCutCounts5[2][NUM_Q_CUTS] = {{0}};
  int qualityCutUniqueCounts5[2][NUM_Q_CUTS] = {{0}};
  int qualityCutOrderedCounts5[2][NUM_Q_CUTS] = {{0}};
  int analysisCutCounts5[2][NUM_A_CUTS] = {{0}};
  int analysisCutUniqueCounts5[2][NUM_A_CUTS] = {{0}};
  int analysisCutOrderedCounts5[2][NUM_A_CUTS] = {{0}};
  int eventCounts5[2] = {0};
  int qEventCounts5[2] = {0};
  int qualityTotalOrderedCount5[2] = {0};
  int analysisTotalOrderedCount5[2] = {0};
  */

  int numEvents = 0;
  float qualityCutCounts[2][NUM_Q_CUTS] = {{0}};
  float qualityCutUniqueCounts[2][NUM_Q_CUTS] = {{0}};
  float qualityCutOrderedCounts[2][NUM_Q_CUTS] = {{0}};
  float analysisCutCounts[2][NUM_A_CUTS] = {{0}};
  float analysisCutUniqueCounts[2][NUM_A_CUTS] = {{0}};
  float analysisCutOrderedCounts[2][NUM_A_CUTS] = {{0}};
  float eventCounts[2] = {0};
  float qEventCounts[2] = {0};
  float aEventCounts[2] = {0};
  float qualityTotalOrderedCount[2] = {0};
  float analysisTotalOrderedCount[2] = {0};
  int peakFailCounts[2] = {0};
  
  float qualityCutCounts5[2][NUM_Q_CUTS] = {{0}};
  float qualityCutUniqueCounts5[2][NUM_Q_CUTS] = {{0}};
  float qualityCutOrderedCounts5[2][NUM_Q_CUTS] = {{0}};
  float analysisCutCounts5[2][NUM_A_CUTS] = {{0}};
  float analysisCutUniqueCounts5[2][NUM_A_CUTS] = {{0}};
  float analysisCutOrderedCounts5[2][NUM_A_CUTS] = {{0}};
  float eventCounts5[2] = {0};
  float qEventCounts5[2] = {0};
  float qualityTotalOrderedCount5[2] = {0};
  float analysisTotalOrderedCount5[2] = {0};
  
  // defs for the preQCuts 
  
  float preQualityCutCounts[2][NUM_PRE_Q_CUTS] = {{0}};
  float preQualityCutUniqueCounts[2][NUM_PRE_Q_CUTS] = {{0}};
  float preQualityCutOrderedCounts[2][NUM_PRE_Q_CUTS] = {{0}};
  float preQualityTotalOrderedCount[2] = {0};
  
  float preQualityCutCounts5[2][NUM_PRE_Q_CUTS] = {{0}};
  float preQualityCutUniqueCounts5[2][NUM_PRE_Q_CUTS] = {{0}};
  float preQualityCutOrderedCounts5[2][NUM_PRE_Q_CUTS] = {{0}};
  float preQualityTotalOrderedCount5[2] = {0};
  
  float preQEventCounts[2] = {0};
  float preQEventCounts5[2] = {0};
  
  vector<int> analysisCutEventList[2][NUM_A_CUTS] = {{vector<int>(0)}}; 
  vector<int> passingEvents[2] = {vector<int>(0)}; 
  vector<int> passingEventRuns[2] = {vector<int>(0)}; 
  //for (int p=0; p<2; ++p) {
    //passingEvents[p] = vector<int>(0);
    //passingEventRuns[p] = vector<int>(0);
    //for (int k=0; k<NUM_A_CUTS; ++k) analysisCutEventList[p][k] = vector<int>(0);
  //}
  vector<int> misRecs[2] = {vector<int>(0)};

  int redundancyCounts[2][NUM_A_CUTS][NUM_A_CUTS] = {{{0}}};
  int peakViolationCount = 0;

  int startEntry = 0;
  if (singleEventNumber > 0 || useEventList) {
    resultTree->BuildIndex("eventNumber", "circPol");    
  }
  if (singleEventNumber > 0) {
    startEntry = resultTree->GetEntryNumberWithIndex(singleEventNumber, 0);
    if (startEntry < 0) {
     lprintf("error: invalid single event entry number %i, terminating \n", startEntry);
      return 1;
    } else {printf("single event entry number is %i \n", startEntry);}
  }
  // these are delcared here for post-loop processing for single event run
  double plEa=0, plNo=0;
  double srcEa[2]={0}, srcNo[2]={0};
  double plSrcPhi[2] = {0};
  vector<float> errEllipseParms0[2];  for (int p=0; p<2; ++p) {errEllipseParms0[p] = {0, 0, 0};}
  vector<float> errEllipseParms1[2];   for (int p=0; p<2; ++p) {errEllipseParms1[p] = {0, 0, 0, 0};}
  vector<vector<double> > errHexVertices[2];  for (int p=0; p<2; ++p) {errHexVertices[p] = vector<vector<double> >(0);}
  vector<vector<double> > errHexPixels[2];  for (int p=0; p<2; ++p) {errHexPixels[p] = vector<vector<double> >(0);}
  vector<vector<double> > errHexRayTrace[2] = {vector<vector<double> >(0)};
  vector<vector<double> > plMapPoints(0);
  float calPulseDist[2] = {0};
  // loop through the interferometry results tree and process events
  int runNumber = -1;
  //bool plotPl = false;
  
  map<int, float> hpEventCounts[2];
  int hpEvents[2] = {0};
  TGraph* plPathGr = new TGraph(0);
  plPathGr->SetFillStyle(0);
  plPathGr->SetLineColor(kRed);
  // TODO support option to read event list from a file.
  vector<int> inputEventList(0);
  int thisEvent;
  char thisPol;
  int thisRun;
  int maxE = maxEvents;
  if (maxE == 0) maxE = resultTree->GetEntries();
  if (useEventList) {
    char eventListFilepath[1024];
    sprintf(eventListFilepath, "%s/%s.txt", inputDirName, eventListFilename);
    lprintf(" input event list filepath is %s \n", eventListFilepath);
    ifstream eventListFile(eventListFilepath);
    bool done = false;
    while (eventListFile.good() && !done) {
      eventListFile >> thisEvent >> thisPol >> thisRun;
      inputEventList.push_back(thisEvent);  
      if (thisPol == 'V') done = true;
    }
    maxE = inputEventList.size();
    lprintf("event list size is %li \n", inputEventList.size());  
  } else if (singleEventNumber > 0) {
    maxE = startEntry+1;
  }
  //lprintf ("50 current directory is %s \n", gDirectory->GetName());  
  lprintf("looping through events: maximum is %i \n", maxE);
  int swapCount = 0;
  fflush(stdout);  
  // loop through events    ---------------------------------------------------------------------------------
  for (int e=startEntry; e<maxE; e += 1 + 2*(dbStride-1)) {
    vbprintf(" entry number = %i \n", e);
    //vbprintf ("50 current directory is %s \n", gDirectory->GetName());    
    resultTree->SetBranchAddress("header", &header);
    resultTree->SetBranchAddress("eventSummary", &eventSummary);
    if (simulationMode) {
      resultTree->SetBranchAddress("weight", &weight);
    }
  
    if (useEventList) {     // TODO either do this better of deprecate it
      resultTree->GetEntryWithIndex(inputEventList[e], 0);
    } else {
      vbprintf(" getting entry number = %i \n", e);
      resultTree->GetEntry(e);
      vbprintf(" entry number = %i gotten \n", e);
      vbprintf(" event weight = %f \n", weight);
    }
    vbprintf("event %i run %i retrieved, entry number %i, circpol flag is %i \n", header->eventNumber, eventSummary->run, e, circPol);
    vbprintf("  trigger time is %i:%09i \n", header->triggerTime, header->triggerTimeNs);
    UsefulAdu5Pat* gps = new UsefulAdu5Pat(gpsRaw);
    //cout << gps->pitch << " " << gps->roll << endl;
    if (eventSummary->run != runNumber) {
      lprintf("run number %i \n", eventSummary->run); fflush(stdout);
      runNumber = eventSummary->run;
      plMapPoints.push_back({gps->latitude, gps->longitude});
      double ea1, no1;
      bedmap->LonLattoEaNo(gps->longitude, gps->latitude, ea1, no1);
      plPathGr->SetPoint(plPathGr->GetN(), ea1/1000, no1/1000);
    }
    int numACuts[2] = {0};
    int numQCuts[2] = {0};
    int numACuts5[2] = {0};
    int numQCuts5[2] = {0};
    int numPreQCuts[2] = {0};
    int numPreQCuts5[2] = {0};

    int aCuts[2][NUM_A_CUTS] = {{0}};
    float peakRatio[2] = {0};
    float linPolFrac[2] = {0};
    float deadTimeFrac[2] = {0};
    float circPolFrac[2] = {0};
    float anyPolFrac[2] = {0};
    float filterPowerFrac_0[2] = {0};
    float filterPowerFrac_1[2] = {0};
    float ldCutVal[2] = {0};
    float minCircPeakVal[2] = {0};
    
    vector<int> maxPeakIndex[2] = {vector<int>(0)};
    bool keep[2] = {true};                          //this marks events that pass/fail the preQCuts OR the qCuts
    bool keepPreQ[2] = {true};                      //this marks spicifically events that pass/fail the preQCuts
    float eventWeight[2] = {0};
    int maxWfAnt[2] = {-1};
    float maxWfSnr[2] = {0};
    //if ((circPol && doCircPol==1) || (!circPol && doCircPol==0)) {
    //vbprintf ("60 current directory is %s \n", gDirectory->GetName());  
    
    if (circPol==0) {
      AnitaEventSummary::SourceHypothesis calSource = eventSummary->wais;      
      int netTrigTime = header->triggerTime - runStartTime;
      if (strcmp(calPulser, "WAIS") == 0) {
        calSource = eventSummary->wais;
      } else if (strcmp(calPulser, "LDB") == 0) {
        calSource = eventSummary->ldb;
      }
      //vbprintf ("60.5 current directory is %s \n", gDirectory->GetName());          
      int sourceDelay = (int)(gps->getTriggerTimeNsFromSource(ksLat, ksLon, ksAlt));
      //int nsDiff = (int)((float)calSource.distance/0.3 - (float)header->triggerTimeNs + pulserNsTime);
      int nsDiff = sourceDelay - header->triggerTimeNs + pulserNsTime;
      vbprintf("   cal pulser offset is %i, time delay is %i, differential is %i \n", sourceDelay, pulserNsTime, nsDiff);
      bool calPulseEvent = (nsDiff>-1000 && nsDiff<1500 && runNumber>320 && runNumber<365);
      calPulseEvent |= (nsDiff>-2000 && nsDiff<1500 && runNumber>140 && runNumber<155);    // TODO check these numbers
      if (useTrigNsWindow && !calPulseEvent) continue; // intentional structure violation
      //vbprintf ("60.7 current directory is %s \n", gDirectory->GetName());          

      vbprintf("\nevent %i ------------------------------------------\n", header->eventNumber);
      for (int p=0; p<2; ++p) {
        vbprintf("    qCuts[%i]: ", p); 
        for (int k=0; k<NUM_Q_CUTS; ++k) {vbprintf("%i ", qCuts[p][k]);}
        vbprintf("\n");
      }
      //vbprintf ("61 current directory is %s \n", gDirectory->GetName());          
      
      // loop through polarizations and process the event ----------------------------------------------------
      //bool firstPol = true;
      for (int p=0; p<2; ++p) {
        //vbprintf ("62 current directory is %s \n", gDirectory->GetName());          
      
        keep[p] = true;
        keepPreQ[p] = true;
        //vbprintf(" top-of-loop keep-it value[%i] is %i \n", p, keep[p]);
        filterPowerFrac_0[p] = 0;
        filterPowerFrac_1[p] = 0;
        //  TODO calculate power loss from the individual WaveformInfo's        
        uint16_t l3Pattern = circPol ? header->l3TrigPattern || header->l3TrigPatternH : 
                p==0 ? header->l3TrigPatternH : header->l3TrigPattern;
        float sa = 0;
        float ca = 0;
        vbprintf("--------------- %cpol l3 trigger pattern is %1x:  ", pChars[p], l3Pattern);
        for (int s=0; s<16; ++s) {
          if (l3Pattern & (1 << s)) {
            sa += sin(s*M_PI/8.0);
            ca += cos(s*M_PI/8.0);
            vbprintf(" %2i", s);
          //} else {
            //vbprintf(" 0");
          }
        }
        vbprintf("\n");
        float totalPowerRaw = 0;
        float totalPowerFiltered = 0;
        for (int ant = 0; ant < NUM_SEAVEYS; ++ant) {
          if (header->isInL3Pattern(geom->getPhiFromAnt(ant), (AnitaPol::AnitaPol_t)p)) {
            //printf("event %i  %f  %f \n", header->eventNumber, eventSummary->inputWfRaw[p][ant].totalPower, 
                    //eventSummary->inputWfFiltered[p][ant].totalPower);
            totalPowerRaw += eventSummary->inputWfRaw[p][ant].totalPower;
            totalPowerFiltered += eventSummary->inputWfFiltered[p][ant].totalPower;
          }
        }
        if (totalPowerRaw > 0) {filterPowerFrac_0[p] = totalPowerFiltered/totalPowerRaw;}
        if (eventSummary->flags.medianPower[0] > 0) {filterPowerFrac_1[p] = eventSummary->flags.medianPowerFiltered[0]/eventSummary->flags.medianPower[0];}
        vbprintf("event %i filter_power_ratio_0 (%f / %f) = %f \n", header->eventNumber, totalPowerFiltered, totalPowerRaw, filterPowerFrac_0[p]);
        vbprintf("event %i filter_power_ratio_1 (%f / %f) = %f \n", header->eventNumber, eventSummary->flags.medianPowerFiltered[0],
                eventSummary->flags.medianPower[0], filterPowerFrac_1[p]);
        float triggerAngle = atan2(sa, ca) * 180.0 / M_PI;
        while (triggerAngle < 0) triggerAngle += 360;
        vbprintf("   trigger angle is %8.6f \n", triggerAngle);

        //maxPeakIndex[p] = vector<int>(0);
        for (int k=0; k<eventSummary->nPeaks[p]; ++k) {maxPeakIndex[p].push_back(k);}
        //vbprintf ("64 current directory is %s \n", gDirectory->GetName());          

        // traverse the array for peak order violations and sort if flag is true
        bool sortPeaks = false;
        for (int j=0; j<eventSummary->nPeaks[p]; ++j) {
          for (int k=0; k<eventSummary->nPeaks[p]-j-1; ++k) {
            if (eventSummary->peak[p][maxPeakIndex[p][k]].value < eventSummary->peak[p][maxPeakIndex[p][k+1]].value) {
              if (pChars[p]=='H' ) {
                //printf("  event %i %cpol comparing peaks %i (%f, %f, %f) and %i (%f, %f, %f)   ratio is %f  \n", 
                //         eventSummary->eventNumber, pChars[p], k, 
                //          eventSummary->peak[p][maxPeakIndex[p][k]].phi, eventSummary->peak[p][maxPeakIndex[p][k]].theta, 
                //          eventSummary->peak[p][maxPeakIndex[p][k]].value, 
                //          k+1, eventSummary->peak[p][maxPeakIndex[p][k+1]].phi, eventSummary->peak[p][maxPeakIndex[p][k+1]].theta, 
                //          eventSummary->peak[p][maxPeakIndex[p][k+1]].value,
                //          eventSummary->peak[p][maxPeakIndex[p][k]].value / eventSummary->peak[p][maxPeakIndex[p][k+1]].value);
                //printf("    ----- (not really) swapping peaks \n");
                ++swapCount;
                if (sortPeaks) {  
                  int temp = maxPeakIndex[p][k];
                  maxPeakIndex[p][j] = maxPeakIndex[p][k+1];
                  maxPeakIndex[p][k] = maxPeakIndex[p][k+1];
                  maxPeakIndex[p][k+1] = temp;
                }
              }
            }
          }
        }
        
        // calculate stuff needed before applying cuts  
        if (eventSummary->nPeaks[p] < 2) {
          lprintf("warning: <2 interferometry peaks found for event %i %cpol, (n=%i) skipping \n", eventSummary->eventNumber, pChars[p], eventSummary->nPeaks[p]);
          keep[p] = false;
          ++peakFailCounts[p];
          //return 3;
          continue; 
        }  
        //vbprintf ("66 current directory is %s \n", gDirectory->GetName());          

        //AnitaEventSummary::PointingHypothesis mainPeak = eventSummary->peak[p][maxPeakIndex[p][0]];
        AnitaEventSummary::PointingHypothesis mainPeak = eventSummary->peak[p][0];
        eventWeight[p] = (simulationMode) ? weight : 1.0;
        vbprintf(" event weight = %f \n", eventWeight[p]);
        // TODO get index of waveform with maximum peak-to-peak difference
        float maxPeakVal = 0;
        for (int ant=0; ant<NUM_SEAVEYS_ANITA3; ++ant) {
          if (mainPeak.value > maxPeakVal) {
            maxPeakVal = mainPeak.value;
            maxWfAnt[p] = ant;
          }
        }
        maxWfSnr[p] = eventSummary->inputWfFiltered[p][maxWfAnt[p]].snr2;
        float myHwAngle = triggerAngle - mainPeak.phi;
        while (myHwAngle <= -180) myHwAngle +=360;
        while (myHwAngle > 180) myHwAngle -=360;
        vbprintf(" pointing result is (%6.3f,%5.3f),  cal pulser is at (%6.3f, %5.3f), hwAngle=%6.3f, myHwAngle=%6.3f\n", 
                mainPeak.phi, mainPeak.theta, calSource.phi, 
                calSource.theta, mainPeak.hwAngle, myHwAngle);

        // finagle localization to continent
        double thetaFinag = mainPeak.theta + 0.0;
        //double newLon, newLat, newAlt, newAdj;
        //printf("old lat/lon/alt is %f, %f, %f \n", mainPeak.latitude, mainPeak.longitude, mainPeak.altitude);
        //printf(" phi = %f, theta = %f \n" , mainPeak.theta,mainPeak.phi);
        if (!gps->traceBackToContinent(mainPeak.phi *M_PI/180.0, thetaFinag *M_PI/180.0, &newLon[p], &newLat[p], &newAlt[p], &newAdj[p])) {
          newLon[p] = -9999;
          newLat[p] = -9999;
          newAlt[p] = -9999;
        };
        //printf(" phi = %f, theta = %f, lat = %f, lon = %f, alt = %f \n", mainPeak.theta, mainPeak.phi, newLat[p], newLon[p], newAlt[p]); 
        //printf("new lat/lon/alt is %f, %f, %f \n", newLat[p], newLon[p], newAlt[p]);
        mainPeak.latitude = newLat[p];
        mainPeak.longitude = newLon[p];
        mainPeak.altitude = newAlt[p];
        // easting/northing of payload and source
        //bedmap->SurfaceLonLattoEN(mainPeak.longitude, mainPeak.latitude, ec, nc);
        //bedmap->ENtoEaNo(ec, nc, srcEa[p], srcNo[p]);
        bedmap->LonLattoEaNo(mainPeak.longitude, mainPeak.latitude, srcEa[p], srcNo[p]);
        //bedmap->SurfaceLonLattoEN(gps->longitude, gps->latitude, ec, nc);
        //bedmap->ENtoEaNo(ec, nc, plEa[p], plNo[p]);
        bedmap->LonLattoEaNo(gps->longitude, gps->latitude, plEa, plNo);
        plSrcPhi[p] = atan2(srcNo[p]-plNo, srcEa[p]-plEa);
        float plSrcDist = sqrt((srcNo[p]-plNo)*(srcNo[p]-plNo) + (srcEa[p]-plEa)*(srcEa[p]-plEa));

        // get the event cohSnr2, just used for some extra print out statements
        float cohSnr2  = eventSummary->coherent[p][maxPeakIndex[p][0]].snr2;
        // also... we need the cohSnr2 for both pols, so... get that too. (BUT ONLY IF p == 1!!)

        float cohSnr2H = 0;
        float cohSnr2V = 0;

        if (p == 0) {
          cohSnr2H = eventSummary->coherent[0][maxPeakIndex[0][0]].snr2; 
        } else if (p == 1) {
          cohSnr2H = eventSummary->coherent[0][maxPeakIndex[0][0]].snr2;
          cohSnr2V = eventSummary->coherent[1][maxPeakIndex[1][0]].snr2;
        }


        // error ellipse
        // phi is the angle of the payload-to-source vector
        // TODO if we can't construct the error ellipse, fail the reconstruct-to-continent cut
        vbprintf("  main peak latitude is %f \n", mainPeak.latitude);
        if (mainPeak.latitude > -9999) {
          float effTheta = mainPeak.theta - asin(plSrcDist/EARTH_POLAR_RADIUS)*180.0/M_PI;
          vbprintf("    Pl-to-Source dist = %f, Pl-phi= %f, Pl-theta= %f, effective theta = %f \n", plSrcDist, mainPeak.phi, mainPeak.theta, effTheta);
          errEllipseParms0[p] = surfaceEllipseParms0(effTheta, plSrcDist, gps->altitude, pointingErr[0], pointingErr[1]);
          errEllipseParms1[p] = surfaceEllipseParms1(mainPeak.phi, mainPeak.theta, plSrcDist, pointingErr[0], pointingErr[1], gps, bedmap);
          //vbprintf("     Error ellipse parameters (method 0): centerD=%f, rSemi=%f, thSemi=%f \n", errEllipseParms0[p][0], errEllipseParms0[p][1], errEllipseParms0[p][2]);
          errHexRayTrace[p] = errorPixelsRayTrace(mainPeak.phi, mainPeak.theta, pointingErr[0], pointingErr[1], gps, bedmap);
          if (errEllipseParms1[p][0] > 0) {
            //vbprintf("     Source (Ea,No) = %f, %f", srcEa[p], srcNo[p]);
            //vbprintf("     Payload (Ea,No) = %f, %f", plEa, plNo);
            //printf("     Error ellipse parameters (method 1): rSemi=%f, thSemi=%f, center=%f,%f \n", 
            //        errEllipseParms1[p][0], errEllipseParms1[p][1], errEllipseParms1[p][2], errEllipseParms1[p][3]);
            errHexVertices[p] = hexagonVertices(errEllipseParms1[p][2], errEllipseParms1[p][3], errEllipseParms1[p][0], errEllipseParms1[p][1], plSrcPhi[p]);
            //vbprintf("     Error hexagon vertices :"); for (int k=0; k<6; ++k) {vbprintf("       %f, %f \n", errHexVertices[p][0][k], errHexVertices[p][1][k]);}            
            errHexPixels[p] = ellipsePixels(errEllipseParms1[p][2], errEllipseParms1[p][3], errEllipseParms1[p][0], errEllipseParms1[p][1], plSrcPhi[p]);
            float pdfSum = 0;
            //vbprintf("     Error pixel vertices : \n"); 
            for (int k=0; k<errHexPixels[p][0].size(); ++k) {
              //vbprintf("       %f, %f, pdf=%f \n", errHexPixels[p][0][k], errHexPixels[p][1][k], errHexPixels[p][2][k]);
              pdfSum += errHexPixels[p][2][k];
            }
            vbprintf("    pdf total =         %.2f \n", pdfSum);
          } else {
            //vbprintf("    no error ellipse for event %i \n", header->eventNumber);
            //vbprintf("          main peak latitude is %f \n", mainPeak.latitude);
            errHexPixels[p]=vector<vector<double> >(3);
            for (int k=0; k<3; ++k) {errHexPixels[p][k] = vector<double>(0);}
          }
        }
        // set up for cut processing
        //vbprintf ("70 current directory is %s \n", gDirectory->GetName());          
        vbprintf(" initial keep-it value[%i] is %i \n", p, keep[p]);

        ///////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // go into preQCuts first, then after, go into the Q_Cuts.  These each of these cuts is done here, manually.
        //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////       

        // get a fresh start before doing these cuts.
        for(int k=0; k<preQCutOrder.size(); k++)
        {
          preQCuts[p][k] = 0;
        }

        // RF Trigger Cut
        // -This is just the '0' qcut
        
        // if this cut is on, and event fails it
        if (    ( (preQualityCutMask[p][PRE_Q_CUT_TRIGGER_TYPE] == 1 || qualityCutMask[p][Q_CUT_TRIGGER_TYPE]== 1) && (qCuts[p][Q_CUT_TRIGGER_TYPE]==1))
             || ( (preQualityCutMask[p][PRE_Q_CUT_TRIGGER_TYPE] ==-1 || qualityCutMask[p][Q_CUT_TRIGGER_TYPE]==-1) && (qCuts[p][Q_CUT_TRIGGER_TYPE]==0))
           ) 
        {
          keep[p] = false;
          preQCuts[p][PRE_Q_CUT_TRIGGER_TYPE] = 1;
        }


        // Payload Blast (as an OR) Cut
        // -This is qcut 2, but it fails if it failed in either pol

        // if this cut is on, and the event fails, OR this cut is accepting only events that fail the cut?
        if (    (     (preQualityCutMask[p][PRE_Q_CUT_PAYLOAD_BLAST]== 1 || qualityCutMask[p][Q_CUT_PAYLOAD_BLAST]== 1) 
                   && ((qCuts[0][Q_CUT_PAYLOAD_BLAST]==1) || (qCuts[1][Q_CUT_PAYLOAD_BLAST]==1)) )
             || (     (preQualityCutMask[p][PRE_Q_CUT_PAYLOAD_BLAST]==-1 || qualityCutMask[p][Q_CUT_PAYLOAD_BLAST]==-1) 
                   && ((qCuts[0][Q_CUT_PAYLOAD_BLAST]==0) && (qCuts[1][Q_CUT_PAYLOAD_BLAST]==0)) )
           )
        {
          keep[p] = false;
          preQCuts[p][PRE_Q_CUT_PAYLOAD_BLAST] = 1;
        }


        // PRE_Q_CUT_WF_SATURATION, as an OR

        if (    (    (preQualityCutMask[p][PRE_Q_CUT_WF_SATURATION]== 1 || qualityCutMask[p][Q_CUT_WF_SATURATION]== 1) 
                  && ((qCuts[0][Q_CUT_WF_SATURATION]==1) || (qCuts[1][Q_CUT_WF_SATURATION]==1)) )
             || (    (preQualityCutMask[p][PRE_Q_CUT_WF_SATURATION]==-1 || qualityCutMask[p][Q_CUT_WF_SATURATION]==-1) 
                  && ((qCuts[0][Q_CUT_WF_SATURATION]==0) && (qCuts[1][Q_CUT_WF_SATURATION]==0)) )
           )
        {
          keep[p] = false;
          preQCuts[p][PRE_Q_CUT_WF_SATURATION] = 1;
        }


        // PRE_Q_CUT_WF_DC_OFFSET, as an OR

        if (    (    (preQualityCutMask[p][PRE_Q_CUT_WF_DC_OFFSET]== 1 || qualityCutMask[p][Q_CUT_WF_DC_OFFSET]== 1) 
                  && ((qCuts[0][Q_CUT_WF_DC_OFFSET]==1) || (qCuts[1][Q_CUT_WF_DC_OFFSET]==1)) )
             || (    (preQualityCutMask[p][PRE_Q_CUT_WF_DC_OFFSET]==-1 || qualityCutMask[p][Q_CUT_WF_DC_OFFSET]==-1) 
                  && ((qCuts[0][Q_CUT_WF_DC_OFFSET]==0) && (qCuts[1][Q_CUT_WF_DC_OFFSET]==0)) )
           )
        {
          keep[p] = false;
          preQCuts[p][PRE_Q_CUT_WF_DC_OFFSET] = 1;
        }


        // PRE_Q_CUT_SHORT_WF, as an OR

        if (    (    (preQualityCutMask[p][PRE_Q_CUT_SHORT_WF]== 1 || qualityCutMask[p][Q_CUT_SHORT_WF]== 1) 
                  && ((qCuts[0][Q_CUT_SHORT_WF]==1) || (qCuts[1][Q_CUT_SHORT_WF]==1)) )
             || (    (preQualityCutMask[p][PRE_Q_CUT_SHORT_WF]==-1 || qualityCutMask[p][Q_CUT_SHORT_WF]==-1) 
                  && ((qCuts[0][Q_CUT_SHORT_WF]==0) && (qCuts[1][Q_CUT_SHORT_WF]==0)) )
           )
        {
          keep[p] = false;
          preQCuts[p][PRE_Q_CUT_SHORT_WF] = 1;
        }


        //cout << header->eventNumber << " " << p << " " << qCuts[0][Q_CUT_PAYLOAD_BLAST] << " " << qCuts[1][Q_CUT_PAYLOAD_BLAST] << endl;

        // Cal Pulser Cut
        // -This is done later in this code (and now also here), and is a a1 cut

        bool preWaisPulse = UCorrelator::isWAISHPol(gps, header);
        //vbprintf ("85.086 current directory is %s \n", gDirectory->GetName());
        if (    ( preQualityCutMask[p][PRE_Q_CUT_CAL_PULSER]== 1 &&  preWaisPulse )
             || ( preQualityCutMask[p][PRE_Q_CUT_CAL_PULSER]==-1 && !preWaisPulse )  
           ) 
        {
          keep[p] = false;
          preQCuts[p][PRE_Q_CUT_CAL_PULSER] = 1;
        }


        // Satellite Stripe Cut
        // -Added in from runAnalysisStage02.cxx
        // -Normally all the cPol stuff happens later, so we need to... shove it in here now
        //  I guess, and hopefully it wont break anything.


        // The data is listed as follows:
        //  EventNumber Cpol
        //  1           0
        //  1           1
        //  2           0
        //  2           1
        //  3           0
        //  3           1
        //  ect...

        // get the cpol stuff for the next three cuts, for the next section of 
        // code, the trees will point at cPol data, not linPol data        



        float cPolMinDist = 0;
        float cPolPeak[2] = {0};
        float cPolPhi[2] = {0};
        float cPolTheta[2] = {0};
        int thisEventNum = header->eventNumber;        

        if (e != maxE) {
          resultTree->GetEntry(e+1);

          // check that stuff is okay
          if (circPol!=1 || header->eventNumber != thisEventNum ) {
            lprintf("warning: no circPol entry found for event %i (next event=%i circpol=%i) \n", thisEventNum, header->eventNumber, circPol); 
          } else { 
            int k = p; // dumb short cut to fix a copy paste
            cPolPeak[0] = eventSummary->peak[0][0].value;
            cPolPeak[1] = eventSummary->peak[1][0].value;
            cPolPhi[0] = eventSummary->peak[0][0].phi;
            cPolPhi[1] = eventSummary->peak[1][0].phi;
            cPolTheta[0] = eventSummary->peak[0][0].theta;
            cPolTheta[1] = eventSummary->peak[1][0].theta;
                  
            // find the minimum distance between first L(R) peak and each R(L) peak
            vector<float> peakDists(0);
            AnitaEventSummary::PointingHypothesis peak0 = eventSummary->peak[0][0];
            for (int kr=0; kr<cPolNumPeaks; ++kr) {
              AnitaEventSummary::PointingHypothesis peak1 = eventSummary->peak[1][kr];
              float cPolDist = vectorAngle(peak0.phi, peak1.phi, peak0.theta, peak1.theta);
              peakDists.push_back(cPolDist);
              vbprintf(" kr=%i, lPeak=(%f,%f) rPeak=(%f,%f) cPolDist = %f \n", kr, eventSummary->peak[0][0].phi, eventSummary->peak[0][0].theta,
                         eventSummary->peak[1][kr].phi, eventSummary->peak[1][kr].theta, cPolDist);
            }
            AnitaEventSummary::PointingHypothesis peak1 = eventSummary->peak[1][0];
            for (int kl=1; kl<cPolNumPeaks; ++kl) {
              AnitaEventSummary::PointingHypothesis peak1 = eventSummary->peak[0][kl];
              float cPolDist = vectorAngle(peak0.phi, peak1.phi, peak0.theta, peak1.theta);
              peakDists.push_back(cPolDist);
              vbprintf(" kl=%i, lPeak=(%f,%f) rPeak=(%f,%f) cPolDist = %f \n", kl, eventSummary->peak[0][kl].phi, eventSummary->peak[0][kl].theta,
                         eventSummary->peak[1][0].phi, eventSummary->peak[1][0].theta, cPolDist);
            }
            cPolMinDist = *min_element(peakDists.begin(), peakDists.end());
          }
        }

        // leaving weird little mini CPol loop.  setting resultTree back to normal!
        resultTree->GetEntry(e);   // check that this works?


        float modLPhi = fmod((cPolPhi[0] - gps->heading + 360), 360);
        float mapPos = gps->longitude - modLPhi;
        while (mapPos > 160) { mapPos -= 360; } //skewed limits because of our logic checks.
        while (mapPos < -200) { mapPos += 360; }
        float plLat = gps->latitude;
        float plLon = gps->longitude;
        float plAlt = (gps->altitude); //in m?
        float sAlt = 35786000; //in m
        float sLat = -10.0;
        
        if (  ( preQualityCutMask[p][PRE_Q_CUT_SAT_STRIPE]== 1 &&
                (   (-193.4150 < mapPos && mapPos < -171.8150 && canANITASeeStripe(plLon, plLat, sLat, -182.6150, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 1.7 )
                 || (-106.6350 < mapPos && mapPos <  -89.4350 && canANITASeeStripe(plLon, plLat, sLat,  -91.0350, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 2.2 )
                 || ( -46.3350 < mapPos && mapPos <  -33.4350 && canANITASeeStripe(plLon, plLat, sLat,  -39.8850, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 2.2 )
                 || (  59.3100 < mapPos && mapPos <   72.1100 && canANITASeeStripe(plLon, plLat, sLat,   65.7100, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 1.7 )
                 || (  52.4950 < mapPos && mapPos <   58.9850 && canANITASeeStripe(plLon, plLat, sLat,   55.7400, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 2.0 )
                 || (  69.4850 < mapPos && mapPos <   87.6850 && canANITASeeStripe(plLon, plLat, sLat,   78.5850, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 2.0 )
                 || ( 104.3225 < mapPos && mapPos <  109.8225 && canANITASeeStripe(plLon, plLat, sLat,  107.0725, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 2.0 )
                 || ( 124.5225 < mapPos && mapPos <  134.1225 && canANITASeeStripe(plLon, plLat, sLat,  129.3225, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 2.0 )
                )
              )
            ||( preQualityCutMask[p][PRE_Q_CUT_SAT_STRIPE]==-1 &&
               !(   (-193.4150 < mapPos && mapPos < -171.8150 && canANITASeeStripe(plLon, plLat, sLat, -182.6150, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 1.7 )
                 || (-106.6350 < mapPos && mapPos <  -89.4350 && canANITASeeStripe(plLon, plLat, sLat,  -91.0350, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 2.2 )
                 || ( -46.3350 < mapPos && mapPos <  -33.4350 && canANITASeeStripe(plLon, plLat, sLat,  -39.8850, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 2.2 )
                 || (  59.3100 < mapPos && mapPos <   72.1100 && canANITASeeStripe(plLon, plLat, sLat,   65.7100, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 1.7 )
                 || (  52.4950 < mapPos && mapPos <   58.9850 && canANITASeeStripe(plLon, plLat, sLat,   55.7400, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 2.0 )
                 || (  69.4850 < mapPos && mapPos <   87.6850 && canANITASeeStripe(plLon, plLat, sLat,   78.5850, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 2.0 )
                 || ( 104.3225 < mapPos && mapPos <  109.8225 && canANITASeeStripe(plLon, plLat, sLat,  107.0725, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 2.0 )
                 || ( 124.5225 < mapPos && mapPos <  134.1225 && canANITASeeStripe(plLon, plLat, sLat,  129.3225, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 2.0 )
                )
              )
           )
        {
          keep[p] = false;
          preQCuts[p][PRE_Q_CUT_SAT_STRIPE] = 1;
        }

        //cout << "Pol: " << p << " Phi: " <<  eventSummary->peak[p][0].phi << " PLLon: " << plLon << " Heading: " << gps->heading << endl;


        // C-Pol Peak Seperation Cut
        // -Added in from runAnalysisStage02.cxx

        float cPolPeakSepThres = 46.0;
        //float cPolPeakStrThres = 0.015;
        float cPolDist = cPolMinDist;
        //float minCircPeakVal = 

        if (  ( preQualityCutMask[p][PRE_Q_CUT_CPOL_SEP]== 1 && cPolDist >  cPolPeakSepThres) 
            ||( preQualityCutMask[p][PRE_Q_CUT_CPOL_SEP]==-1 && cPolDist <= cPolPeakSepThres)
           ) 
        {
          keep[p] = false;
          preQCuts[p][PRE_Q_CUT_CPOL_SEP] = 1;
        }

      
        // C-Pol Peak Strength Cut
        // -Added in from runAnalysisStage02.cxx

        // ACTUALLY this cut DOES depend on H and V pol, alot, actually.  so this one will come later!
        /*

        if (  ( preQualityCutMask[p][PRE_Q_CUT_CPOL_STR]== 1 && minCircPeakVal < cPolPeakStrThres)
            ||( preQualityCutMask[p][PRE_Q_CUT_CPOL_STR]==-1 && minCircPeakVal >=cPolPeakStrThres)
           )
        {
          keep[p] = false;
          preQCuts[p][PRE_Q_CUT_CPOL_STR] = 1;
        }

        */

        // L3 Trigger (OR)
        // -This is the '1' qCut, make into an or

        if (    (    (qualityCutMask[p][PRE_Q_CUT_NUM_PHISECTOR]== 1 || qualityCutMask[p][Q_CUT_NUM_PHISECTOR]== 1) 
                  && ( (qCuts[0][Q_CUT_NUM_PHISECTOR]==1) && (qCuts[1][Q_CUT_NUM_PHISECTOR]==1) ))
             || (    (qualityCutMask[p][PRE_Q_CUT_NUM_PHISECTOR]==-1 || qualityCutMask[p][Q_CUT_NUM_PHISECTOR]==-1) 
                  && ( (qCuts[0][Q_CUT_NUM_PHISECTOR]==0) || (qCuts[1][Q_CUT_NUM_PHISECTOR]==0) ))
           )
        {
          keep[p] = false;
          preQCuts[p][PRE_Q_CUT_NUM_PHISECTOR] = 1;
        }


        // now the preQ cuts are applied, lets count up values for 
        //  -count
        //  -orderedCount
        //  -uniqueCount
        //  -totalCount

        for (int k=0; k<NUM_PRE_Q_CUTS; ++k) {
          if (preQCuts[p][k]>0) { preQualityCutCounts[p][k] += eventWeight[p]; ++numPreQCuts[p];}
        }
        for (int k=0; k<NUM_PRE_Q_CUTS; ++k) {
          if (preQCuts[p][k]>0 && numPreQCuts[p]==1) { preQualityCutUniqueCounts[p][k] += eventWeight[p];}
        }
        if (simulationMode && (cohSnr2H > 5.0 || cohSnr2V > 5.0) ) {
          for (int k=0; k<NUM_PRE_Q_CUTS; ++k) {
            if (preQCuts[p][k]>0) { preQualityCutCounts5[p][k] += eventWeight[p]; ++numPreQCuts5[p];}
          }
          for (int k=0; k<NUM_PRE_Q_CUTS; ++k) {
            if (preQCuts[p][k]>0 && numPreQCuts5[p]==1) { preQualityCutUniqueCounts5[p][k] += eventWeight[p];}
          }
        }
        bool done = false;
        for (int k=0; k<preQCutOrder.size() && !done; ++k) {
          int n = preQCutOrder[k];
          if (preQCuts[p][n]>0) {
            preQualityCutOrderedCounts[p][n] += eventWeight[p];
            preQualityTotalOrderedCount[p] += eventWeight[p];
            done = true;}
        }
        if (simulationMode && (cohSnr2H > 5.0 || cohSnr2V > 5.0) ) {
          done = false;
          for (int k=0; k<preQCutOrder.size() && !done; ++k) {
            int n = preQCutOrder[k];
            if (preQCuts[p][n]>0) {
              preQualityCutOrderedCounts5[p][n] += eventWeight[p];
              preQualityTotalOrderedCount5[p] += eventWeight[p];
              done = true;}
          }
        }

        /*
        cout << "Event: " << header->eventNumber << ", PreQCuts: ";
        for(int k=0; k<preQCutOrder.size(); k++)
        {
          cout << preQCuts[p][k] << " ";
        }
        cout << endl;
        */

        // After preQCuts are done, we need to apply the L3 Trigger 
        // cut seperatly to each pol, and save as the number of events
        // passing our preQCuts, preQEventCounts. 
        keepPreQ[p] = keep[p];

        if (    ((qualityCutMask[p][Q_CUT_NUM_PHISECTOR] ==  1) && (qCuts[p][Q_CUT_NUM_PHISECTOR] == 1))
             || ((qualityCutMask[p][Q_CUT_NUM_PHISECTOR] == -1) && (qCuts[p][Q_CUT_NUM_PHISECTOR] == 0))
           )
        {
          keepPreQ[p] = false;
          keep[p] = false;
        }
        

        if (keep[p]) {  preQEventCounts[p] += eventWeight[p]; }        
        if (keep[p] && simulationMode && cohSnr2 > 5.0) { preQEventCounts5[p] += eventWeight[p]; }


        for (int k=0; k<NUM_Q_CUTS; ++k) {  
          //vcout << " qcut value [" << p << "][" << k << "] is " << qCuts[p][k] << endl;
          //vcout << "   qcut mask value is " << qualityCutMask[p][k] << endl;
          keep[p] &= ((qualityCutMask[p][k]==0) || ((qualityCutMask[p][k]==1) && (qCuts[p][k]==0)) || ((qualityCutMask[p][k]==-1) && (qCuts[p][k]==1)));
        }
    
        for (int k=0; k<NUM_Q_CUTS; ++k) {        
          if (qCuts[p][k]>0 && keepPreQ[p]) { qualityCutCounts[p][k] += eventWeight[p]; ++numQCuts[p];}
        }
        for (int k=0; k<NUM_Q_CUTS; ++k) {
          if (qCuts[p][k]>0 && numQCuts[p]==1 && keepPreQ[p]) { qualityCutUniqueCounts[p][k] += eventWeight[p];}
        }  
        if (simulationMode && cohSnr2 > 5.0) {
          for (int k=0; k<NUM_Q_CUTS; ++k) {
            if (qCuts[p][k]>0 && keepPreQ[p]) { qualityCutCounts5[p][k] += eventWeight[p]; ++numQCuts5[p];}
          }
          for (int k=0; k<NUM_Q_CUTS; ++k) {
            if (qCuts[p][k]>0 && numQCuts5[p]==1 && keepPreQ[p]) { qualityCutUniqueCounts5[p][k] += eventWeight[p];}
          }
        }
        done = false;  
        for (int k=0; k<qCutOrder.size() && !done; ++k) {
          int n = qCutOrder[k];
          if (qCuts[p][n]>0 && keepPreQ[p]) {
            qualityCutOrderedCounts[p][n] += eventWeight[p];
            qualityTotalOrderedCount[p] += eventWeight[p]; 
            done = true;}
        }  
        if (simulationMode && cohSnr2 > 5.0) {
          done = false;
          for (int k=0; k<qCutOrder.size() && !done; ++k) {
            int n = qCutOrder[k];
            if (qCuts[p][n]>0 && keepPreQ[p]) {
              qualityCutOrderedCounts5[p][n] += eventWeight[p];
              qualityTotalOrderedCount5[p] += eventWeight[p]; 
              done = true;}
          }
        }

        vbprintf("  keep-it value[%i] is %i \n", p, keep[p]);
        //vbprintf ("80 current directory is %s \n", gDirectory->GetName());          

        // look at the peaks and warn if not in order by strength   
        /*
        bool peakOrderViol = false;      
        if (false) {    
          for (int i=0; i<eventSummary->nPeaks[p]; ++i) {
            vcout << " peak " << pChars[p] << ":" << i << " " << eventSummary->peak[p][i].value << endl;
            if (i>0)  {
              if (eventSummary->peak[p][i].value > eventSummary->peak[p][i-1].value) {
                //cout << "warning: event " << header->eventNumber << "   peak value order violation at peak " << i << "," << i-1 << "   " 
                //        << eventSummary->peak[p][i].value << ", " << eventSummary->peak[p][i-1].value << endl;
                if (p==0) {
                  lprintf("warning: event %8i %c  peak value order violation at peak %1i,%1i   (%5.2f,%5.2f):%4f           (%5.2f,%5.2f):%4f      wais is (%5.2f,%5.2f) \n", 
                        header->eventNumber, pChars[p], i-1, i, 
                        eventSummary->peak[p][i-1].phi, eventSummary->peak[p][i-1].theta, eventSummary->peak[p][i-1].value,
                        eventSummary->peak[p][i].phi, eventSummary->peak[p][i].theta, eventSummary->peak[p][i].value,
                        eventSummary->wais.phi, eventSummary->wais.theta);
                }
                peakOrderViol = true;
              }
            }
          }
          if (peakOrderViol) ++peakViolationCount;
        }
        */
        // analysis cuts
        //vbprintf ("82 current directory is %s \n", gDirectory->GetName());          
        float ptDiffPhi = 0;
        float ptDiffTheta = 0;
        if (useTrigNsWindow) {
          ptDiffPhi = mainPeak.phi - calSource.phi;
          while (ptDiffPhi > 180) ptDiffPhi -= 360;
          while (ptDiffPhi <= -180) ptDiffPhi += 360;
          ptDiffTheta = mainPeak.theta - calSource.theta;
          //if (p==0) {
          if (true) {
            vbprintf("event %8i:%c pointing result is (%6.3f,%6.3f)     cal pulser is at (%6.3f,%6.3f)      difference (%6.3f,%6.3f) \n",
                    header->eventNumber, pChars[p], mainPeak.phi, mainPeak.theta, calSource.phi, calSource.theta, ptDiffPhi, ptDiffTheta);
          }
        }                      

        //int* thisACuts = (p==0) ? &aCuts[0][0] : &aCuts[1][0];
        vcout << "    aCuts[" << p << "] "; for (int k=0; k<NUM_A_CUTS; ++k) {vcout << aCuts[p][k] << " ";} vcout << "  total " << numACuts[p] << endl;
        for (int k=0; k<NUM_A_CUTS; ++k) {
          if (aCuts[p][k]>0 && numACuts[p]==1 && keep[p]) {analysisCutUniqueCounts[p][k] += eventWeight[p];}
          if (simulationMode && cohSnr2 > 5.0) {
            if (aCuts[p][k]>0 && numACuts5[p]==1 && keep[p]) {analysisCutUniqueCounts5[p][k] += eventWeight[p];}
          }
        }
        
        // arbitrary circle cut on pl-north sky
        
        float pnPhi = fmod((mainPeak.phi-gps->heading+360),360);
        
        // solar direct cut
        float sunDistPhi = mainPeak.phi - eventSummary->sun.phi;        
        float sunDistTheta = -mainPeak.theta - (-eventSummary->sun.theta);   
        while (sunDistPhi > 180) sunDistPhi -= 360;
        while (sunDistPhi < -180) sunDistPhi += 360;
        vbprintf(" peak position %6.3f, %6.3f \n", mainPeak.phi, mainPeak.theta);
        vbprintf(" sun position %6.3f, %6.3f \n", eventSummary->sun.phi, eventSummary->sun.theta);
        aCuts[p][A_CUT_SOLAR_DIRECT] = ((abs(sunDistPhi) < sunCutDistPhi) && (abs(sunDistTheta) <  sunCutDistTheta)) ? 1 : 0;  
        vbprintf("  solar direct cut %i \n", aCuts[p][A_CUT_SOLAR_DIRECT]);
        vcout << "    aCuts[" << p << "] "; for (int k=0; k<NUM_A_CUTS; ++k) {vcout << aCuts[p][k] << " ";} vcout << "  total " << numACuts[p] << endl;
        for (int k=0; k<NUM_A_CUTS; ++k) {
          if (aCuts[p][k]>0 && numACuts[p]==1 && keep[p]) {analysisCutUniqueCounts[p][k] += eventWeight[p];}
          if (simulationMode && cohSnr2 > 5.0) {
            if (aCuts[p][k]>0 && numACuts5[p]==1 && keep[p]) {analysisCutUniqueCounts5[p][k] += eventWeight[p];}
          }
        }

        // solar reflected cut
        float reflTheta = reflectionAngle(gps->altitude/1000.0, -mainPeak.theta);
        float reflPhi = mainPeak.phi;
        float sunReflDistPhi = eventSummary->sun.phi - reflPhi;        
        float sunReflDistTheta = -eventSummary->sun.theta - reflTheta;
        while (sunReflDistPhi > 180) sunReflDistPhi -= 360;
        while (sunReflDistPhi <= -180) sunReflDistPhi += 360;
        cout << "Phi Dist from sun: " << sunReflDistPhi << " Theta Dist from sun: " << sunReflDistTheta << endl;
        aCuts[p][A_CUT_SOLAR_REFL] = ((abs(sunReflDistPhi) < sunReflCutDistPhi) && (abs(sunReflDistTheta) < sunReflCutDistTheta)) ? 1 : 0;  
        vcout << "  solar reflection cut " << aCuts[p][A_CUT_SOLAR_REFL] << endl;
        vcout << "    aCuts[" << p << "] "; for (int k=0; k<NUM_A_CUTS; ++k) {vcout << aCuts[p][k] << " ";} vcout << "  total " << numACuts[p] << endl;
        for (int k=0; k<NUM_A_CUTS; ++k) {
          if (aCuts[p][k]>0 && numACuts[p]==1 && keep[p]) {analysisCutUniqueCounts[p][k] += eventWeight[p];}
          if (simulationMode && cohSnr2 > 5.0) {
            if (aCuts[p][k]>0 && numACuts5[p]==1 && keep[p]) {analysisCutUniqueCounts5[p][k] += eventWeight[p];}
          }
        }
                
        // reconstruct to continent cut
        //  TODO distinguish sea ice from continent  

        aCuts[p][A_CUT_CONTINENT] = (mainPeak.latitude == -9999) ? 1 : 0;
        vcout << "    reconstruction latitude = " << mainPeak.latitude << endl;
        vcout << "  reconstruct to continent cut: " << aCuts[p][A_CUT_CONTINENT] << endl;
        vcout << "    aCuts[" << p << "] "; for (int k=0; k<NUM_A_CUTS; ++k) {vcout << aCuts[p][k] << " ";} vcout << "  total " << numACuts[p] << endl;
        for (int k=0; k<NUM_A_CUTS; ++k) {
          if (aCuts[p][k]>0 && numACuts[p]==1 && keep[p]) {analysisCutUniqueCounts[p][k] += eventWeight[p];}
          if (simulationMode && cohSnr2 > 5.0) {
            if (aCuts[p][k]>0 && numACuts5[p]==1 && keep[p]) {analysisCutUniqueCounts5[p][k] += eventWeight[p];}
          }
        }
        // elevation angle of interferometery peak
        if ( eventSummary->eventNumber == 15717147 ) {
          cout << " theta: " << -mainPeak.theta << " minAngle: " << minPeakTheta << endl;
          aCuts[p][A_CUT_THETA] = (-mainPeak.theta < minPeakTheta-10) || (-mainPeak.theta > maxPeakTheta) ? 1 : 0;
        } else {
          aCuts[p][A_CUT_THETA] = (-mainPeak.theta < minPeakTheta) || (-mainPeak.theta > maxPeakTheta) ? 1 : 0;
        }
        vcout << "  peak elevation angle cut: " << aCuts[p][A_CUT_THETA] << endl;
        vcout << "    aCuts[" << p << "] "; for (int k=0; k<NUM_A_CUTS; ++k) {vcout << aCuts[p][k] << " ";} vcout << "  total " << numACuts[p] << endl;
        for (int k=0; k<NUM_A_CUTS; ++k) {
          if (aCuts[p][k]>0 && numACuts[p]==1 && keep[p]) {analysisCutUniqueCounts[p][k] += eventWeight[p];}
          if (simulationMode && cohSnr2 > 5.0) {
            if (aCuts[p][k]>0 && numACuts5[p]==1 && keep[p]) {analysisCutUniqueCounts5[p][k] += eventWeight[p];}
          }
        }
        // angle from peak to trigger center
        aCuts[p][A_CUT_HW_ANGLE] = abs(mainPeak.hwAngle) > maxHwAngle ? 1 : 0;
        vcout << "  hardware angle cut: " << aCuts[p][A_CUT_HW_ANGLE] << "  (" << mainPeak.hwAngle << ")" << endl;
        vcout << "    aCuts[" << p << "] "; for (int k=0; k<NUM_A_CUTS; ++k) {vcout << aCuts[p][k] << " ";} vcout << "  total " << numACuts[p] << endl;
        for (int k=0; k<NUM_A_CUTS; ++k) {
          if (aCuts[p][k]>0 && numACuts[p]==1 && keep[p]) {analysisCutUniqueCounts[p][k] += eventWeight[p];}
          if (simulationMode && cohSnr2 > 5.0) {
            if (aCuts[p][k]>0 && numACuts5[p]==1 && keep[p]) {analysisCutUniqueCounts5[p][k] += eventWeight[p];}
          }
        }
        // L3 triggering direction
        bool foundIt = false;
        int peakPhisec = (int)fmod((mainPeak.phi + 45 + 11.25), 360)/22.5;
        vbprintf(" peak phi-sector is %i \n", peakPhisec);
        for (int psd = -trigPhisecProx; psd <= trigPhisecProx && !foundIt; ++psd) {
          int thisPs = (peakPhisec + psd) % 16;
          if (thisPs < 0) thisPs += 16;
          if(header->isInL3Pattern(thisPs, (AnitaPol::AnitaPol_t)p)) {
            foundIt = true;
          }
        }        
        aCuts[p][A_CUT_L3_TRIG_DIR] = foundIt ? 0 : 1;
        vcout << "  l3 trigger direction cut: " << aCuts[p][A_CUT_L3_TRIG_DIR] << endl;
        vcout << "    aCuts[" << p << "] "; for (int k=0; k<NUM_A_CUTS; ++k) {vcout << aCuts[p][k] << " ";} vcout << "  total " << numACuts[p] << endl;
        for (int k=0; k<NUM_A_CUTS; ++k) {
          if (aCuts[p][k]>0 && numACuts[p]==1 && keep[p]) {analysisCutUniqueCounts[p][k] += eventWeight[p];}
          if (simulationMode && cohSnr2 > 5.0) {
            if (aCuts[p][k]>0 && numACuts5[p]==1 && keep[p]) {analysisCutUniqueCounts5[p][k] += eventWeight[p];}
          }
        }
        // fraction linear polarization
        //AnitaEventSummary::WaveformInfo reconWfI = eventSummary->coherent[p][maxPeakIndex[p][0]];
        AnitaEventSummary::WaveformInfo reconWfI = eventSummary->coherent[p][0];
        double sI, sQ, sU, sV;
        sI = reconWfI.I;
        if (circPol) {
          sQ = reconWfI.U;
          sU = -reconWfI.V;
          sV = -reconWfI.Q;
        } else {
          sQ = reconWfI.Q;
          sU = reconWfI.U;
          sV = reconWfI.V;
        }
        //vbprintf ("84 current directory is %s \n", gDirectory->GetName());          

        float thisLinPolAll = 0;
        float thisLinPolTrig = 0;
        float thisCircPolAll = 0;
        float thisCircPolTrig = 0;
        float numAntsTrig = 0;
        for (int a=0; a<NUM_SEAVEYS; ++a) {
          float stokesI = 0;
          float stokesQ = 0;
          float stokesU = 0;
          float stokesV = 0;        
          AnitaEventSummary::WaveformInfo thisWFI = eventSummary->inputWfRaw[p][a];
          stokesI = thisWFI.I;
          if (!circPol) {
            stokesQ = thisWFI.Q;
            stokesU = thisWFI.U;
            stokesV = thisWFI.V;
          } else {
            stokesQ = thisWFI.U;
            stokesU = -thisWFI.V;
            stokesV = -thisWFI.Q;
          }
          if (header->isInL3Pattern(geom->getPhiFromAnt(a), (AnitaPol::AnitaPol_t)p)) {
            thisLinPolTrig += sqrt(stokesQ*stokesQ + stokesU*stokesU)/stokesI;
            thisCircPolTrig += sqrt(stokesV*stokesV)/stokesI;
            numAntsTrig += 1.0;
          }
          thisLinPolAll += sqrt(stokesQ*stokesQ + stokesU*stokesU)/stokesI;
          thisCircPolAll += sqrt(stokesV*stokesV)/stokesI;
        }
        thisLinPolAll /= (float)NUM_SEAVEYS;
        thisLinPolTrig /= numAntsTrig;
        thisCircPolAll /= (float)NUM_SEAVEYS;
        thisCircPolTrig /= numAntsTrig;
        //printf("event %i I=%f, Q=%f, U=%f, V=%f \n", header->eventNumber, stokesI, stokesQ, stokesU, stokesV);

        // careful now airstream driver  : calculate linear polarization from reconstruction or from raw waveform average?
        if (polFracMethod==0) {
          linPolFrac[p] = sqrt(sQ*sQ + sU*sU)/sI;
          circPolFrac[p] = sqrt(sV*sV)/sI;
        } else if (polFracMethod==1) {
          linPolFrac[p] = thisLinPolAll;
          circPolFrac[p] = thisCircPolAll;
        } else if (polFracMethod==2) {
          linPolFrac[p] = thisLinPolTrig;
          circPolFrac[p] = thisCircPolTrig;
        }
        anyPolFrac[p] = sqrt(sQ*sQ + sU*sU + sV*sV)/sI;
        //vbprintf ("85 current directory is %s \n", gDirectory->GetName());          

        aCuts[p][A_CUT_LINPOL_FRAC] = (linPolFrac[p] < minLinPolFrac) ? 1 : 0;
        vbprintf("linear pol fraction cut: %i \n", aCuts[p][A_CUT_LINPOL_FRAC]);
        //vbprintf ("85.01 current directory is %s \n", gDirectory->GetName());          
        // fraction any polarization
        aCuts[p][A_CUT_ANYPOL_FRAC] = (anyPolFrac[p] < minAnyPolFrac) ? 1 : 0;
        vbprintf("any pol fraction cut: %i \n", aCuts[p][A_CUT_ANYPOL_FRAC]);
        //vbprintf ("85.02 current directory is %s \n", gDirectory->GetName());          
        // fraction circ polarization
        aCuts[p][A_CUT_CIRCPOL_FRAC] = (circPolFrac[p] > maxCircPolFrac) ? 1 : 0;
        vbprintf("circ pol fraction cut: %i \n", aCuts[p][A_CUT_CIRCPOL_FRAC]);
        //vbprintf ("85.03 current directory is %s \n", gDirectory->GetName());          
        // peak ratio cut TODO make this accommodate out-of-order peaks
        //peakRatio[p] = eventSummary->peak[p][maxPeakIndex[p][1]].value / eventSummary->peak[p][maxPeakIndex[p][0]].value;
        peakRatio[p] = eventSummary->peak[p][1].value / eventSummary->peak[p][0].value;
        aCuts[p][A_CUT_PEAK_RATIO] = (peakRatio[p] > peakRatioThreshold) ? 1 : 0;
        vbprintf("1st/2nd peak ratio cut: %i \n", aCuts[p][A_CUT_PEAK_RATIO]);
        //vbprintf ("85.04 current directory is %s \n", gDirectory->GetName());          
        // snr cut
        aCuts[p][A_CUT_SNR] = (reconWfI.snr2 < minSnr) ? 1 : 0;
        vbprintf("snr fraction cut: %i \n", aCuts[p][A_CUT_SNR]);        
        // vbprintf ("85.05 current directory is %s \n", gDirectory->GetName());          
        // correlation peak value
        aCuts[p][A_CUT_CORR_PEAK] = (mainPeak.value < corrPeakThresh) ? 1 : 0;
        vbprintf("correlation peak cut: %i (%f) (%f) \n", aCuts[p][A_CUT_CORR_PEAK], mainPeak.value, corrPeakThresh);
        //vbprintf ("85.06 current directory is %s \n", gDirectory->GetName());          
        // hilbert peak value
        aCuts[p][A_CUT_HILBERT_PEAK] = (reconWfI.peakHilbert < hilbPeakThresh) ? 1 : 0;;
        vbprintf("hilbert peak cut: %i \n", aCuts[p][A_CUT_HILBERT_PEAK]);
        //vbprintf ("85.07 current directory is %s \n", gDirectory->GetName());          
        // corr_hilbert linear discriminant cut
        float diagVal = -corrPeakDiagThresh/hilbPeakDiagThresh*reconWfI.peakHilbert + corrPeakDiagThresh - mainPeak.value;
        aCuts[p][A_CUT_DIAGONAL] = (diagVal > 0) ? 1 : 0;
        vbprintf("linear discriminant cut: %i \n", aCuts[p][A_CUT_DIAGONAL]);
        //vbprintf ("85.08 current directory is %s \n", gDirectory->GetName());          
        // corr_snr linear discriminant cut (just calculate for now)
        ldCutVal[p] = -ldSlope*mainPeak.value + reconWfI.snr2;
        // WAIS pulser cut
        //int nsDiff = (int)((float)eventSummary->wais.distance/0.3 - (float)header->triggerTimeNs);
        //vcout << "  WAIS nsDiff = " << nsDiff << endl;
        //vbprintf ("85.085 current directory is %s \n", gDirectory->GetName());          
        bool waisPulse = UCorrelator::isWAISHPol(gps, header);
        vbprintf ("85.086 current directory is %s \n", gDirectory->GetName());          
        aCuts[p][A_CUT_WAIS] = waisPulse;
        vcout << "WAIS pulser cut " << aCuts[p][A_CUT_WAIS] << endl;
        // TODO LDB pulser cut
        bool ldbPulse = UCorrelator::isLDB(header);
        vbprintf ("85.087 current directory is %s \n", gDirectory->GetName());          
        aCuts[p][A_CUT_LDB] = ldbPulse;
        vcout << "LDB pulser cut " << aCuts[p][A_CUT_LDB] << endl;
        // TODO HiCal pulser cut?
        //vbprintf ("85.1 current directory is %s \n", gDirectory->GetName());          
        
        // south cut     fail if further-than-threshold away from south
        aCuts[p][A_CUT_SOUTH] = abs(180.0 - pnPhi) > southCutDist ? 1 : 0;
        vbprintf("south cut: %i \n", aCuts[p][A_CUT_SOUTH]);
        
        // wedge cut (for obtaining satellite-free (ha) noise samples
        //   why am I using sky events?  Won't the ground give more characteristic thermal?
        float yLim = (pnPhi>180) ? wedgeCutMid + (wedgeCutEdge-wedgeCutMid)/180.0 * (pnPhi - 180.0) : 
                wedgeCutEdge - (wedgeCutEdge-wedgeCutMid)/180.0 * (pnPhi);
        aCuts[p][A_CUT_WEDGE] = (-mainPeak.theta > yLim) ? 1 : 0;
        vbprintf("wedge cut: %i \n", aCuts[p][A_CUT_WEDGE]);
        
        //printf("sky rectangle cut: event is at (%6.3f,%6.3f), rectangle is (%6.3f,%6.3f),(%6.3f,%6.3f) \n",
        //          pnPhi, -mainPeak.theta, skyRectCut[0], skyRectCut[1], skyRectCut[2], skyRectCut[3]);
        // rectangle cut, for generally studying regions
        aCuts[p][A_CUT_RECT] = (pnPhi >= skyRectCut[0] && pnPhi <= skyRectCut[2] && -mainPeak.theta >= skyRectCut[1] 
                && -mainPeak.theta <= skyRectCut[3]) ? 0 : 1;
        vbprintf("rectangle cut: %i \n", aCuts[p][A_CUT_RECT]);
        
        deadTimeFrac[p] = (float)(header->deadTime) / 65536.0;
        //printf("  deadtime from header %i  frac %f \n", header->deadTime, deadTimeFrac[p]);
        aCuts[p][A_CUT_DEADTIME] = (deadTimeFrac[p] > deadTimeThreshold) ? 1 : 0;
        vbprintf("deadtime cut: %i \n", aCuts[p][A_CUT_DEADTIME]);
        //vbprintf ("85.2 current directory is %s \n", gDirectory->GetName());          
        
        // continent circle cut
        calPulseDist[p] = -1;
        if (aCuts[p][A_CUT_CONTINENT]==0) {
          //double evNo, evEa;
          calPulseDist[p] = sqrt((ksEa-plEa)*(ksEa-plEa) + (ksNo-plNo)*(ksNo-plNo));
          vbprintf("dist from payload to cal pulser is %f \n", calPulseDist[p]);
          double circEa, circNo;
          bedmap->LonLattoEaNo(HOT_SPOT_1[1], HOT_SPOT_1[0], circEa, circNo);
          float dist = sqrt((srcEa[p]-circEa)*(srcEa[p]-circEa)+(srcNo[p]-circNo)*(srcNo[p]-circNo));
          vbprintf("dist from event to cut circle center is %f \n", dist);
          aCuts[p][A_CUT_CONT_CIRC] = (dist < 50000) ? 1: 0;
        }
        //vbprintf ("86 current directory is %s \n", gDirectory->GetName());          

        // simulated neutrino energy
        aCuts[p][A_CUT_SIM_E_NU] = 0;
        if (simulationMode) {
          aCuts[p][A_CUT_SIM_E_NU] = (pnu < simENuCut[0] || pnu > simENuCut[1]) ? 1 : 0;
        }     
        vbprintf("simulated neutrino energy cut: %i \n", aCuts[p][A_CUT_SIM_E_NU]);
        /////////////////////////   end of cut calculations   /////////////
        /*
        // accumulate for snr-binned hists
        //int snrIndex = (int)(reconWfI.snr);
        int snrIndex = (int)(maxWfSnr[p]);
        if (snrIndex > 19) snrIndex = 19;
        if (snrIndex < 0) snrIndex = 0;
        errPhiBySnr[p][snrIndex] += (float)pow(ptDiffPhi, 2);
        errThetaBySnr[p][snrIndex] += (float)pow(ptDiffTheta, 2);
        float dist = calSource.distance/1000.0;
        eventCountsSnr[p][snrIndex] += eventWeight[p];
        //printf("filling distSnrAll: %f, %f, %f \n", reconWfI.snr, dist, eventWeight[p]);
        //distSnrAll[p]->Fill(reconWfI.snr, dist, eventWeight[p]);
        distSnrAll[p]->Fill(maxWfSnr[p], dist, eventWeight[p]);
        //totSnrPFFrac[p]->Fill(reconWfI.snr, filterPowerFrac_0[p], eventWeight[p]);
        totSnrPFFrac[p]->Fill(maxWfSnr[p], filterPowerFrac_1[p], eventWeight[p]);
        */
        //firstPol = false;
      }
      //vbprintf ("90 current directory is %s \n", gDirectory->GetName());          

      // now's the time on schprockets when we process the circPol entry, if requested
      //float cPolDist = 0;
      float cPolLesserPeak = 0;
      float cPolMinDist = 0;
      float cPolPeak[2] = {0};
      float cPolWeight = 0; 
      float cPolPhi[2] = {0};
      float cPolTheta[2] = {0};
      float cPolSnr2[2] = {0};
      if (doCircPol) {   // CPol processing was requested
        resultTree->SetBranchAddress("header", &cPolHeader);
        resultTree->SetBranchAddress("eventSummary", &cPolEventSummary);
        if (simulationMode) {
          resultTree->SetBranchAddress("weight", &cPolWeight);
        }
        //vbprintf ("100 current directory is %s \n", gDirectory->GetName());          
      
                // retain the info from the linpol entry
        // read the next entry; if it's not a circPol, warn
        int thisEventNum = header->eventNumber;
        if (e >= resultTree->GetEntries()-1) {
          lprintf("warning: no circPol entry found for event %i \n", thisEventNum);
        } else {
          resultTree->GetEntry(e+1);
          if (circPol!=1 || header->eventNumber != thisEventNum ) {
            lprintf("warning: no circPol entry found for event %i (next event=%i circpol=%i) \n", thisEventNum, header->eventNumber, circPol);   // I HATE repeating code, even one line
          } else {  // we are good to go for CPol analysis stage
            if (cPolEventSummary->nPeaks[0] <2 || cPolEventSummary->nPeaks[1] < 2) {
              printf("warning: <2 peaks in circPol interferometry result for event %i: skipping \n", thisEventNum);
              continue;
            } else {   // really do circ pol
              vbprintf("-----doing circ pol analysis \n");
            // apply the distance-between-peaks cut
              //AnitaEventSummary::PointingHypothesis peak0 = cPolEventSummary->peak[0][maxPeakIndex[0][0]];
              //AnitaEventSummary::PointingHypothesis peak1 = cPolEventSummary->peak[1][maxPeakIndex[1][0]];

              AnitaEventSummary::WaveformInfo cPolWfi[2]; for (int k=0; k<2; ++k) {cPolWfi[k] = cPolEventSummary->coherent[k][maxPeakIndex[k][0]];}
              for (int k=0; k<2; ++k) {
                cPolPeak[k] = cPolEventSummary->peak[k][0].value;
                cPolPhi[k] = cPolEventSummary->peak[k][0].phi;
                cPolTheta[k] = cPolEventSummary->peak[k][0].theta;
                cPolSnr2[k] = cPolWfi[k].snr2;
              }
              
              // find the minimum distance between first L(R) peak and each R(L) peak 
              vector<float> peakDists(0);
              AnitaEventSummary::PointingHypothesis peak0 = cPolEventSummary->peak[0][0];
              for (int kr=0; kr<cPolNumPeaks; ++kr) {
                AnitaEventSummary::PointingHypothesis peak1 = cPolEventSummary->peak[1][kr];
                float cPolDist = vectorAngle(peak0.phi, peak1.phi, peak0.theta, peak1.theta);
                peakDists.push_back(cPolDist);
                vbprintf(" kr=%i, lPeak=(%f,%f) rPeak=(%f,%f) cPolDist = %f \n", kr, cPolEventSummary->peak[0][0].phi, cPolEventSummary->peak[0][0].theta, 
                        cPolEventSummary->peak[1][kr].phi, cPolEventSummary->peak[1][kr]. theta, cPolDist);
              }
              AnitaEventSummary::PointingHypothesis peak1 = cPolEventSummary->peak[1][0];
              for (int kl=1; kl<cPolNumPeaks; ++kl) {
                AnitaEventSummary::PointingHypothesis peak1 = cPolEventSummary->peak[0][kl];
                float cPolDist = vectorAngle(peak0.phi, peak1.phi, peak0.theta, peak1.theta);
                peakDists.push_back(cPolDist);
                vbprintf(" kl=%i, lPeak=(%f,%f) rPeak=(%f,%f) cPolDist = %f \n", kl, cPolEventSummary->peak[0][kl].phi, cPolEventSummary->peak[0][kl].theta, 
                        cPolEventSummary->peak[1][0].phi, cPolEventSummary->peak[1][0]. theta, cPolDist);
              }
 
              cPolMinDist = *min_element(peakDists.begin(), peakDists.end());
              vbprintf(" min peak separation is %f \n", cPolMinDist);
              // apply the L/R peak separation cut
              for (int pol=0; pol<2; ++pol) {
                aCuts[pol][A_CUT_CPOL_PEAK_SEPARATION] = (cPolMinDist > cPolPeakSeparationThreshold) ? 1 : 0;
                vcout << "  cpol peak separation: " << aCuts[pol][A_CUT_CPOL_PEAK_SEPARATION] << endl;
              }                
                  
              // apply the linpol/circpol lesser peak cut
              cPolLesserPeak = min(peak0.value, peak1.value) / (peak0.value + peak1.value);
              for (int pol=0; pol<2; ++pol) {
                aCuts[pol][A_CUT_CPOL_LESSER_PEAK] = (cPolLesserPeak < cPolLesserPeakThreshold) ? 1 : 0;
                vcout << "  cpol lesser peak: " << aCuts[pol][A_CUT_CPOL_LESSER_PEAK] << endl;
              }          
              
              // calculate the LPol and RPol peak strengths, relative to the H and V peaks
              
              for (int lp=0; lp<2; ++lp) {
                
              
              }
              // geostat satellite cut on Lpol localization 
              // TODO this is the left-hand wedge only; also need to do the right-hand wedge
              float pnPhi = fmod((peak0.phi-gps->heading+360),360);
              float gs1 = pnPhi/50.0 - 3.0;
              float gs2 = -pnPhi/5.0 + 15.0;
              aCuts[0][A_CUT_GEOSTAT_DIRECT] = (-peak0.theta > gs1 && -peak0.theta < gs2) ? 1 : 0;   // fail if in geostat zone
              aCuts[1][A_CUT_GEOSTAT_DIRECT] = aCuts[0][A_CUT_GEOSTAT_DIRECT];
              // apply the left-circular pol theta cut (stuff that reconstructs Lpol outside the theta region is likely satellite)
              aCuts[0][A_CUT_LPOL_THETA] = (-peak0.theta > minPeakTheta && -peak0.theta < maxPeakTheta) ? 0 : 1;
              // TODO geostationary reflected cut
              //    shoot the reflection and do the direct test on it
              
            }            
            ++e;   // intentional structure violation: this is also incremented in for-loop declaration
          }
        }
      } else {           // CPol was not requested
        // read the next entry; if it's a circPol, skip it
        if (e < resultTree->GetEntries()-1) {
          resultTree->GetEntry(e+1);
          if (circPol==1) ++e;  // intentional structure violation: this is also incremented in for-loop declaration
        }
      }

      for (int p=0; p<2; ++p) {
        //AnitaEventSummary::PointingHypothesis mainPeak = eventSummary->peak[p][maxPeakIndex[p][0]];
        AnitaEventSummary::PointingHypothesis mainPeak = eventSummary->peak[p][0];
        float cohSnr2 = eventSummary->coherent[p][maxPeakIndex[p][0]].snr2;

        // count cuts as requested by the user   
        for (int k=0; k<aCutOrderStage1Limit; ++k) {        
          int n=aCutOrderStage1[k];
          if (aCuts[p][n]>0 && keep[p]) {
            analysisCutCounts[p][n] += eventWeight[p]; 
            analysisCutEventList[p][n].push_back(header->eventNumber);
            ++numACuts[p];
          }
        } 
        if (simulationMode && cohSnr2 > 5.0){
          for (int k=0; k<aCutOrderStage1Limit; ++k) {
            int n=aCutOrderStage1[k];
            if (aCuts[p][n]>0 && keep[p]) {
              analysisCutCounts5[p][n] += eventWeight[p];
              //analysisCutEventList[p][n].push_back(header->eventNumber);
              ++numACuts5[p];
            }
          }
        }
        
        if (keep[p]) {
          vbprintf(" keeping event %i %cpol after quality cuts ----------------- \n", header->eventNumber, pChars[p]);
          qEventCounts[p] += eventWeight[p];
          if (simulationMode && cohSnr2 > 5.0) {qEventCounts5[p] += eventWeight[p]; }
          vbprintf("  count is now %f \n", qEventCounts[p]);
        }        
        // enforce analysis cuts on plots and output as requested by the user      
        vcout << "    aCuts[" << p << "] "; for (int k=0; k<NUM_A_CUTS; ++k) {vcout << aCuts[p][k] << " ";} vcout << "  total " << numACuts[p] << endl;
        for (int k=0; k<NUM_A_CUTS; ++k) {
          if (aCuts[p][k]>0 && numACuts[p]==1 && keep[p]) {analysisCutUniqueCounts[p][k] += eventWeight[p];}
          if (simulationMode && cohSnr2 > 5.0) {         
            if (aCuts[p][k]>0 && numACuts5[p]==1 && keep[p]) {analysisCutUniqueCounts5[p][k] += eventWeight[p];}
          }
        }
        bool done = false; 
        bool passedCuts = true; 
        //for (int k=0; k<aCutOrderStage1.size(); ++k) {
        for (int k=0; k<aCutOrderStage1Limit; ++k) {
          int n = aCutOrderStage1[k];
          if (!done && aCuts[p][n]>0 && keep[p]) {
            analysisCutOrderedCounts[p][n] += eventWeight[p]; 
            analysisTotalOrderedCount[p] += eventWeight[p]; 
            done = true;}
          passedCuts &= (aCuts[p][n]==0);
        }
        done = false;
        if (simulationMode && cohSnr2 > 5.0) {
          for (int k=0; k<aCutOrderStage1Limit; ++k) {
            int n = aCutOrderStage1[k];
            if (!done && aCuts[p][n]>0 && keep[p]) {
              analysisCutOrderedCounts5[p][n] += eventWeight[p]; 
              analysisTotalOrderedCount5[p] += eventWeight[p]; 
              done = true;}
            //passedCuts &= (aCuts[p][n]==0);
          }
        }

        if (passedCuts) {
          float dist = calSource.distance/1000.0;
          //acceptsSnr[p][snrIndex] += eventWeight[p];
          //distSnrEff[p]->Fill(reconWfI.snr, dist, eventWeight[p]);
          distSnrEff[p]->Fill(maxWfSnr[p], dist, eventWeight[p]);
          //effSnrPFFrac[p]->Fill(reconWfI.snr, filterPowerFrac_0[p], eventWeight[p]);
          effSnrPFFrac[p]->Fill(maxWfSnr[p], filterPowerFrac_1[p], eventWeight[p]);
        }
        for (int k=0; k<NUM_A_CUTS; ++k) {
          keep[p] &= ((analysisCutMask[p][k]==0) || ((analysisCutMask[p][k]==1) && (aCuts[p][k]==0)) || ((analysisCutMask[p][k]==-1) && (aCuts[p][k]==1)));
        }

        vcout << "  keep-it value is " << keep[p] << " " << endl;
        if (keep[p]) {
          //newEventNumber = header->eventNumber;
          vbprintf("event %i lat/lon/alt = %f, %f, %f  \n", header->eventNumber, newLat[p], newLon[p], newAlt[p]);
          //moreResultsTree->Fill();
          aEventCounts[p] += eventWeight[p];
          if (listPassingEvents) {
            vbprintf("keeping event %3i:%8i %cPol  localized to (%6.2f,%6.2f)  triggertime %10i:%09i \n", 
                    runNumber, header->eventNumber, pChars[p], mainPeak.latitude, mainPeak.longitude, header->triggerTime, header->triggerTimeNs);
            vbprintf("   cal pulse differential is %i \n", nsDiff);
          }
        }
        // calculate cut redundancies
        for (int j=0; j<NUM_A_CUTS-1; ++j) {
          for (int k=j+1; k<NUM_A_CUTS; ++k) {
            if (aCuts[p][k]>0 && aCuts[p][j]>0 && keep[p]) {
              ++redundancyCounts[p][j][k];
            }
          }
        }         
        vbprintf(" incrementing event count %i \n", p);
        eventCounts[p] += eventWeight[p];
        if (simulationMode && cohSnr2 > 5.0) { eventCounts5[p] += eventWeight[p]; }
        vbprintf(" incremented event count %i to %f \n", p, eventCounts[p]);
      }
      //vbprintf ("110 current directory is %s \n", gDirectory->GetName());  
      if (keep[0] || keep[1]) {
        int ec, nc;
        double ea, no;
        vcout << "  keeping event after quality cuts" << endl;
        AnitaEventSummary::WaveformInfo wfi[2]; for (int p=0; p<2; ++p) {wfi[p] = eventSummary->coherent[p][maxPeakIndex[p][0]];}
        float linPol[2] = {0};
        bool firstPol = true;
        map<int, float> hpWeights[2] = {map<int, float>()};
        for (int p=0; p<2; ++p) {
          linPol[p] = circPol  ?  sqrt(wfi[p].V*wfi[p].V + wfi[p].U*wfi[p].U)/wfi[p].I  :  sqrt(wfi[p].U*wfi[p].U + wfi[p].Q*wfi[p].Q)/wfi[p].I;
          if (keep[p]) {
            vcout << "  keeping " << pChars[p] << " event after cuts" << endl;
            passingEvents[p].push_back(eventSummary->eventNumber);
            passingEventRuns[p].push_back(eventSummary->run);
            //AnitaEventSummary::PointingHypothesis thisPeak = eventSummary->peak[p][maxPeakIndex[p][0]];
            AnitaEventSummary::PointingHypothesis thisPeak = eventSummary->peak[p][0];
            skyMapHist[p]->Fill(fmod((thisPeak.phi-gps->heading+360),360), -thisPeak.theta);
            skyMapPlHist[p]->Fill(thisPeak.phi, -thisPeak.theta);
            hwAngleHist[p]->Fill(abs(thisPeak.hwAngle));
            corrPeakHist[p]->Fill(thisPeak.value);
            hilbPeakHist[p]->Fill(wfi[p].peakHilbert);
            corrHilbHist[p]->Fill(wfi[p].peakHilbert, thisPeak.value);
            corrSnrHist[p]->Fill(wfi[p].snr2, thisPeak.value);
            peakRatioHist[p]->Fill(peakRatio[p]);
            cm3Hist[p]->Fill(wfi[p].peakMoments[2]);
            if (aCuts[p][A_CUT_CONTINENT]==0) {
              //bedmap->SurfaceLonLattoEN(thisPeak.longitude, thisPeak.latitude, ec, nc);
              bedmap->SurfaceLonLattoEN(newLon[p], newLat[p], ec, nc);
              bedmap->ENtoEaNo(ec, nc, ea, no);
              //cout << "filling localization map at " << ea << "," << no << endl;
              eventHist[p]->Fill(ea/1000, no/1000);
            }
            
            // WATCH THOSE SIGNS!!!
            float ptDiffPhi = thisPeak.phi - calSource.phi;
            while (ptDiffPhi > 180) ptDiffPhi -= 360;
            while (ptDiffPhi <= -180) ptDiffPhi += 360;
            float ptDiffTheta = -(thisPeak.theta-calSource.theta);
            float sunDiffPhi = thisPeak.phi - eventSummary->sun.phi;
            while (sunDiffPhi > 180) sunDiffPhi -= 360;
            while (sunDiffPhi <= -180) sunDiffPhi += 360;
            float sunDiffTheta = -(thisPeak.theta - eventSummary->sun.theta);

            // WATCH THOSE SIGNS!!!
            float reflTheta = reflectionAngle(gps->altitude/1000.0, -thisPeak.theta);
            float reflPhi = thisPeak.phi;
            float sunReflDiffPhi = eventSummary->sun.phi - reflPhi;        
            while (sunReflDiffPhi > 180) sunReflDiffPhi -= 360;
            while (sunReflDiffPhi <= -180) sunReflDiffPhi += 360;
            float sunReflDiffTheta = reflTheta + eventSummary->sun.theta;

            pointingHist[p]->Fill(ptDiffPhi, ptDiffTheta);                
            pointingHistF[p]->Fill(ptDiffPhi, ptDiffTheta);   
            timeHist[p]->Fill(netTrigTime);
            ptTimeHistPhi[p]->SetPoint(ptTimeHistPhi[p]->GetN(), netTrigTime, ptDiffPhi);
            //if (netTrigTime < 50) {
              //printf(" event %i, net trig time is %i, phi error is %f \n", header->eventNumber, netTrigTime, ptDiffPhi);
            //}
            //timeNsHist[p]->Fill(header->triggerTimeNs - calPulseDist[p]/1000.0);
            timeNsHist[p]->Fill(nsDiff);
            calPulseDistHist[p]->Fill(netTrigTime, calPulseDist[p]/1000.0);
            ptHistPhi[p]->Fill(ptDiffPhi);
            ptHistTheta[p]->Fill(ptDiffTheta);
            
            if (useTrigNsWindow && p==0) {
              if (abs(ptDiffPhi) > misRecDistPhi || abs(ptDiffTheta) > misRecDistTheta) {
                misRecs[p].push_back(header->eventNumber);
                lprintf("  calibration pulse misreconstruction: run %3i, event %8i %cPol\n", runNumber, header->eventNumber, pChars[p]);
                lprintf("   correct position is (%.2f,%.2f); peak position is (%.2f,%.2f); difference is (%.2f,%.2f) \n",
                        calSource.phi, calSource.theta, thisPeak.phi, thisPeak.theta, ptDiffPhi, ptDiffTheta);
              }
            }             
            //linPolHist[p]->Fill(linPol[p]);
            //linPolHist[p]->Fill(linPolFrac[p]);
            anyPolHist[p]->Fill(anyPolFrac[p]);
            snrHist0[p]->Fill(wfi[p].snr);
            snrHist1[p]->Fill(wfi[p].snr1);
            snrHist2[p]->Fill(wfi[p].snr2);
            snrHist3[p]->Fill(wfi[p].snr3);
                        //snrHist[p]->Fill(maxWfSnr[p]);
            linCircPolHist[p]->Fill(linPolFrac[p], circPolFrac[p]);
            linPolHist[p]->Fill(linPolFrac[p]);
            if (firstPol) {
              //headingHist->Fill(netTrigTime, gps->heading);
              headingHist->SetPoint(headingHist->GetN(), netTrigTime, gps->heading);
            }
            float errDist = sqrt(ptDiffPhi*ptDiffPhi + ptDiffTheta*ptDiffTheta);
            //errSnrHist[p]->Fill(wfi[p].snr, errDist);
            errSnrHist[p]->Fill(maxWfSnr[p], errDist);
            //sunHist[p]->Fill(sunDiffPhi, sunDiffTheta);
            sunHist[p]->Fill(sunDiffPhi, sunDiffTheta);
            sunReflHist[p]->Fill(sunReflDiffPhi, sunReflDiffTheta);
            corrDistHist[p]->Fill(calPulseDist[p]/1000.0, thisPeak.value, 1);
            hilbDistHist[p]->Fill(calPulseDist[p]/1000.0, wfi[p].peakHilbert, 1);
            linPolDistHist[p]->Fill(calPulseDist[p]/1000.0, linPol[p], 1);
            distPFFrac_0[p]->Fill(calPulseDist[p]/1000.0, filterPowerFrac_0[p], eventWeight[p]);
            distPFFrac_1[p]->Fill(calPulseDist[p]/1000.0, filterPowerFrac_1[p], eventWeight[p]);
            linPolPlNPhiHist[p]->Fill(fmod((thisPeak.phi-gps->heading+360),360), linPol[p], eventWeight[p]);
            linPolPlPhiHist[p]->Fill(thisPeak.phi, linPol[p], eventWeight[p]);
            peakTimeHist[p]->Fill(wfi[p].peakTime, eventWeight[p]);
                        
            //if (firstPol) { 
              //stokesIHist->Fill(stokesI, eventWeight[p]);
              //stokesQHist->Fill(stokesQ/stokesI, eventWeight[p]);
              //stokesUHist->Fill(stokesU/stokesI, eventWeight[p]);
              //stokesVHist->Fill(stokesV/stokesI, eventWeight[p]);
            //}
            // count healpix bins
            for (int k=0; k<2; ++k) {oCPolPeak[k] = cPolPeak[k];}
            //for (int k=0; k<NUM_Q_CUTS; ++k) {oQCuts[k] = qCuts[p][k];}
            for (int k=0; k<NUM_A_CUTS; ++k) {oACuts[k] = aCuts[p][k];}
            //for (int k=0; k<4; ++k) {oEllParms[k] = errEllipseParms1[p][k];}
            //printf("writing event entry %i %cPol,  hPeak=%.1f, cPeak=%.4f, snr=%.2f \n", 
            //        oEventNumber, oPol, oHPeak, oCPeak, oCohSnr);
            if (aCuts[p][A_CUT_CONTINENT] ==0) {
              // calculate the weighting functions for each pixel
              hpEvents[p]+=1;
              //printf("counting healpix bins event %i, %i \n", header->eventNumber, hpEvents[p]);
              hpWeights[p] = map<int, float>();
              for (int k=0; k<errHexPixels[p][0].size(); ++k) {
                double thisLat, thisLon;
                bedmap->EaNoToLonLat(errHexPixels[p][0][k], errHexPixels[p][1][k], thisLon, thisLat);
                //vbprintf(" pixel ea/no = %f,%f \n", errHexPixels[p][0][k], errHexPixels[p][1][k]);
                //vbprintf(" pixel lat/lon = %f,%f \n", thisLat, thisLon);

                //include hp offsets??
                // ToDo:: Make this pol dependent?  
                double hpPhiOffset = 0; double hpThOffset = 0;
                //Oindree's new offsets
                //if ( whichPol == 'V' ) {
                  hpPhiOffset = 0.56;
                  hpThOffset = -5.04;
                //}
                //if ( whichPol == 'H' ) {
                //  hpPhiOffset = 3.92;
                //  hpThOffset = -0.00;
                //}

                double thisTheta = (-thisLat+90.0+hpThOffset) * M_PI/180.0;
                double thisPhi = (thisLon+hpPhiOffset) * M_PI/180.0;                

                while (thisPhi < -M_PI) thisPhi += M_PI;
                //vbprintf("  pixel theta/phi = %f,%f \n", thisTheta, thisPhi);
                pointing point(thisTheta, thisPhi);
                int pixelNo = healpix->ang2pix(point);

                vbprintf("   Healpix index is %i \n", pixelNo);
                map<int, float>::iterator thisEntry = hpEventCounts[p].find(pixelNo);
                //float thisSize  = errHexPixels[p][0].size();
                //float thisWeight = 1.0/thisSize;
                //float thisWeight = errHexPixels[p][2][k];    maybe a bug
                float thisWeight = errHexPixels[p][2][k] * eventWeight[p];
                if (thisEntry == hpEventCounts[p].end()) {
                  vbprintf(" creating new entry for pixel %i    weight=%f \n", pixelNo, thisWeight);
                  hpEventCounts[p].insert(pair<int, float>(pixelNo, thisWeight));
                } else {
                  vbprintf(" updating entry for pixel %i    weight=%f    new value=%f \n", pixelNo, thisWeight, hpEventCounts[p][pixelNo]);
                  hpEventCounts[p][pixelNo] += thisWeight;
                }
                thisEntry = hpWeights[p].find(pixelNo);
                if (thisEntry == hpWeights[p].end()) {
                  vbprintf("  creating bin list entry for event %i %cPol bin %i  %8.3f \n", header->eventNumber, pChars[p], pixelNo, thisWeight);
                  hpWeights[p].insert(pair<int, float>(pixelNo, thisWeight));
                } else {
                  vbprintf("  updating bin list entry for event %i %cPol bin %i  %8.3f \n", 
                        header->eventNumber, pChars[p], pixelNo, thisWeight);
                  hpWeights[p][pixelNo] += thisWeight;
                }
              }
              for (pair<int, float> thisEntry : hpWeights[p]) {
                oRunNumber = runNumber;
                oHpBin = thisEntry.first;
                oHpWeight = thisEntry.second;
                oEventNumber = header->eventNumber;
                oPol = pChars[p];
                vbprintf(" writing bin list entry for event %i %cPol bin %i  %8.3f \n", 
                        oEventNumber, oPol, oHpBin, oHpWeight);
                outputTree1->Fill();
              }
            }
          } else {
            for (int k=0; k<3; ++k) {errEllipseParms0[p][k] = 0;}
          }
          firstPol = false;
        }
        // populate circPol maps
        //vbprintf ("120 current directory is %s \n", gDirectory->GetName());  
        for (int p=0; p<2; ++p) {
          // careful now airstream driver: p means R/Lpol here
          //AnitaEventSummary::PointingHypothesis thisPeakC = cPolEventSummary->peak[p][maxPeakIndex[p][0]];            
          AnitaEventSummary::PointingHypothesis thisPeakC = cPolEventSummary->peak[p][0];            
          skyMapCircHist[p]->Fill(fmod((thisPeakC.phi-gps->heading+360),360), -thisPeakC.theta);
          skyMapPlCircHist[p]->Fill(thisPeakC.phi, -thisPeakC.theta);
          // careful now airstream driver: p means H/Vpol here
          //printf(" pol %i snr=%f, lesserPeak=%f, weight=%f \n", p, wfi[p].snr, cPolLesserPeak, eventWeight[p]);
          lesserCircPeakSNRHist[p]->Fill(wfi[p].snr2, cPolLesserPeak, eventWeight[p]); 
          circPeakSepSNRHist[p]->Fill(wfi[p].snr2, cPolMinDist, eventWeight[p]);
          lesserCircPeakLDHist[p]->Fill(ldCutVal[p], cPolLesserPeak, eventWeight[p]); 
          circPeakSepLDHist[p]->Fill(ldCutVal[p], cPolMinDist, eventWeight[p]);
          
          for (int lp=0; lp<2; ++lp) {
            // find the minimum correlation value corresponding to the appropriate linPol
            // p = circpol index; lp = linpol
            //printf(" maxPwVal[0][%i] = %f \n", lp, maxPwVal[0][lp]);
            //printf(" maxPwVal[1][%i] = %f \n", lp, maxPwVal[1][lp]);
            minCircPeakVal[lp] = min(maxPwVal[0][lp], maxPwVal[1][lp]) / eventSummary->peak[lp][0].value;
            //printf(" minval = %f \n", minCircPeakVal[lp]);
            minCircValHist[lp]->Fill(minCircPeakVal[lp]);
          }
          
        }

        float maxEventWeight = max(eventWeight[0], eventWeight[1]);
        //float maxSnr = max(wfi[0].snr, wfi[1].snr);
        peakSepCircHist->Fill(cPolMinDist, maxEventWeight);
        lesserPeakCircHist->Fill(cPolLesserPeak, maxEventWeight);
        
        circPeakHist->Fill(cPolPeak[0], cPolPeak[1], maxEventWeight);
        cPolPeakDistHist->Fill(cPolLesserPeak, cPolMinDist, maxEventWeight);
        // get the hpBin count map entries for the two pols

        //printf("hpWeights[0] has %i entries \n", hpWeights[0].size()); 
        //printf("hpWeights[1] has %i entries \n", hpWeights[1].size()); 
        //float thisWeight0 = 0;
        //float thisWeight1 = 0;
        /*
        map<int, float>::iterator thisEntry0 = hpWeights[0].find(ourHpBin);
        map<int, float>::iterator thisEntry1 = hpWeights[1].find(ourHpBin);
        if (thisEntry0 != hpWeights[0].end()) {
          //printf("found bin %i 0Pol \n", ourHpBin);
          thisWeight0 = hpWeights[0][ourHpBin];
          lprintf("found bin %i 0Pol wgt=%f \n", ourHpBin, thisWeight0);
        }
        if (thisEntry1 != hpWeights[1].end()) {
          thisWeight1 = hpWeights[1][ourHpBin];
          lprintf("found bin %i 1Pol wgt=%f \n", ourHpBin, thisWeight1);
        }
        */
        if (maxEventWeight > 0) {
          float minPeakVal = min(eventSummary->peak[0][0].value, eventSummary->peak[1][0].value);
          //float maxPeakVal = max(eventSummary->peak[0][0].value, eventSummary->peak[1][0].value);
          //float peakNorm = 1.414 / sqrt(minPeakVal*minPeakVal + maxPeakVal*maxPeakVal);
          corrPeakPolHist->Fill(eventSummary->peak[0][0].value, eventSummary->peak[1][0].value, maxEventWeight);
          corrPeakLesserPolHist->Fill(minPeakVal, maxEventWeight);
        }
        //vbprintf ("130 current directory is %s \n", gDirectory->GetName());
        // output to root file.  
        for (int p=0; p<2; ++p) {
          if (keep[p]) {
            oRunNumber = runNumber;
            oEventNumber = header->eventNumber;
            oPol = pChars[p];
            oEventWeight = eventWeight[p];
            oHPeak = wfi[p].peakHilbert;
            oCPeak = eventSummary->peak[p][0].value;
            oCohSnr = wfi[p].snr;
            oCohSnr1 = wfi[p].snr1;
            oCohSnr2 = wfi[p].snr2;
            oCohSnr3 = wfi[p].snr3;
            oMaxSnr = maxWfSnr[p];
            oPfFrac_0 = filterPowerFrac_0[p];
            oPfFrac_1 = filterPowerFrac_1[p];
            oLinPolFrac = linPolFrac[p];
            oDeadTimeFrac = deadTimeFrac[p];
            oPeakRatio = peakRatio[p];
            oLdCutVal = ldCutVal[p];
            oCalPulseDist = calPulseDist[p]/1000;
            oCPolDist = cPolMinDist;
            oCPolLesserPeak = cPolLesserPeak;
            oMinCircPeakVal = minCircPeakVal[p];

            //Jacob!
            //for (int k=0; k<4; ++k) {oEllParms[k] = errEllipseParms1[p][k];}
            oTheta = eventSummary->peak[p][0].theta;
            oPhi = eventSummary->peak[p][0].phi;
            oEa = srcEa[p];
            oNo = srcNo[p]; 
            oLat = newLat[p];
            oLon = newLon[p];
            oAlt = newAlt[p];
            oGps = gps;  
            oGpsRaw = gpsRaw;

            oLTheta = cPolTheta[0]; 
            oRTheta = cPolTheta[1];
            oLPhi = cPolPhi[0];
            oRPhi = cPolPhi[1];
            oLSnr2 = cPolSnr2[0];
            oRSnr2 = cPolSnr2[1];

            oTruthPhi = truthPhi;
            oTruthTheta = truthTheta;

            //printf("ophi = %f,otheta = %f,olat = %f,olon = %f,oalt = %f \n", oTheta, oPhi, oLat, oLon, oAlt);
            //sleep(10);

            outputTree0->Fill();          
          }
        }

        ++numEvents;
        if (numEvents % 100 == 0) lprintf("%i events processed \n", numEvents);
      }
      delete gps; gps=0;
    } // else ++e;
    
  }
  // TODO change this to a giant string that can then be output wherever we want
  

  lprintf("%i events processed \n", numEvents);
  lprintf("%i peak order mismatches ocurred \n", swapCount);
  lprintf("interferometry peak failures (%c,%c)  (%i,%i) \n", pChars[0], pChars[1], peakFailCounts[0], peakFailCounts[1]);




  char cutTableStr[16384];
  //lprintf ("140 current directory is %s \n", gDirectory->GetName());  

  sprintf (cutTableStr, "\nPre Quality cuts: \n");
  for (int p=1; p<2; ++p) {
    sprintf(cutTableStr+strlen(cutTableStr), " both polarizations \n");
    sprintf(cutTableStr+strlen(cutTableStr), "   total events processed       %8.3f\n", eventCounts[p]);
    sprintf(cutTableStr+strlen(cutTableStr), "  cut-id       description              as first cut        as ordered cut        as last cut \n");
    sprintf(cutTableStr+strlen(cutTableStr), "                                      number  fraction     number  fraction     number  fraction \n");
    for (int k=0; k<preQCutOrder.size(); ++k) {
      //  lprintf("  %2i %28s %8i    fraction %8.5f               if last cut:  %8.3f   fraction %8.5f\n",
      sprintf(cutTableStr+strlen(cutTableStr), "  %2i %28s   %8.3f  %8.5f   %8.3f  %8.5f   %8.3f  %8.5f \n",
              preQCutOrder[k], PRE_Q_CUT_DESC[preQCutOrder[k]], preQualityCutCounts[p][preQCutOrder[k]], (float)preQualityCutCounts[p][preQCutOrder[k]]/eventCounts[p],
              preQualityCutOrderedCounts[p][preQCutOrder[k]], (float)preQualityCutOrderedCounts[p][preQCutOrder[k]]/eventCounts[p],
              preQualityCutUniqueCounts[p][preQCutOrder[k]], (float)preQualityCutUniqueCounts[p][preQCutOrder[k]]/eventCounts[p]);
    }
    sprintf(cutTableStr+strlen(cutTableStr), "  total events cut:                                      %8.3f  %8.5f \n",
            preQualityTotalOrderedCount[p], (float)preQualityTotalOrderedCount[p]/eventCounts[p]);
    sprintf(cutTableStr+strlen(cutTableStr), "  surviving events:                                      %8.3f  %8.5f \n",
            eventCounts[p]-preQualityTotalOrderedCount[p], (float)(eventCounts[p]-preQualityTotalOrderedCount[p])/eventCounts[p]);
  }

  if(simulationMode) {
    sprintf (cutTableStr+strlen(cutTableStr), "\nPre Quality cuts (for events with > 5 SNR): \n");
    for (int p=1; p<2; ++p) {
      sprintf(cutTableStr+strlen(cutTableStr), " both polarizations \n");
      sprintf(cutTableStr+strlen(cutTableStr), "   total events processed       %8.3f\n", eventCounts5[p]);
      sprintf(cutTableStr+strlen(cutTableStr), "  cut-id       description              as first cut        as ordered cut        as last cut \n");
      sprintf(cutTableStr+strlen(cutTableStr), "                                      number  fraction     number  fraction     number  fraction \n");
      for (int k=0; k<preQCutOrder.size(); ++k) {
        //  lprintf("  %2i %28s %8i    fraction %8.5f               if last cut:  %.38i   fraction %8.5f\n",
        sprintf(cutTableStr+strlen(cutTableStr), "  %2i %28s   %8.3f  %8.5f   %8.3f  %8.5f   %8.3f  %8.5f \n",
                preQCutOrder[k], PRE_Q_CUT_DESC[preQCutOrder[k]], preQualityCutCounts5[p][preQCutOrder[k]], (float)preQualityCutCounts5[p][preQCutOrder[k]]/eventCounts5[p],
                preQualityCutOrderedCounts5[p][preQCutOrder[k]], (float)preQualityCutOrderedCounts5[p][preQCutOrder[k]]/eventCounts5[p],
                preQualityCutUniqueCounts5[p][preQCutOrder[k]], (float)preQualityCutUniqueCounts5[p][preQCutOrder[k]]/eventCounts5[p]);
      }
      sprintf(cutTableStr+strlen(cutTableStr), "  total events cut:                                      %8.3f  %8.5f \n",
              preQualityTotalOrderedCount5[p], (float)preQualityTotalOrderedCount5[p]/eventCounts5[p]);
      sprintf(cutTableStr+strlen(cutTableStr), "  surviving events:                                      %8.3f  %8.5f \n",
              eventCounts5[p]-preQualityTotalOrderedCount5[p], (float)(eventCounts5[p]-preQualityTotalOrderedCount5[p])/eventCounts5[p]);
    }
  }


  sprintf (cutTableStr+strlen(cutTableStr), "\nQuality cuts: \n");
  for (int p=0; p<2; ++p) {
    sprintf(cutTableStr+strlen(cutTableStr), " polarization %c\n", pChars[p]);
    sprintf(cutTableStr+strlen(cutTableStr), "   total events processed       %8.3f\n", preQEventCounts[p]);
    sprintf(cutTableStr+strlen(cutTableStr), "  cut-id       description              as first cut        as ordered cut        as last cut \n");
    sprintf(cutTableStr+strlen(cutTableStr), "                                      number  fraction     number  fraction     number  fraction \n");
    for (int k=0; k<qCutOrder.size(); ++k) {
      //  lprintf("  %2i %28s %8i    fraction %8.5f               if last cut:  %8.3f   fraction %8.5f\n",
      sprintf(cutTableStr+strlen(cutTableStr), "  %2i %28s   %8.3f  %8.5f   %8.3f  %8.5f   %8.3f  %8.5f \n",
              qCutOrder[k], Q_CUT_DESC[qCutOrder[k]], qualityCutCounts[p][qCutOrder[k]], (float)qualityCutCounts[p][qCutOrder[k]]/preQEventCounts[p],
              qualityCutOrderedCounts[p][qCutOrder[k]], (float)qualityCutOrderedCounts[p][qCutOrder[k]]/preQEventCounts[p],    
              qualityCutUniqueCounts[p][qCutOrder[k]], (float)qualityCutUniqueCounts[p][qCutOrder[k]]/preQEventCounts[p]);      
    }
    sprintf(cutTableStr+strlen(cutTableStr), "  total events cut:                                      %8.3f  %8.5f \n", 
            qualityTotalOrderedCount[p], (float)qualityTotalOrderedCount[p]/preQEventCounts[p]);
    sprintf(cutTableStr+strlen(cutTableStr), "  surviving events:                                      %8.3f  %8.5f \n", 
            preQEventCounts[p]-qualityTotalOrderedCount[p], (float)(preQEventCounts[p]-qualityTotalOrderedCount[p])/preQEventCounts[p]);
  }
 
  // for sim, also show cut efficenct for >5 SNR
  if(simulationMode) {
    sprintf (cutTableStr+strlen(cutTableStr), "\nQuality cuts (for events with > 5 SNR): \n");
    for (int p=0; p<2; ++p) {
      sprintf(cutTableStr+strlen(cutTableStr), " polarization %c\n", pChars[p]);
      sprintf(cutTableStr+strlen(cutTableStr), "   total events processed       %8.3f\n", preQEventCounts5[p]);
      sprintf(cutTableStr+strlen(cutTableStr), "  cut-id       description              as first cut        as ordered cut        as last cut \n");
      sprintf(cutTableStr+strlen(cutTableStr), "                                      number  fraction     number  fraction     number  fraction \n");
      for (int k=0; k<qCutOrder.size(); ++k) {
        //  lprintf("  %2i %28s %8i    fraction %8.5f               if last cut:  %.38i   fraction %8.5f\n",
        sprintf(cutTableStr+strlen(cutTableStr), "  %2i %28s   %8.3f  %8.5f   %8.3f  %8.5f   %8.3f  %8.5f \n",
                qCutOrder[k], Q_CUT_DESC[qCutOrder[k]], qualityCutCounts5[p][qCutOrder[k]], (float)qualityCutCounts5[p][qCutOrder[k]]/preQEventCounts5[p],
                qualityCutOrderedCounts5[p][qCutOrder[k]], (float)qualityCutOrderedCounts5[p][qCutOrder[k]]/preQEventCounts5[p],
                qualityCutUniqueCounts5[p][qCutOrder[k]], (float)qualityCutUniqueCounts5[p][qCutOrder[k]]/preQEventCounts5[p]);
      }
      sprintf(cutTableStr+strlen(cutTableStr), "  total events cut:                                      %8.3f  %8.5f \n",
              qualityTotalOrderedCount5[p], (float)qualityTotalOrderedCount5[p]/preQEventCounts5[p]);
      sprintf(cutTableStr+strlen(cutTableStr), "  surviving events:                                      %8.3f  %8.5f \n",
              preQEventCounts5[p]-qualityTotalOrderedCount5[p], (float)(preQEventCounts5[p]-qualityTotalOrderedCount5[p])/preQEventCounts5[p]);
    }
  }


  
  sprintf(cutTableStr+strlen(cutTableStr), "\nAnalysis cuts, stage 1: \n");
  for (int p=0; p<2; ++p) {  
    sprintf(cutTableStr+strlen(cutTableStr), " polarization %c\n", pChars[p]);
    sprintf(cutTableStr+strlen(cutTableStr), "   total events processed       %8.3f\n", qEventCounts[p]);
    sprintf(cutTableStr+strlen(cutTableStr), "  cut-id       description              as first cut        as ordered cut        as last cut \n");
    sprintf(cutTableStr+strlen(cutTableStr), "                                      number  fraction     number  fraction     number  fraction \n");
    for (int k=0; k<aCutOrderStage1Limit; ++k) {
      sprintf(cutTableStr+strlen(cutTableStr), "  %2i %28s   %8.3f  %8.5f   %8.3f  %8.5f   %8.3f  %8.5f \n",
              aCutOrderStage1[k], A_CUT_DESC[aCutOrderStage1[k]] ,analysisCutCounts[p][aCutOrderStage1[k]], (float)analysisCutCounts[p][aCutOrderStage1[k]]/qEventCounts[p],
              analysisCutOrderedCounts[p][aCutOrderStage1[k]], (float)analysisCutOrderedCounts[p][aCutOrderStage1[k]]/qEventCounts[p],    
              analysisCutUniqueCounts[p][aCutOrderStage1[k]], (float)analysisCutUniqueCounts[p][aCutOrderStage1[k]]/qEventCounts[p]);      
    }
    sprintf(cutTableStr+strlen(cutTableStr), "  total events cut:                                      %8.3f  %8.5f \n", 
            analysisTotalOrderedCount[p], (float)analysisTotalOrderedCount[p]/qEventCounts[p]);
    sprintf(cutTableStr+strlen(cutTableStr), "  surviving events:                                      %8.3f  %8.5f \n", 
            qEventCounts[p]-analysisTotalOrderedCount[p], (float)(qEventCounts[p]-analysisTotalOrderedCount[p])/qEventCounts[p]);    
  }

  if(simulationMode) {
    sprintf(cutTableStr+strlen(cutTableStr), "\nAnalysis cuts, stage 1 (SNR > 5): \n");
    for (int p=0; p<2; ++p) {
      sprintf(cutTableStr+strlen(cutTableStr), " polarization %c\n", pChars[p]);
      sprintf(cutTableStr+strlen(cutTableStr), "   total events processed       %8.3f\n", qEventCounts5[p]);
      sprintf(cutTableStr+strlen(cutTableStr), "  cut-id       description              as first cut        as ordered cut        as last cut \n");
      sprintf(cutTableStr+strlen(cutTableStr), "                                      number  fraction     number  fraction     number  fraction \n");
      for (int k=0; k<aCutOrderStage1Limit; ++k) {
        sprintf(cutTableStr+strlen(cutTableStr), "  %2i %28s   %8.3f  %8.5f   %8.3f  %8.5f   %8.3f  %8.5f \n",
                aCutOrderStage1[k], A_CUT_DESC[aCutOrderStage1[k]] ,
                analysisCutCounts5[p][aCutOrderStage1[k]], (float)analysisCutCounts5[p][aCutOrderStage1[k]]/qEventCounts5[p],
                analysisCutOrderedCounts5[p][aCutOrderStage1[k]], (float)analysisCutOrderedCounts5[p][aCutOrderStage1[k]]/qEventCounts5[p],
                analysisCutUniqueCounts5[p][aCutOrderStage1[k]], (float)analysisCutUniqueCounts5[p][aCutOrderStage1[k]]/qEventCounts5[p]);
      }
      sprintf(cutTableStr+strlen(cutTableStr), "  total events cut:                                      %8.3f  %8.5f \n",
              analysisTotalOrderedCount5[p], (float)analysisTotalOrderedCount5[p]/qEventCounts5[p]);
      sprintf(cutTableStr+strlen(cutTableStr), "  surviving events:                                      %8.3f  %8.5f \n",
              qEventCounts5[p]-analysisTotalOrderedCount5[p], (float)(qEventCounts5[p]-analysisTotalOrderedCount5[p])/qEventCounts5[p]);
    }
  }


  lprintf(cutTableStr);

  if (saveOutput) {
    char filepath[1024];
    sprintf(filepath, "%s/cutTable.txt", outputDir);
    FILE* cutTableFile = fopen(filepath, "w");
    fprintf(cutTableFile, cutTableStr);
    fclose(cutTableFile);
  }
  
  
  if (saveOutput) {
    // TODO get rid of this TEX bullshit
    char filepath[1024];
    sprintf(filepath, "%s/cutTableQ_tex.txt", outputDir);
    FILE* cutTableFile = fopen(filepath, "w");
    fprintf (cutTableFile, "\\begin{footnotesize}\n");
    fprintf (cutTableFile, "\nQuality cuts: \\\\ \n");
    for (int p=0; p<2; ++p) {  
      fprintf(cutTableFile, "polarization %c \\\\ \n", pChars[p]);
      fprintf(cutTableFile, "total events processed       %8.3f \\\\ \n", eventCounts[p]);
      fprintf(cutTableFile, "\\begin{tabular}{ l l c c c c c c } \n");
      fprintf(cutTableFile, "cut-id & description & \\multicolumn{2}{c}{as first cut} & \\multicolumn{2}{c}{as ordered cut} & \\multicolumn{2}{c}{as last cut} \\\\ \n");
      fprintf(cutTableFile, "         &                    &       number & fraction   &  number & fraction  &   number & fraction \\\\ \n");
      for (int k=0; k<qCutOrder.size(); ++k) {
        //  lprintf("  %2i %28s %8.3f    fraction %8.5f               if last cut:  %8.3f   fraction %8.5f\n",
        fprintf(cutTableFile, "  %2i & %28s &  %8.3f & %8.5f &  %8.3f & %8.5f &  %8.3f & %8.5f \\\\ \n",
                qCutOrder[k], Q_CUT_DESC[qCutOrder[k]], qualityCutCounts[p][qCutOrder[k]], (float)qualityCutCounts[p][qCutOrder[k]]/eventCounts[p],
                qualityCutOrderedCounts[p][qCutOrder[k]], (float)qualityCutOrderedCounts[p][qCutOrder[k]]/eventCounts[p],    
                qualityCutUniqueCounts[p][qCutOrder[k]], (float)qualityCutUniqueCounts[p][qCutOrder[k]]/eventCounts[p]);      
      }
      fprintf(cutTableFile, "  & total events cut:    & & &             %8.3f & %8.5f \\\\ \n", 
              qualityTotalOrderedCount[p], (float)qualityTotalOrderedCount[p]/eventCounts[p]);
      fprintf(cutTableFile, "  & surviving events:      & & &                 %8.3f & %8.5f \\\\ \n \\\\ \n", 
              eventCounts[p]-qualityTotalOrderedCount[p], (float)(eventCounts[p]-qualityTotalOrderedCount[p])/eventCounts[p]);
      fprintf(cutTableFile, "\\end{tabular} \\par \n");
    }
    fprintf (cutTableFile, "\\end{footnotesize}\n");
    fprintf(cutTableFile, "\n");
    fclose(cutTableFile);
    sprintf(filepath, "%s/cutTableA_tex.txt", outputDir);
    cutTableFile = fopen(filepath, "w");
    
    fprintf (cutTableFile, "\\begin{footnotesize}\n");
    fprintf (cutTableFile, "\nAnalysis cuts (set 1): \\\\ \n");
    for (int p=0; p<2; ++p) {  
      fprintf(cutTableFile, "polarization %c \\\\ \n", pChars[p]);
      fprintf(cutTableFile, "total events processed       %8.3f \\\\ \n", qEventCounts[p]);
      fprintf(cutTableFile, "\\begin{tabular}{ l l c c c c c c } \n");
      fprintf(cutTableFile, "cut-id & description & \\multicolumn{2}{c}{as first cut} & \\multicolumn{2}{c}{as ordered cut} & \\multicolumn{2}{c}{as last cut} \\\\ \n");
      fprintf(cutTableFile, "       &                     &      number & fraction  &   number & fraction  &   number & fraction \\\\ \n");
      //for (int k=0; k<aCutOrderStage1.size(); ++k) {
      for (int k=0; k<aCutOrderStage1Limit; ++k) {
        fprintf(cutTableFile, "  %2i & %28s &  %8.3f & %8.5f &  %8.3f & %8.5f &  %8.3f & %8.5f \\\\ \n",
                aCutOrderStage1[k], A_CUT_DESC[aCutOrderStage1[k]] ,analysisCutCounts[p][aCutOrderStage1[k]], (float)analysisCutCounts[p][aCutOrderStage1[k]]/qEventCounts[p],
                analysisCutOrderedCounts[p][aCutOrderStage1[k]], (float)analysisCutOrderedCounts[p][aCutOrderStage1[k]]/qEventCounts[p],    
                analysisCutUniqueCounts[p][aCutOrderStage1[k]], (float)analysisCutUniqueCounts[p][aCutOrderStage1[k]]/qEventCounts[p]);      
      }
      fprintf(cutTableFile, "  & total events cut:               & & &                      %8.3f & %8.5f \\\\ \n", 
              analysisTotalOrderedCount[p], (float)analysisTotalOrderedCount[p]/qEventCounts[p]);
      fprintf(cutTableFile, "  & surviving events:               & & &                      %8.3f & %8.5f \\\\ \n \\\\ \n", 
              qEventCounts[p]-analysisTotalOrderedCount[p], (float)(qEventCounts[p]-analysisTotalOrderedCount[p])/qEventCounts[p]);    
      fprintf(cutTableFile, "\\end{tabular} \\par \n");

    }
    fprintf(cutTableFile, "\\end{footnotesize}\n");
    fclose(cutTableFile);
  }  
  
  cout << "\nCut redundancy table" << endl;
  for (int p=0; p<2; ++p) {
    lprintf(" polarization %c\n",pChars[p]);
    lprintf("%34s ", ""); for (int k=0; k<aCutOrderStage1Limit; ++k) {printf("%2i ", aCutOrderStage1[k]);} lprintf("\n");
    lprintf("%34s ", ""); for (int k=0; k<aCutOrderStage1Limit; ++k) {printf("%8s ", "--");} lprintf("\n");    
    for (int j=0; j<aCutOrderStage1Limit; ++j) {    
      int m = aCutOrderStage1[j];
      lprintf("%8i %26s", m, A_CUT_DESC[m]); 
      for (int k=0; k<aCutOrderStage1Limit; ++k) {
        int n = aCutOrderStage1[k];
        if (k>j) {  
          lprintf("%8i ", redundancyCounts[p][m][n]);
        } else if (k<j) {
          lprintf("%8i ", redundancyCounts[p][n][m]);
        } else {
          lprintf("%8.3f ", analysisCutCounts[p][m]);
        }
      }
      lprintf("\n");
    }
  }

  lprintf("%i events with peak value order violations \n", peakViolationCount); 


  if (listPassingEvents) {
    lprintf("Events passing plotting cuts: \n");
    char filename[1024];
    sprintf(filename, "%s/eventList_%03i_%03i_%sPol.txt", outputDir, startRunNumber, endRunNumber, pChars);
    ofstream eventListOutFile(filename);
    for (int p : whichPols) {
      lprintf(" %cPol \n", pChars[p]);
      for (int k=0; k<passingEvents[p].size(); ++k) {
        lprintf(" %3i  %8i \n", passingEventRuns[p][k], passingEvents[p][k]);
        eventListOutFile << passingEvents[p][k] << " " << pChars[p] << " " << passingEventRuns[p][k]<< endl;
      }
      eventListOutFile.close();
    }
  }


  for (int p : whichPols) {
    if (useTrigNsWindow) {printf("%cPol: %lu calibration pulse misreconstructions \n", pChars[p], misRecs[p].size());}
    if (useTrigNsWindow || simulationMode) {
      for (int snrIndex = 0; snrIndex<20; ++snrIndex) {
        float thisRms = eventCountsSnr[p][snrIndex] > 0 ? sqrt(errPhiBySnr[p][snrIndex]/eventCountsSnr[p][snrIndex]) : 0;
        lprintf("  snr %i  rms=%f \n", snrIndex, thisRms);
        ptHistPhiSnr[p]->Fill(snrIndex, thisRms);
        thisRms = eventCountsSnr[p][snrIndex] > 0 ? sqrt(errThetaBySnr[p][snrIndex]/eventCountsSnr[p][snrIndex]) : 0;
        ptHistThetaSnr[p]->Fill(snrIndex, thisRms);
        lprintf(" snr=%i, keeps=%.3f, total events=%.3f \n", snrIndex, acceptsSnr[p][snrIndex], eventCountsSnr[p][snrIndex]);
        float thisAcceptance = (eventCountsSnr[p][snrIndex]>0) ? (float)acceptsSnr[p][snrIndex]/(float)eventCountsSnr[p][snrIndex] : 0;
        acceptHistSnr[p]->Fill(snrIndex, thisAcceptance);
      }
      for (int b=0; b<=acceptHistSnr[p]->GetNbinsX(); ++b) {
        if (acceptHistSnr[p]->GetBinContent(b)>0) {
          acceptHistSnr[p]->SetBinError(b, 0.01);
        }
      }
    }
  }

  // write out pixel event counts
  lprintf("Healpix bin event counts \n");
  for (int p : whichPols) {
    lprintf(" %cPol: %i events \n", pChars[p], hpEvents[p]);
    float thisCount = 0;
    for (pair<int, float> thisEntry : hpEventCounts[p]) {
      pointing thisPoint = healpix->pix2ang(thisEntry.first);
      float thisLat = 90.0 - (thisPoint.theta * 180.0/M_PI) ;
      float thisLon = thisPoint.phi * 180.0/M_PI;
      lprintf("   bin %8i: %10.3f events.      theta=%4.3f, phi=%4.3f, lat=%6.3f, lon=%6.3f\n", 
              thisEntry.first, thisEntry.second, thisPoint.theta, thisPoint.phi, thisLat, thisLon);
      thisCount += thisEntry.second;
    }
    lprintf(" total  %.3f \n", thisCount);
  }
  lprintf("filling the healpix map histogram \n");
  // populate the Healpix maps of the continent
  for (int p : whichPols) {
    int count = 0;
    lprintf("healpix map histogram has %i bins \n", hpCountsHist[p]->GetNcells());
    for (int gBin = 0; gBin <hpCountsHist[p]->GetNcells()-1; ++gBin) {
    int xBin, yBin, zBin;
    hpCountsHist[p]->GetBinXYZ(gBin, xBin, yBin, zBin);
    float ea = hpCountsHist[p]->GetXaxis()->GetBinCenter(xBin) * 1000;
    float no = hpCountsHist[p]->GetYaxis()->GetBinCenter(yBin) * 1000;
      double thisLat, thisLon;
      bedmap->EaNoToLonLat(ea, no, thisLon, thisLat);
      double thisTheta = (-thisLat+90.0) * M_PI/180.0;
      double thisPhi = thisLon * M_PI/180.0;
      //while (thisPhi < -M_PI) thisPhi += 2*M_PI;
      pointing point(thisTheta, thisPhi);
      int pixelNo = healpix->ang2pix(point);
      map<int, float>::iterator thisEntry = hpEventCounts[p].find(pixelNo);
      float eventCount = 0;
      if (thisEntry != hpEventCounts[p].end()) {
        eventCount = hpEventCounts[p][pixelNo];
      }
      //printf(" ea=%8.0f, no=%8.0f, lat=%6.3f, lon=%6.3f,  theta=%4.3f, phi=%4.3f,  pixelNo=%i, count=%8.1f \n", 
      //        ea, no, thisLat, thisLon, thisTheta, thisPhi, pixelNo, eventCount);
      hpCountsHist[p]->Fill(ea/1000, no/1000, eventCount);
      ++count;
      //if (count%100 == 0) {printf("  %8i entries processed", count);}
    }
  } 
  
  if (saveOutput) {
    outputFile->cd();
    //TFile* outputFile = new TFile(outputFilename, "RECREATE");
    lprintf ("result tree 0 has %llu entries \n", outputTree0->GetEntries());
    lprintf ("result tree 1 has %llu entries \n", outputTree1->GetEntries());
    lprintf ("150 current directory is %s \n", gDirectory->GetName());
    outputTree0->SetDirectory(outputFile);
    outputTree1->SetDirectory(outputFile);
    outputTree0->Write();
    outputTree1->Write();
    outputFile->Save();
    //outputFile->Close();
  }
    
  //lprintf ("160 current directory is %s \n", gDirectory->GetName());

  TCanvas* anyPolCanv = 0;
  TCanvas* snrCanv0 = 0;
  TCanvas* snrCanv1 = 0;
  TCanvas* snrCanv2 = 0;
  TCanvas* snrCanv3 = 0;
  TCanvas* sunCanv = 0;
  TCanvas* sunReflCanv = 0;
  TCanvas* timeCanv = 0;
  TCanvas* ptTimeCanvTh = 0;
  TCanvas* calPulseDistCanv = 0;
  TCanvas* timeNsCanv = 0;
  TCanvas* linCircPolCanv = 0;
  TCanvas* linPolCanv = 0;
  
  TCanvas* skyMapPlCanv = 0;
  TCanvas* skyMapCanv = 0;
  TCanvas* skyMapPlCircCanv = 0;
  TCanvas* skyMapCircCanv = 0;
  TCanvas* corrPeakCanv = 0;
  TCanvas* hilbPeakCanv = 0;
  TCanvas* corrHilbCanv = 0;
  TCanvas* corrSnrCanv = 0;
  TCanvas* hwAngleCanv = 0;
  TCanvas* conMapCanv = 0;
  TCanvas* conMapCanvF = 0;
  TCanvas* peakRatioCanv = 0;
  TCanvas* cm3Canv = 0;
  TCanvas* hpCountsCanv = 0;
  TCanvas* pointingCanvC = 0;
  TCanvas* pointingCanvF = 0;
  TCanvas* errSnrCanv = 0;
  TCanvas* errSnrCanv1 = 0;
  TCanvas* pointingCanv1 = 0;
  TCanvas* acceptSnrCanv = 0;
  TCanvas* distSnrEffCanv = 0;
  TCanvas* corrDistCanv = 0;
  TCanvas* hilbDistCanv = 0;
  TCanvas* linPolDistCanv = 0;
  TCanvas* effSnrPFFracCanv = 0;
  TCanvas* distPFFracCanv = 0;
  TCanvas* stokesCanv = 0;
  TCanvas* linPolPlPhiCanv = 0;
  TCanvas* headingCanv = 0;
  TCanvas* corrPeakPolCanv = 0;
  TCanvas* corrPeakLesserPolCanv = 0;
  TCanvas* peakSepCircCanv = 0;
  TCanvas* lesserPeakCircCanv = 0;
  TCanvas* circPeakCanv = 0;
  TCanvas* cPolPeakDistCanv = 0;
  TCanvas* lesserCircPeakSNRCanv = 0;
  TCanvas* circPeakSepSNRCanv = 0;
  TCanvas* lesserCircPeakLDCanv = 0;
  TCanvas* circPeakSepLDCanv = 0;
  TCanvas* minCircValCanv = 0;
  TCanvas* peakTimeCanv = 0;

  
  
                      //anyPolCanv = new TCanvas("anyPolCanv", "Total polarization fraction", 500, numPols*275);
                      //anyPolCanv->Divide(1,numPols);
    
  snrCanv0 = new TCanvas("snrCanv0", "SNR0", 500, numPols*275);
  snrCanv0->Divide(1,numPols);
  snrCanv1 = new TCanvas("snrCanv1", "SNR1", 500, numPols*275);
  snrCanv1->Divide(1,numPols);
  snrCanv2 = new TCanvas("snrCanv2", "SNR2", 500, numPols*275);
  snrCanv2->Divide(1,numPols);
  snrCanv3 = new TCanvas("snrCanv3", "SNR3", 500, numPols*275);
  snrCanv3->Divide(1,numPols);
  
  sunCanv = new TCanvas("sunCanv", "Sun coordinates", 500, numPols*275);
  sunCanv->Divide(1,numPols);
  
  
  sunReflCanv = new TCanvas("sunReflCanv", "Sun reflection coordinates", 500, numPols*275);
  sunReflCanv->Divide(1,numPols);
  
                      
                      //timeCanv = new TCanvas("timeCanv", "Time-resolved", 500, numPols*275);
                      //timeCanv->Divide(1,numPols);
                      
  
  ptTimeCanvTh = new TCanvas("ptTimeCanvTh", "Time-resolved pointing", 500, numPols*275);
  ptTimeCanvTh->Divide(1,numPols);
  
  timeNsCanv = new TCanvas("timeNsCanv", "ns trigger time", 500, numPols*275);
  timeNsCanv->Divide(1,numPols);
  
  linCircPolCanv = new TCanvas("linCircPolCanv", "Linear/circular polarization fraction", 450, numPols*450);
  linCircPolCanv->Divide(1,numPols);
  
  linPolCanv = new TCanvas("linPolCanv", "Linear polarization fraction", 450, numPols*450);
  linPolCanv->Divide(1,numPols);
  
  skyMapPlCanv = new TCanvas("skyMapPlCanv", "Reconstruction sky maps", 600, numPols*300);
  skyMapPlCanv->Divide(1,numPols);
  
  skyMapCanv = new TCanvas("skyMapCanv", "Reconstruction sky maps", 600, numPols*300);
  skyMapCanv->Divide(1,numPols);
  
  skyMapPlCircCanv = new TCanvas("skyMapPlCircCanv", "Reconstruction sky maps", 600, numPols*300);
  skyMapPlCircCanv->Divide(1,numPols);
  
  skyMapCircCanv = new TCanvas("skyMapCircCanv", "Reconstruction sky maps", 600, numPols*300);
  skyMapCircCanv->Divide(1,numPols);
  
  corrPeakCanv = new TCanvas("corrPeakCanv", "Correlation peak value", 600, numPols*300);
  corrPeakCanv->Divide(1,numPols);
  
  hilbPeakCanv = new TCanvas("hilbPeakCanv", "Hilbert peak value", 600, numPols*300);
  hilbPeakCanv->Divide(1,numPols);
  
  corrHilbCanv = new TCanvas("corrHilbCanv", "Correlation peak vs. Hilbert peak value", 450, numPols*450);
  corrHilbCanv->Divide(1,numPols);

  corrSnrCanv = new TCanvas("corrSnrCanv", "Correlation peak vs. SNR", 450, numPols*450);
  corrSnrCanv->Divide(1,numPols);
                        
                      //hwAngleCanv = new TCanvas("hwAngleCanv", "Peak-to-trigger angle", 600, numPols*300);
                      //hwAngleCanv->Divide(1,numPols);
                      
  conMapCanv = new TCanvas("conMapCanv", "Event localizations", 455, numPols*450);
  conMapCanv->Divide(1, numPols);
  
  conMapCanvF = new TCanvas("conMapCanvF", "Event localization", 900, 900);
  //conMapCanvF->Divide(1, numPols);
  
  peakRatioCanv = new TCanvas("peakRatioCanv", "Coherent sum peak ratio", 450, numPols*450);
  peakRatioCanv->Divide(1, numPols);
                      
                      //cm3Canv = new TCanvas("cm3PeakCanv", "3rd central moment", 600, numPols*300);
                      //cm3Canv->Divide(1,numPols);
                      
  hpCountsCanv = new TCanvas("hpCountsCanv", "Healpix binning counts", 555, numPols*550);
  hpCountsCanv->Divide(1,numPols);

  linPolPlPhiCanv = new TCanvas("linPolPlNPhiCanv", "Linear polarization fraction vs azimuth (pl-north)", 800, numPols*400);
  linPolPlPhiCanv->Divide(2, numPols);
  
  headingCanv = new TCanvas("headingCanv", "Payload heading vs trigger time", 600, 300);

                      //stokesCanv = new TCanvas("stokesCanv", "stokesCanv", 1200, 400);
                      //stokesCanv->Divide(3, 1);
  corrPeakPolCanv = new TCanvas("corrPeakPolCanv", "corr peak by pol", 500, 500);
  //corrPeakLesserPolCanv = new TCanvas("corrPeakLesserPolCanv", "lesser corr peak by pol", 500, 500);
  peakSepCircCanv = new TCanvas("peakSepCircCanv", "cPol peak separation", 500, 500);
  lesserPeakCircCanv = new TCanvas("lesserPeakCircCanv", "cPol lesser peak", 500, 500);
  circPeakCanv = new TCanvas("circPeakCanv", "cPol peak", 500, 500);
  cPolPeakDistCanv = new TCanvas("cPolPeakDistCanv", "cPol peak vs distance", 500, 500);
  
  lesserCircPeakSNRCanv = new TCanvas("lesserCircPeakSNRCanv", "lesserCircPeakSNRCanv", 400, numPols*400);
  lesserCircPeakSNRCanv->Divide(1, 2);
  circPeakSepSNRCanv = new TCanvas("circPeakSepSNRCanv", "circPeakSepSNRCanv", 400, numPols*400);
  circPeakSepSNRCanv->Divide(1, 2);
  lesserCircPeakLDCanv = new TCanvas("lesserCircPeakLDCanv", "lesserCircPeakLDCanv", 400, numPols*400);
  lesserCircPeakLDCanv->Divide(1, 2);
  circPeakSepLDCanv = new TCanvas("circPeakSepLDCanv", "circPeakSepLDCanv", 400, numPols*400);
  circPeakSepLDCanv->Divide(1, 2);
  minCircValCanv = new TCanvas("minCircValCanv", "minCircValCanv", 400, numPols*400);
  minCircValCanv->Divide(1, 2);

  if (useTrigNsWindow) {
    calPulseDistCanv = new TCanvas("calPulseDistCanv", "Distance from cal pulser", 500, numPols*275);
    calPulseDistCanv->Divide(1,numPols);
    pointingCanvC = new TCanvas("pointingCanvC", "Calibration pulse pointing", 700, numPols*350);
    pointingCanvC->Divide(1, numPols);
    pointingCanvF = new TCanvas("pointingCanvF", "Calibration pulse pointing", 700, numPols*350);
    pointingCanvF->Divide(1, numPols);
    pointingCanv1 = new TCanvas("pointingCanv1", "Calibration pulse pointing", numPols*700, 700);
    pointingCanv1->Divide(numPols, 2);
    errSnrCanv = new TCanvas("errSnrCanv", "Pointing error vs SNR", 350, numPols*350);
    errSnrCanv->Divide(1,numPols);
    errSnrCanv1 = new TCanvas("errSnrCanv1", "Calibration pulse pointing by SNR", numPols*350, 700);
    errSnrCanv1->Divide(numPols, 2);
    corrDistCanv = new TCanvas("corrDistCanv", "Correlation peak vs distance ", 400, numPols*400);
    corrDistCanv->Divide(1, numPols);
    hilbDistCanv = new TCanvas("hilbDistCanv", "Hilbert peak vs distance ", 400, numPols*400);
    hilbDistCanv->Divide(1, numPols);
    linPolDistCanv = new TCanvas("linPolDistCanv", "Lin pol frac vs distance ", 400, numPols*400);
    linPolDistCanv->Divide(1, numPols);
    distPFFracCanv = new TCanvas("distPFFracCanv", "Filtered power frac vs distance ", 800, numPols*400);
    distPFFracCanv->Divide(2, numPols);
  }
  if (useTrigNsWindow || simulationMode) {
    acceptSnrCanv = new TCanvas("acceptSnrCanv", "Acceptance by SNR", 350, numPols*350);
    acceptSnrCanv->Divide(1, numPols);
    distSnrEffCanv = new TCanvas("distSnrEffCanv", "Efficiency vs. SNR, Distance", 500, numPols*500);
    distSnrEffCanv->Divide(1, numPols);
    effSnrPFFracCanv = new TCanvas("effSnrPFFracCanv", "Acceptance by SNR, filter power fraction", 700, numPols*350);     
    effSnrPFFracCanv->Divide(2, numPols);     
  }
  
  vector<TMarker*> conMapMarkers(0);
  for (vector<double> thisPoint : conMapPoints) {
    int ec, nc; double ea, no;
    bedmap->SurfaceLonLattoEN(thisPoint[1], thisPoint[0], ec, nc);
    bedmap->ENtoEaNo(ec, nc, ea, no);
    bedmap->LonLattoEaNo(thisPoint[1], thisPoint[0], ea, no);
    TMarker* thisMarker = new TMarker(ea/1000, no/1000, kPlus);
    thisMarker->SetMarkerColor(kBlack);
    conMapMarkers.push_back(thisMarker);
  }

  // output cut fail event lists 
  //   for every analysis cut in the cut order
  if (false) {    
    lprintf("\nAnalysis cut event lists \n");
    for (int p=0; p<2; ++p) {
      lprintf("Polarization %c \n", pChars[p]);
      //for (int k=0; k<aCutOrderStage1.size(); ++k) {
      for (int k=0; k<aCutOrderStage1Limit; ++k) {
        int m = aCutOrderStage1[k];
        lprintf("%2i %s, %8lu events \n", m, A_CUT_DESC[m], analysisCutEventList[p][m].size());
        for (int thisEvent : analysisCutEventList[p][m]) {
          lprintf("  %i \n", thisEvent);
        }
      }
    }
  }

  // make the plots
  TLine* cutLine = 0; // Yes, I am a bad, bad man, leaker of TLine's.  But this code is only executed once (well, twice) and I don't want to keep a bunch of named variables.    
  float skyMapMax = 0;
  float skyMapMaxCirc = 0;
  float skyMapPlMax = 0;
  float skyMapPlMaxCirc = 0;
  float continentMax = 0;
  int ec, nc; double ea, no;
  bedmap->SurfaceLonLattoEN(0, -80, ec, nc);
  bedmap->ENtoEaNo(ec, nc, ea, no);  
  TEllipse* deg80Circ = new TEllipse(0, 0, no/1000);
  deg80Circ->SetLineColor(kBlack);
  deg80Circ->SetLineWidth(1);
  deg80Circ->SetLineStyle(2);
  deg80Circ->SetFillStyle(0);
  bedmap->SurfaceLonLattoEN(0, -70, ec, nc);
  bedmap->ENtoEaNo(ec, nc, ea, no);  
  TEllipse* deg70Circ = new TEllipse(0, 0, no/1000);
  deg70Circ->SetLineColor(kBlack);
  deg70Circ->SetLineWidth(1);
  deg70Circ->SetLineStyle(2);
  deg70Circ->SetFillStyle(0);
  //lprintf ("current directory is %s \n", gDirectory->GetName());
  bool firstPol = true;
  for (int w=0; w<whichPols.size(); ++w) {
    int p = whichPols[w];
    if (skyMapCanv) {
      skyMapCanv->cd(w+1);
      gStyle->SetOptStat("e");
      gPad->SetLogz();
      skyMapHist[p]->Draw("COLZ");
      int maxBinC = skyMapHist[p]->GetMaximumBin();
      float zPeakC = skyMapHist[p]->GetBinContent(maxBinC);
      skyMapMax = max(skyMapMax, zPeakC);
    }
    if (skyMapPlCanv) {
      skyMapPlCanv->cd(w+1);
      gStyle->SetOptStat("e");
      skyMapPlHist[p]->Draw("COLZ");
      int maxBinC = skyMapPlHist[p]->GetMaximumBin();
      float zPeakC = skyMapPlHist[p]->GetBinContent(maxBinC);
      skyMapPlMax = max(skyMapPlMax, zPeakC);
    }
    
    if (skyMapCircCanv) {
      skyMapCircCanv->cd(w+1);
      gStyle->SetOptStat("e");
      gPad->SetLogz();
      skyMapCircHist[p]->Draw("COLZ");
      int maxBinC = skyMapCircHist[p]->GetMaximumBin();
      float zPeakC = skyMapCircHist[p]->GetBinContent(maxBinC);
      skyMapMaxCirc = max(skyMapMaxCirc, zPeakC);
    }
    if (skyMapPlCircCanv) {
      skyMapPlCircCanv->cd(w+1);
      gStyle->SetOptStat("e");
      skyMapPlCircHist[p]->Draw("COLZ");
      int maxBinC = skyMapPlCircHist[p]->GetMaximumBin();
      float zPeakC = skyMapPlCircHist[p]->GetBinContent(maxBinC);
      skyMapPlMaxCirc = max(skyMapPlMaxCirc, zPeakC);
    }
    
    if (corrPeakCanv) {         
      corrPeakCanv->cd(w+1);
      gStyle->SetOptStat("emr");
      gPad->SetLogy();
      corrPeakHist[p]->Draw();
      corrPeakCanv->Update();
      cutLine = 0; cutLine = new TLine(corrPeakThresh, pow(10,gPad->GetUymin()), corrPeakThresh, pow(10,gPad->GetUymax()));
      cutLine->SetLineColor(kRed);
      cutLine->Draw();
    }
    if (hilbPeakCanv) {
      hilbPeakCanv->cd(w+1);
      gStyle->SetOptStat("emr");
      gPad->SetLogy();
      hilbPeakHist[p]->Draw();
      hilbPeakCanv->Update();
      cutLine = 0; cutLine = new TLine(hilbPeakThresh, pow(10,gPad->GetUymin()), hilbPeakThresh, pow(10,gPad->GetUymax()));
      cutLine->SetLineColor(kRed);
      cutLine->Draw();
    }
    if (corrHilbCanv) {
      corrHilbCanv->cd(w+1);
      gStyle->SetOptStat("emr");
      corrHilbHist[p]->SetLabelSize(0.04, "xyz");
      gPad->SetLogz();
      corrHilbHist[p]->Draw("COLZ");
      corrHilbCanv->Update();
      cutLine = 0; cutLine = new TLine(hilbPeakThresh,gPad->GetUymin(),hilbPeakThresh,gPad->GetUymax());
      cutLine->SetLineColor(kRed);
      cutLine->Draw();
      cutLine = 0; cutLine = new TLine(gPad->GetUxmin(),corrPeakThresh,gPad->GetUxmax(), corrPeakThresh);
      cutLine->SetLineColor(kRed);
      cutLine->Draw();
      cutLine = 0; cutLine = new TLine(0, corrPeakDiagThresh, hilbPeakDiagThresh, 0);
      cutLine->SetLineColor(kRed);
      //cutLine->Draw();
    }
    if (corrSnrCanv) {
      corrSnrCanv->cd(w+1);
      gStyle->SetOptStat("emr");
      corrSnrHist[p]->SetLabelSize(0.04, "xyz");
      gPad->SetLogz();
      corrSnrHist[p]->Draw("COLZ");
      corrSnrCanv->Update();
      cutLine = 0; cutLine = new TLine(hilbPeakThresh,gPad->GetUymin(),hilbPeakThresh,gPad->GetUymax());
      cutLine->SetLineColor(kRed);
      cutLine->Draw();
      cutLine = 0; cutLine = new TLine(gPad->GetUxmin(), snrThreshold, gPad->GetUxmax(), snrThreshold);
      cutLine->SetLineColor(kRed);
      cutLine->Draw();
      //cutLine = 0; cutLine = new TLine(0, corrPeakDiagThresh, hilbPeakDiagThresh, 0);
      //cutLine->SetLineColor(kRed);
      //cutLine->Draw();
    }
    if (hwAngleCanv) {  
      hwAngleCanv->cd(w+1);
      gStyle->SetOptStat("emr");
      hwAngleHist[p]->Draw();
      hwAngleCanv->Update();
      cutLine = 0; cutLine = new TLine(maxHwAngle, pow(10, gPad->GetUymin()), maxHwAngle, pow(10, gPad->GetUymax()));
      cutLine->SetLineColor(kRed);
      cutLine->Draw();
    }
    if (peakRatioCanv) {      
      peakRatioCanv->cd(w+1);
      gStyle->SetOptStat("emr");
      peakRatioHist[p]->SetLabelSize(0.04, "xyz");
      peakRatioHist[p]->Draw();
      peakRatioCanv->Update();
      cutLine = 0; cutLine = new TLine(peakRatioThreshold, gPad->GetUymin(), peakRatioThreshold, gPad->GetUymax());
      cutLine->SetLineColor(kRed);
      cutLine->Draw();
    }
    if (conMapCanv) {
      conMapCanv->cd(w+1);
      gStyle->SetOptStat("e");
      gPad->SetRightMargin(0.15);
      eventHist[p]->SetLabelSize(0.04, "xyz");
      gPad->SetLogz();
      eventHist[p]->GetYaxis()->SetRangeUser(-3000, 3000);
      eventHist[p]->Draw("COLZ");
      conMapCanv->Update();
      int maxBinC = eventHist[p]->GetMaximumBin();
      float zPeakC = eventHist[p]->GetBinContent(maxBinC);
      continentMax = max(continentMax, zPeakC);
    }
    if (anyPolCanv) {
      anyPolCanv->cd(w+1);
      gStyle->SetOptStat("emr");
      anyPolHist[p]->Draw();
      anyPolCanv->Update();
      cutLine = 0; cutLine = new TLine(minAnyPolFrac, gPad->GetUymin(), minAnyPolFrac, gPad->GetUymax());
      cutLine->SetLineColor(kRed);
      cutLine->Draw();
    }
    if (snrCanv0) {
      snrCanv0->cd(w+1);
      gStyle->SetOptStat("emrou");
      snrHist0[p]->Draw();
      snrCanv0->Update();
      cutLine = 0; cutLine = new TLine(minSnr, gPad->GetUymin(), minSnr, gPad->GetUymax());
      cutLine->SetLineColor(kRed);
      cutLine->Draw();
    }
    if (snrCanv1) {
      snrCanv1->cd(w+1);
      gStyle->SetOptStat("emrou");
      snrHist1[p]->Draw();
      snrCanv1->Update();
      cutLine = 0; cutLine = new TLine(minSnr, gPad->GetUymin(), minSnr, gPad->GetUymax());
      cutLine->SetLineColor(kRed);
      cutLine->Draw();
    }
    if (snrCanv2) {
      snrCanv2->cd(w+1);
      gStyle->SetOptStat("emrou");
      snrHist2[p]->Draw();
      snrCanv2->Update();
      cutLine = 0; cutLine = new TLine(minSnr, gPad->GetUymin(), minSnr, gPad->GetUymax());
      cutLine->SetLineColor(kRed);
      cutLine->Draw();
    }
    if (snrCanv3) {
      snrCanv3->cd(w+1);
      gStyle->SetOptStat("emrou");
      snrHist3[p]->Draw();
      snrCanv3->Update();
      cutLine = 0; cutLine = new TLine(minSnr, gPad->GetUymin(), minSnr, gPad->GetUymax());
      cutLine->SetLineColor(kRed);
      cutLine->Draw();
    }
    if (timeCanv) {
      timeCanv->cd(w+1);
      gStyle->SetOptStat("emr");
      timeHist[p]->Draw();
      timeCanv->Update();
    }
    if (timeNsCanv) {
      timeNsCanv->cd(w+1);
      gStyle->SetOptStat("emr");
      timeNsHist[p]->Draw();
      timeNsCanv->Update();
    }
    if (linCircPolCanv) {
      linCircPolCanv->cd(w+1);
      gStyle->SetOptStat("emr");
      linCircPolHist[p]->SetLabelSize(0.04, "xyz");
      linCircPolHist[p]->Draw("COLZ");
      linCircPolCanv->Update();
      cutLine = 0; cutLine = new TLine(minLinPolFrac, gPad->GetUymin(), minLinPolFrac, gPad->GetUymax());
      cutLine->SetLineColor(kRed);
      cutLine->Draw();
      cutLine = 0; cutLine = new TLine(gPad->GetUxmin(), maxCircPolFrac, gPad->GetUxmax(), maxCircPolFrac);
      cutLine->SetLineColor(kRed);
      cutLine->Draw();
    }
    if (linPolCanv) {
      linPolCanv->cd(w+1);
      gStyle->SetOptStat("emr");
      linPolHist[p]->SetLabelSize(0.04, "xyz");
      linPolHist[p]->Draw("COLZ");
      linPolCanv->Update();      
    }
    if (headingCanv && firstPol) {
      headingCanv->cd(0);
      headingHist->SetMarkerStyle(1);
      headingHist->SetTitle("heading vs time");
      headingHist->Draw("AP");
      gStyle->SetOptStat("");
      headingCanv->Update();            
    }
    if (cm3Canv) {
      cm3Canv->cd(w+1);
      gStyle->SetOptStat("emr");
      cm3Hist[p]->Draw("COLZ");
      cm3Canv->Update();
      //cutLine = 0; cutLine = new TLine(minLinPolFrac, gPad->GetUymin(), minLinPolFrac, gPad->GetUymax());
      //cutLine->SetLineColor(kRed);
      //cutLine->Draw();
    }
    if (sunReflCanv) {
      sunReflCanv->cd(w+1);
      gStyle->SetOptStat("emr");
      sunHist[p]->Draw("COLZ");
      sunReflCanv->cd(w+1);
      sunReflHist[p]->Draw("COLZ");
    }
    if (sunCanv) {
      sunCanv->cd(w+1);
      gStyle->SetOptStat("emr");
      sunHist[p]->Draw("COLZ");
      sunCanv->cd(w+1);
      sunHist[p]->Draw("COLZ");
    }
    if (pointingCanvC) {
      pointingCanvC->cd(w+1);
      gPad->SetLogz();
      pointingHist[p]->Draw("COLZ");
      gStyle->SetOptStat("");
      pointingCanvC->Update();
    }
    if (pointingCanvF) {
      pointingCanvF->cd(w+1);
      gPad->SetLogz();
      pointingHistF[p]->Draw("COLZ");
      gStyle->SetOptStat("");
      pointingCanvF->Update();
    }
    if (errSnrCanv) {
      errSnrCanv->cd(w+1);
      gStyle->SetOptStat("e");
      errSnrHist[p]->Draw("COLZ");
      errSnrCanv->Update();
      cutLine = 0; cutLine = new TLine(minSnr, gPad->GetUymin(), minSnr, gPad->GetUymax());
      cutLine->SetLineColor(kRed);
      cutLine->Draw();
    }
    if (ptTimeCanvTh) {
      ptTimeCanvTh->cd(w+1);
      //gStyle->SetOptStat("emr");
      ptTimeHistPhi[p]->SetMarkerStyle(1);
      ptTimeHistPhi[p]->SetTitle("Phi pointing vs time");
      ptTimeHistPhi[p]->Draw("AP");
      ptTimeCanvTh->Update();
    }
    if (calPulseDistCanv) {
      calPulseDistCanv->cd(w+1);
      gStyle->SetOptStat("e");
      calPulseDistHist[p]->Draw();
      calPulseDistCanv->Update();
    }
    if (pointingCanv1) {
      pointingCanv1->cd(2*w+1);
      gStyle->SetOptStat("emrou");
      gPad->SetLogy();
      ptHistPhi[p]->Draw();
      pointingCanv1->Update();
      pointingCanv1->cd(2*w+2);
      gPad->SetLogy();
      ptHistTheta[p]->Draw();
      pointingCanv1->Update();
    }
    if (errSnrCanv1) {
      errSnrCanv1->cd(2*w+1);
      gStyle->SetOptStat("");
      gPad->SetLogy();
      ptHistPhiSnr[p]->Draw("hist");
      errSnrCanv1->Update();
      errSnrCanv1->cd(2*w+2);
      gPad->SetLogy();
      ptHistThetaSnr[p]->Draw("hist");
      errSnrCanv1->Update();
    }
    if (acceptSnrCanv) {
      acceptSnrCanv->cd(w+1);
      acceptHistSnr[p]->SetMarkerStyle(kFullTriangleUp);
      acceptHistSnr[p]->SetMarkerSize(0.6);
      acceptHistSnr[p]->Draw();
      acceptHistSnr[p]->GetYaxis()->SetRangeUser(0, 1.05);
      acceptHistSnr[p]->Draw();
      gStyle->SetOptStat("");
      acceptSnrCanv->Update();
    }
    if (distSnrEffCanv) {
      distSnrEffCanv->cd(w+1);
      gStyle->SetOptStat("");
      //distSnrEff[p]->SetMarkerStyle(kFullTriangleUp);
      //distSnrEff[p]->SetMarkerSize(0.6);
      distSnrEff[p]->Divide(distSnrAll[p]);
      distSnrEff[p]->Draw("COLZ");
      distSnrEffCanv->Update();
    }
    if (effSnrPFFracCanv) {
      effSnrPFFracCanv->cd(2*w+1);
      gStyle->SetOptStat("");
      effSnrPFFrac[p]->Divide(totSnrPFFrac[p]);
      totSnrPFFrac[p]->Draw("COLZ");
      effSnrPFFracCanv->cd(2*w+2);
      effSnrPFFrac[p]->Draw("COLZ");
      effSnrPFFracCanv->Update();
    }
    if (corrDistCanv) {
      corrDistCanv->cd(w+1);
      gStyle->SetOptStat("eou");
      corrDistHist[p]->Draw("COLZ");
      corrDistCanv->Update();
    }
    if (hilbDistCanv) {
      hilbDistCanv->cd(w+1);
      gStyle->SetOptStat("");
      hilbDistHist[p]->Draw("COLZ");
      hilbDistCanv->Update();
    }
    if (linPolDistCanv) {
      linPolDistCanv->cd(w+1);
      gStyle->SetOptStat("");
      linPolDistHist[p]->Draw("COLZ");
      linPolDistCanv->Update();
    }
    if (distPFFracCanv) {
      distPFFracCanv->cd(2*w+1);
      distPFFrac_0[p]->Draw("COLZ");
      gStyle->SetOptStat("");
      distPFFracCanv->cd(2*w+2);
      distPFFrac_1[p]->Draw("COLZ");
      gStyle->SetOptStat("");
      distPFFracCanv->Update();
    }

    if (linPolPlPhiCanv) {
      linPolPlPhiCanv->cd(2*w+1);
      linPolPlPhiHist[p]->Draw("COLZ"); 
      linPolPlPhiCanv->cd(2*w+2);
      linPolPlNPhiHist[p]->Draw("COLZ");
    }
    if (lesserCircPeakSNRCanv) {
      lesserCircPeakSNRCanv->cd(w+1);
      lesserCircPeakSNRHist[p]->Draw("COLZ");
      gStyle->SetOptStat("");
      lesserCircPeakSNRCanv->Draw();    
    }
    if (circPeakSepSNRCanv) {
      circPeakSepSNRCanv->cd(w+1);
      circPeakSepSNRHist[p]->Draw("COLZ");
      gStyle->SetOptStat("");
      circPeakSepSNRCanv->Draw();        
    }
    if (lesserCircPeakLDCanv) {
      lesserCircPeakLDCanv->cd(w+1);
      lesserCircPeakLDHist[p]->Draw("COLZ");
      gStyle->SetOptStat("");
      lesserCircPeakLDCanv->Draw();    
    }
    if (circPeakSepLDCanv) {
      circPeakSepLDCanv->cd(w+1);
      circPeakSepLDHist[p]->Draw("COLZ");
      gStyle->SetOptStat("");
      circPeakSepLDCanv->Draw();        
    }
    
    if (minCircValCanv) {
      minCircValCanv->cd(w+1);
      minCircValHist[p]->Draw();
      //gStyle->SetOptStat("");
      //minCircValCanv->Draw();        
    }
        
    if (peakTimeCanv) {
      peakTimeCanv->cd(w+1);
      peakTimeHist[p]->Draw("COLZ");
      gStyle->SetOptStat("");
      peakTimeCanv->Draw();        
    }

    if (firstPol) {
      if (stokesCanv) {
        int cNum = 0;
        //stokesCanv->cd(++cNum);
        //stokesIHist->Draw();
        gStyle->SetOptStat("emrou");
        stokesCanv->cd(++cNum);
        stokesQHist->Draw();
        gStyle->SetOptStat("emrou");
        stokesCanv->cd(++cNum);
        stokesUHist->Draw();
        gStyle->SetOptStat("emrou");
        stokesCanv->cd(++cNum);
        stokesVHist->Draw();
        gStyle->SetOptStat("emrou");
      }
    }
    firstPol = false;
  }
  // cout << "skyMapMax = " << skyMapMax << endl;
  //lprintf ("current directory is %s \n", gDirectory->GetName());

  if (corrPeakPolCanv) {
    corrPeakPolCanv->cd(0);
    corrPeakPolHist->Draw("COLZ");
    cutLine = new TLine(0, 0, 0.3, 0.3);
    cutLine->SetLineColor(kRed);
    cutLine->Draw();
    gStyle->SetOptStat("");
    corrPeakCanv->Draw();  
  }
  
  if (corrPeakLesserPolCanv) {
    corrPeakLesserPolCanv->cd(0);
    corrPeakLesserPolHist->Sumw2(false);  
    corrPeakLesserPolHist->Draw();
    corrPeakLesserPolCanv->Draw();
  }
  
  if (peakSepCircCanv) {
    peakSepCircCanv->cd(0);
    gPad->SetLogy();
    peakSepCircHist->Draw();
    gStyle->SetOptStat("");
    peakSepCircCanv->Draw();
  }
  
  if (lesserPeakCircCanv) {
    lesserPeakCircCanv->cd(0);
    gPad->SetLogy();
    lesserPeakCircHist->Draw("");
    gStyle->SetOptStat("");
    lesserPeakCircCanv->Draw();
  }
  
  if (circPeakCanv) {
    circPeakCanv->cd(0);
    circPeakHist->Draw("COLZ");
    gStyle->SetOptStat("");
    circPeakCanv->Draw();
  }
  
  if (cPolPeakDistCanv) {
    cPolPeakDistCanv->cd(0);
    cPolPeakDistHist->Draw("COLZ");
    gStyle->SetOptStat("");
    cPolPeakDistCanv->Draw();
  }
  

  // redraw maps to common scale, where needed
  for (int w=0; w<whichPols.size(); ++w) {
    int p = whichPols[w];
    if (skyMapCanv) {
      skyMapCanv->cd(w+1);
      gStyle->SetOptStat("e");
      skyMapHist[p]->GetZaxis()->SetRangeUser(0, skyMapMax);
      skyMapHist[p]->Draw("COLZ");
      //skyMapHist[p]->Draw("");
    }
    if (skyMapPlCanv) {
      skyMapPlCanv->cd(w+1);
      gStyle->SetOptStat("e");
      skyMapPlHist[p]->GetZaxis()->SetRangeUser(0, skyMapMax);
      skyMapPlHist[p]->Draw("COLZ");
    }
    if (conMapCanv) {
      conMapCanv->cd(w+1);
      gStyle->SetOptStat("e");
      coastLineGr->SetMarkerColor(14);
      coastLineGr->SetFillColor(18);
      coastLineGr->SetLineWidth(1);
      coastLineGr->SetFillStyle(1001);
      eventHist[p]->GetZaxis()->SetRangeUser(0, continentMax);
      eventHist[p]->Draw("COLZ SAME");
      coastLineGr->Draw("C");
    }
        //int ec, nc;
    double ea, no;
    //(new TLine(-1000000, -1000000, 1000000, 1000000))->Draw();   // for squaring the graph
    //TMarker* thisMarker;
    for (TMarker* thisMarker : conMapMarkers) {thisMarker->Draw();}
    
    if (plMapPoints.size() > 1) {
      plPathGr->Draw("C");
      /*
      for (int k=0; k<plMapPoints.size()-1; ++k) {
        double ea0, no0, ea1, no1;
        bedmap->LonLattoEaNo(plMapPoints[k][1], plMapPoints[k][0], ea0, no0);
        bedmap->LonLattoEaNo(plMapPoints[k+1][1], plMapPoints[k+1][0], ea1, no1);
        TLine* thisLine = new TLine(ea0, no0, ea1, no1);
        thisLine->SetLineColor(kRed);
        thisLine->Draw();      
      } 
      */
    } else {
      bedmap->LonLattoEaNo(plMapPoints[0][1], plMapPoints[0][0], ea, no);      
      TMarker* thisMarker = new TMarker(ea, no, kPlus);
      thisMarker->SetMarkerColor(kRed);
      thisMarker->Draw();
    }
    //TEllipse* thisCirc = new TEllipse(0, 0, no);
    deg80Circ->Draw();
    deg70Circ->Draw();

    TGaxis* scale = bedmap->distanceScale(-2000000, 0, -2000000, -2000000);
    scale->SetLabelSize(0.02);
    scale->SetTitleSize(0.02);
    scale->SetLabelFont(42);
    scale->SetTitleFont(42);
    scale->Draw();
    if (hpCountsCanv) {
      hpCountsCanv->cd(p+1);
      gPad->SetLogz();
      gPad->SetRightMargin(0.15);
      hpCountsHist[p]->SetLabelSize(0.04, "xyz");    
      hpCountsHist[p]->Draw("COLZ");
      coastLineGr->SetLineWidth(1);
      coastLineGr->Draw("C");
      deg80Circ->Draw();
      deg70Circ->Draw();
      plPathGr->Draw("C");
      for (TMarker* thisMarker : conMapMarkers) {thisMarker->Draw();}
      // write the healpix bin numbers
      // iterate through the healpix count bins and write a text to each occupied bin
      for (pair<int, float> thisEntry : hpEventCounts[p]) {
        pointing thisPoint = healpix->pix2ang(thisEntry.first);
        double thisLat = 90.0 - (thisPoint.theta * 180.0/M_PI) ;
        double thisLon = thisPoint.phi * 180.0/M_PI;
        double ea, no;
        bedmap->LonLattoEaNo(thisLon, thisLat,	ea, no);
        ea /= 1000;
        no /= 1000;
        lprintf("hp bin %i lat,lon is %f,%f  ea,no is %f,%f \n", thisEntry.first, thisLat, thisLon, ea, no);
        char hpBinText[8]; sprintf(hpBinText, "%i", thisEntry.first);
        //TPaveText* thisText = new TPaveText(ea, no, ea+100, no+20, "NB NDC");
        lprintf("drawing hp bin number at %f, %f, %f, %f \n", ea-125, no, ea+125, no+100);
        bool posQuad = (ea*no > 0);
        double textEa = posQuad ? ea : ea-75;
        double textNo = posQuad ? no-125 : no+100;
        TText* thisText = new TText(textEa, textNo, hpBinText);
        thisText->SetTextSize(0.02);
        thisText->SetTextAngle(posQuad ? 55 : -55);
        thisText->SetTextFont(42);
        //thisText->SetBorderSize(0);
        //thisText->SetFillColorAlpha(0, 0.5);
        //thisText->AddText(hpBinText);
        thisText->Draw();
        //printf("   bin %8i: %10.3f events.      theta=%4.3f, phi=%4.3f, lat=%6.3f, lon=%6.3f\n", 
        //        thisEntry.first, thisEntry.second, thisPoint.theta, thisPoint.phi, thisLat, thisLon);
      }
    }    
    if ((singleEventNumber > 0 || maxEvents==1) && p==0) {
      gStyle->SetPalette(100, pdfPalette);
      conMapCanvF->cd(0);
      gStyle->SetOptStat("");
      gPad->Range(srcEa[p]-errEllipseParms1[p][0]*1.4, srcNo[p]-errEllipseParms1[p][0]*1.4, 
              srcEa[p]+errEllipseParms1[p][0]*1.4, srcNo[p]+errEllipseParms1[p][0]*1.4);
      
      if (errEllipseParms1[p][0]>0) {
        float eCtrEa = srcEa[p] + errEllipseParms0[p][0]*cos(plSrcPhi[p]);
        float eCtrNo = srcNo[p] + errEllipseParms0[p][0]*sin(plSrcPhi[p]);
        vbprintf(" error ellipse center is %f, %f \n", eCtrEa, eCtrNo);
        TEllipse* errorEllipse1 = new TEllipse(eCtrEa, eCtrNo, errEllipseParms0[p][1], errEllipseParms0[p][2], 0, 360, plSrcPhi[p]*180.0/M_PI);
        errorEllipse1->SetFillStyle(0);
        errorEllipse1->SetLineColor(kRed);
        //errorEllipse1->Draw("A");

        vbprintf("  ellipse parameters: rAxis=%f, thAxis=%f, centerX=%f, centerY=%f \n", errEllipseParms1[p][0], 
                errEllipseParms1[p][1], errEllipseParms1[p][2], errEllipseParms1[p][3]);
        TEllipse* errorEllipse2 = new TEllipse(errEllipseParms1[p][2], errEllipseParms1[p][3], 
                errEllipseParms1[p][0], errEllipseParms1[p][1], 0, 360, plSrcPhi[p]*180.0/M_PI);
        errorEllipse2->SetFillStyle(0);
        errorEllipse2->SetLineColor(kRed);
        errorEllipse2->Draw();
        (new TLine(plEa, plNo, srcEa[p], srcNo[p]))->Draw();
        
        float errHexX[7], errHexY[7];
        for (int k=0; k<6; ++k) {errHexX[k] = errHexVertices[p][0][k]; errHexY[k] = errHexVertices[p][1][k];} 
        errHexX[6] = errHexX[0]; errHexY[6] = errHexY[0];
        TPolyLine* errorHex = new TPolyLine(7, errHexX, errHexY);
        errorHex->SetLineColor(kBlue);
        errorHex->SetFillStyle(0);
        //errorHex->Draw();
        
        vbprintf("drawing error ellipse plot %cPol \n", pChars[p]);
        vbprintf(" old way \n");
        for (int k=0; k<errHexPixels[p][0].size(); ++k) {vbprintf("   %f, %f, %f \n", errHexPixels[p][0][k], errHexPixels[p][1][k], errHexPixels[p][2][k]);}
        vbprintf(" new way \n");
        for (int k=0; k<errHexRayTrace[p][0].size(); ++k) {vbprintf("   %f, %f, %f \n", errHexRayTrace[p][0][k], errHexRayTrace[p][1][k], errHexRayTrace[p][2][k]);}
        bool errorEllipseNewWay = true;
        if (errorEllipseNewWay) {
          for (int k=0; k<errHexRayTrace[p][0].size(); ++k) {
            TMarker* pixelMarker = new TMarker(errHexRayTrace[p][0][k], errHexRayTrace[p][1][k], kFullDotLarge);
            int colorIndex = 100 * errHexRayTrace[p][2][k]/0.03;
            vbprintf("   %f, %f, %f   colorIndex=%i \n", errHexRayTrace[p][0][k], errHexRayTrace[p][1][k], errHexRayTrace[p][2][k], colorIndex);
            colorIndex = (colorIndex>100 ? 100 : colorIndex);
            pixelMarker->SetMarkerColor(pdfPalette[colorIndex]);
            pixelMarker->Draw();
          } 
        } else {
          for (int k=0; k<errHexPixels[p][0].size(); ++k) {
            TMarker* pixelMarker = new TMarker(errHexPixels[p][0][k], errHexPixels[p][1][k], kPlus);
            pixelMarker->SetMarkerColor(kBlue);
            pixelMarker->Draw();
          } 
        }
        TGaxis* scale1 = bedmap->distanceScale(srcEa[p], srcEa[p]+errEllipseParms1[p][0]*0.6, srcNo[p]-errEllipseParms1[p][0]*0.5, srcNo[p]-errEllipseParms1[p][0]*0.5);
        scale1->SetLabelSize(0.02);
        scale1->SetTitleSize(0.02);
        scale1->SetLabelFont(42);
        scale1->SetTitleFont(42);
        scale1->Draw();
      }
    }
    /*
    bedmap->SurfaceLonLattoEN(MCMURDO[1], MCMURDO[0], ec, nc);
    bedmap->ENtoEaNo(ec, nc, ea, no);
    (new TMarker(ea, no, kPlus))->Draw();
    bedmap->SurfaceLonLattoEN(WAIS_DIVIDE[1], WAIS_DIVIDE[0], ec, nc);
    bedmap->ENtoEaNo(ec, nc, ea, no);
    (new TMarker(ea, no, kPlus))->Draw();
    bedmap->SurfaceLonLattoEN(VOSTOK[1], VOSTOK[0], ec, nc);
    bedmap->ENtoEaNo(ec, nc, ea, no);
    (new TMarker(ea, no, kPlus))->Draw();
    bedmap->SurfaceLonLattoEN(SOUTH_POLE[1], SOUTH_POLE[0], ec, nc);
    bedmap->ENtoEaNo(ec, nc, ea, no);
    (new TMarker(ea, no, kPlus))->Draw();
    */
  }
  //vbprintf("single event number is %i, entry number is %i \n", singleEventNumber, startEntry);
  /*
  TFile* outputFile = new TFile(outputFilename, "RECREATE");
  lprintf ("result tree 0 has %llu entries \n", outputTree0->GetEntries());
  lprintf ("result tree 1 has %llu entries \n", outputTree1->GetEntries());
  lprintf ("current directory is %s \n", gDirectory->GetName());
  
  if (saveOutput) {
    outputTree0->SetDirectory(outputFile);
    outputTree1->SetDirectory(outputFile);
    outputTree0->Write();
    outputTree1->Write();
    outputFile->Save();
    outputFile->Close();
  }
  */
  // save the canvii to files
  if (saveOutput) {
    char canvFilename[1024];
    if (corrPeakCanv) {
      sprintf(canvFilename, "%s/corrPeak_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      corrPeakCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/corrPeak_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      corrPeakCanv->SaveAs(canvFilename);
    }
    if (hilbPeakCanv) {
      sprintf(canvFilename, "%s/hilbPeak_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      hilbPeakCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/hilbPeak_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      hilbPeakCanv->SaveAs(canvFilename);
    }
    if (anyPolCanv) {
      sprintf(canvFilename, "%s/anyPol_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      anyPolCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/anyPol_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      anyPolCanv->SaveAs(canvFilename);
    }
    if (snrCanv0) {
      sprintf(canvFilename, "%s/snr0_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      snrCanv0->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/snr0_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      snrCanv0->SaveAs(canvFilename);
    }
    if (snrCanv1) {
      sprintf(canvFilename, "%s/snr1_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      snrCanv1->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/snr1_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      snrCanv1->SaveAs(canvFilename);
    }
    if (snrCanv2) {
      sprintf(canvFilename, "%s/snr2_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      snrCanv2->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/snr2_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      snrCanv2->SaveAs(canvFilename);
    }
    if (snrCanv3) {
      sprintf(canvFilename, "%s/snr3_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      snrCanv3->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/snr3_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      snrCanv3->SaveAs(canvFilename);
    }
    if (sunCanv) {  
      sprintf(canvFilename, "%s/sun_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      sunCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/sun_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      sunCanv->SaveAs(canvFilename);
    }
    if (timeCanv) {
      sprintf(canvFilename, "%s/time_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      timeCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/time_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      timeCanv->SaveAs(canvFilename);
    }
    if (timeNsCanv) {
      sprintf(canvFilename, "%s/timeNs_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      timeNsCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/timeNs_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      timeNsCanv->SaveAs(canvFilename);
    }
    if (calPulseDistCanv) {
      sprintf(canvFilename, "%s/calPulseDist_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      calPulseDistCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/calPulseDist_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      calPulseDistCanv->SaveAs(canvFilename);
    }
    if (hwAngleCanv) {
      sprintf(canvFilename, "%s/hwAngle_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      hwAngleCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/hwAngle_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      hwAngleCanv->SaveAs(canvFilename);
    }
    if (linCircPolCanv) {
      sprintf(canvFilename, "%s/linCircPol_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      linCircPolCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/linCircPol_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      linCircPolCanv->SaveAs(canvFilename);
    }
    if (skyMapCanv) {
      sprintf(canvFilename, "%s/skyMap_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      skyMapCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/skyMap_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      skyMapCanv->SaveAs(canvFilename);
    }
    if (skyMapPlCanv) {
      sprintf(canvFilename, "%s/skyMapPl_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      skyMapPlCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/skyMapPl_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      skyMapPlCanv->SaveAs(canvFilename);
    }  

    if (skyMapCircCanv) {
      sprintf(canvFilename, "%s/skyMapCirc_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pCharsC);
      skyMapCircCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/skyMapCirc_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pCharsC);
      skyMapCircCanv->SaveAs(canvFilename);
    }
    if (skyMapPlCircCanv) {
      sprintf(canvFilename, "%s/skyMapPlCirc_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pCharsC);
      skyMapPlCircCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/skyMapPlCirc_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pCharsC);
      skyMapPlCircCanv->SaveAs(canvFilename);
    }  

    if (corrHilbCanv) {
      sprintf(canvFilename, "%s/corrHilb_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      corrHilbCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/corrHilb_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      corrHilbCanv->SaveAs(canvFilename);
    }
    if (corrSnrCanv) {
      sprintf(canvFilename, "%s/corrSnr_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      corrSnrCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/corrSnr_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      corrSnrCanv->SaveAs(canvFilename);
    }
    if (conMapCanv) {
      sprintf(canvFilename, "%s/conMap_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      conMapCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/conMap_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      conMapCanv->SaveAs(canvFilename);
    }
    if (conMapCanvF) {
      sprintf(canvFilename, "%s/conMapF_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      conMapCanvF->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/conMapF_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      conMapCanvF->SaveAs(canvFilename);
    }
    if (peakRatioCanv) {
      sprintf(canvFilename, "%s/peakRatio_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      peakRatioCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/peakRatio_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      peakRatioCanv->SaveAs(canvFilename);
    }
    if (cm3Canv) {
      sprintf(canvFilename, "%s/cm3_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      cm3Canv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/cm3_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      cm3Canv->SaveAs(canvFilename);
    }
    if (acceptSnrCanv) {
      sprintf(canvFilename, "%s/acceptSnr_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      acceptSnrCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/acceptSnr_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      acceptSnrCanv->SaveAs(canvFilename);
    }
    if (pointingCanvC) {
      sprintf(canvFilename, "%s/pointingCanvC_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      pointingCanvC->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/pointingCanvC_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      pointingCanvC->SaveAs(canvFilename);
    }
    if (pointingCanvF) {
      sprintf(canvFilename, "%s/pointingCanvF_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      pointingCanvF->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/pointingCanvF_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      pointingCanvF->SaveAs(canvFilename);
    }
    if (errSnrCanv) {
      sprintf(canvFilename, "%s/errSnrCanv_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      errSnrCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/errSnrCanv_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      errSnrCanv->SaveAs(canvFilename);
    }

    if (ptTimeCanvTh) {
      sprintf(canvFilename, "%s/ptTimeCanvTh_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      ptTimeCanvTh->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/ptTimeCanvTh_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      ptTimeCanvTh->SaveAs(canvFilename);
    }
    if (calPulseDistCanv) {
      sprintf(canvFilename, "%s/calPulseDistCanv_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      calPulseDistCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/calPulseDistCanv_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      calPulseDistCanv->SaveAs(canvFilename);
    }
    if (pointingCanv1) {
      sprintf(canvFilename, "%s/pointingCanv1_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      pointingCanv1->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/pointingCanv1_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      pointingCanv1->SaveAs(canvFilename);
    }
    if (errSnrCanv1) {
      sprintf(canvFilename, "%s/errSnrCanv1_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      errSnrCanv1->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/errSnrCanv1_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      errSnrCanv1->SaveAs(canvFilename);
    }
    if (distSnrEffCanv) {
      sprintf(canvFilename, "%s/distSnrEffCanv_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      distSnrEffCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/distSnrEffCanv_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      distSnrEffCanv->SaveAs(canvFilename);
    }
    if (corrDistCanv) {
      sprintf(canvFilename, "%s/corrDistCanv_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      corrDistCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/corrDistCanv_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      corrDistCanv->SaveAs(canvFilename);
    }
    if (hilbDistCanv) {
      sprintf(canvFilename, "%s/hilbDistCanv_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      hilbDistCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/hilbDistCanv_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      hilbDistCanv->SaveAs(canvFilename);
    }
    if (linPolDistCanv) {
      sprintf(canvFilename, "%s/linPolDistCanv_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      linPolDistCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/linPolDistCanv_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      linPolDistCanv->SaveAs(canvFilename);
    }
    if (effSnrPFFracCanv) {
      sprintf(canvFilename, "%s/effSnrPFFracCanv_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      effSnrPFFracCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/effSnrPFFracCanv_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      effSnrPFFracCanv->SaveAs(canvFilename);
    }
    if (distPFFracCanv) {
      sprintf(canvFilename, "%s/distPFFracCanv_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      distPFFracCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/distPFFracCanv_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      distPFFracCanv->SaveAs(canvFilename);
    }
    if (linPolCanv) {
      sprintf(canvFilename, "%s/linPolCanv_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      linPolCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/linPolCanv_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      linPolCanv->SaveAs(canvFilename);      
    }
    if (linPolPlPhiCanv) {
      sprintf(canvFilename, "%s/linPolPlPhiCanv_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      linPolPlPhiCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/linPolPlPhiCanv_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      linPolPlPhiCanv->SaveAs(canvFilename);      
    }
    if (corrPeakPolCanv) {
      sprintf(canvFilename, "%s/corrPeakPolCanv_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      corrPeakPolCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/corrPeakPolCanv_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      corrPeakPolCanv->SaveAs(canvFilename);      
    }    
    if (corrPeakLesserPolCanv) {
      sprintf(canvFilename, "%s/corrPeakLesserPolCanv_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      corrPeakLesserPolCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/corrPeakLesserPolCanv_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      corrPeakLesserPolCanv->SaveAs(canvFilename);      
    }    
    if (hpCountsCanv) {
      sprintf(canvFilename, "%s/hpCountsCanv_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      hpCountsCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/hpCountsCanv_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      hpCountsCanv->SaveAs(canvFilename);      
    }    
    if (peakSepCircCanv) {
      sprintf(canvFilename, "%s/peakSepCircCanv_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      peakSepCircCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/peakSepCircCanv_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      peakSepCircCanv->SaveAs(canvFilename);      
    }    
    if (lesserPeakCircCanv) {
      sprintf(canvFilename, "%s/lesserPeakCircCanv_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      lesserPeakCircCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/lesserPeakCircCanv_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      lesserPeakCircCanv->SaveAs(canvFilename);      
    }    
    if (circPeakCanv) {
      sprintf(canvFilename, "%s/circPeakCanv_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      circPeakCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/circPeakCanv_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      circPeakCanv->SaveAs(canvFilename);      
    }    
    if (cPolPeakDistCanv) {
      sprintf(canvFilename, "%s/cPolPeakDistCanv_%03i_%03i_%s.png", outputDir, startRunNumber, endRunNumber, pChars);
      cPolPeakDistCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/cPolPeakDistCanv_%03i_%03i_%s.root", outputDir, startRunNumber, endRunNumber, pChars);
      cPolPeakDistCanv->SaveAs(canvFilename);      
    }   
    
    if (lesserCircPeakSNRCanv) {
      sprintf(canvFilename, "%s/lesserCircPeakSNRCanv_%03i_%03i.png", outputDir, startRunNumber, endRunNumber);
      lesserCircPeakSNRCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/lesserCircPeakSNRCanv_%03i_%03i.root", outputDir, startRunNumber, endRunNumber);
      lesserCircPeakSNRCanv->SaveAs(canvFilename);      
    }
    if (circPeakSepSNRCanv) {
      sprintf(canvFilename, "%s/circPeakSepSNRCanv_%03i_%03i.png", outputDir, startRunNumber, endRunNumber);
      circPeakSepSNRCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/circPeakSepSNRCanv_%03i_%03i.root", outputDir, startRunNumber, endRunNumber);
      circPeakSepSNRCanv->SaveAs(canvFilename);      
    }
    if (lesserCircPeakLDCanv) {
      sprintf(canvFilename, "%s/lesserCircPeakLDCanv_%03i_%03i.png", outputDir, startRunNumber, endRunNumber);
      lesserCircPeakLDCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/lesserCircPeakLDCanv_%03i_%03i.root", outputDir, startRunNumber, endRunNumber);
      lesserCircPeakLDCanv->SaveAs(canvFilename);      
    }
    if (circPeakSepLDCanv) {
      sprintf(canvFilename, "%s/circPeakSepLDCanv_%03i_%03i.png", outputDir, startRunNumber, endRunNumber);
      circPeakSepLDCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/circPeakSepLDCanv_%03i_%03i.root", outputDir, startRunNumber, endRunNumber);
      circPeakSepLDCanv->SaveAs(canvFilename);      
    }
    if (minCircValCanv) {
      sprintf(canvFilename, "%s/minCircValCanv_%03i_%03i.png", outputDir, startRunNumber, endRunNumber);
      minCircValCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/minCircValCanv_%03i_%03i.root", outputDir, startRunNumber, endRunNumber);
      minCircValCanv->SaveAs(canvFilename);      
    }
    if (peakTimeCanv) {
      sprintf(canvFilename, "%s/peakTimeCanv_%03i_%03i.png", outputDir, startRunNumber, endRunNumber);
      peakTimeCanv->SaveAs(canvFilename);
      sprintf(canvFilename, "%s/peakTimeCanv_%03i_%03i.root", outputDir, startRunNumber, endRunNumber);
      peakTimeCanv->SaveAs(canvFilename);      
    }

  }
  fclose(logFile);
  if (displayGraphics) app->Run();
  
  return 0;
}

// return center (Ea, No), radial and lateral semi-axes
vector<float> surfaceEllipseParms1(float phi, float theta, float plSrcDist, float errPh, float errTh, 
        UsefulAdu5Pat* gps, BedmapReader* bedmap) {
  // use bedmap to obtain ellipse left, right, near and far vertices, easting/northing
  vbprintf("phi=%f, theta=%f, dist=%f, errPhi=%f, errTheta=%f \n", phi, theta, plSrcDist, errPh, errTh);
  double thisLat, thisLon, thisAlt, thisAdj;
  double far[2], near[2]; 
  vector<float> fail = {-9999, -9999, -9999, -9999};
  // pointing peak (just for test)
  gps->traceBackToContinent(phi *M_PI/180.0, theta *M_PI/180.0, &thisLon, &thisLat, &thisAlt, &thisAdj);
  if (thisLat == -9999 || thisLon == -9999 || thisAlt == -9999) {
    vbprintf("trace back failed for peak location \n");
    return fail;   // structure violation yeah yeah
  } 
  // far vertex
  gps->traceBackToContinent(phi *M_PI/180.0, (theta-errTh) *M_PI/180.0, &thisLon, &thisLat, &thisAlt, &thisAdj);
  if (thisLat == -9999 || thisLon == -9999 || thisAlt == -9999) {
    vbprintf("trace back failed (%.2f,%.2f) for far ellipse vertex theta-far=%f, theta-near=%f \n", thisLat, thisLon, theta-errTh, theta+errTh);
    return fail;   // structure violation yeah yeah
  } 
  vbprintf(" far vertex  lat=%.2f, lon=%.2f, alt=%.2f, adj=%.2f \n", thisLat, thisLon, thisAlt, thisAdj);
  bedmap->LonLattoEaNo(thisLon, thisLat, far[0], far[1]);
  // near vertex
  gps->traceBackToContinent(phi *M_PI/180.0, (theta+errTh) *M_PI/180.0, &thisLon, &thisLat, &thisAlt, &thisAdj);
  if (thisLat == -9999 || thisLon == -9999 || thisAlt == -9999) {
    vbprintf("trace back failed for near ellipse vertex theta-far=%f, theta-near=%f \n", theta-errTh, theta+errTh);
    return fail;   // structure violation
  }
  vbprintf(" near vertex lat=%.2f, lon=%.2f, alt=%.2f, adj=%.2f \n", thisLat, thisLon, thisAlt, thisAdj);
  bedmap->LonLattoEaNo(thisLon, thisLat, near[0], near[1]);

  float rAxis = 0.5*sqrt((far[0]-near[0])*(far[0]-near[0]) + (far[1]-near[1])*(far[1]-near[1]));
  float centerX = (far[0]+near[0])*0.5;
  float centerY = (far[1]+near[1])*0.5;
  float thAxis = plSrcDist * errPh * M_PI/180.0;;  // errPh assumed small
  return vector<float>({rAxis, thAxis, centerX, centerY});
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
  else if (keyword.compare("--MIN_PEAK_THETA")==0) {minPeakTheta = stof(parm);}
  else if (keyword.compare("--MAX_HW_ANGLE")==0) {maxHwAngle = stof(parm);}
  else if (keyword.compare("--MIN_LIN_POL_FRAC")==0) {minLinPolFrac = stof(parm);}
  else if (keyword.compare("--MIN_ANY_POL_FRAC")==0) {minAnyPolFrac = stof(parm);}
  else if (keyword.compare("--MAX_CIRC_POL_FRAC")==0) {maxCircPolFrac = stof(parm);}
  else if (keyword.compare("--MIN_SNR")==0) {minSnr = stof(parm);}
  else if (keyword.compare("--CORR_PEAK_THRESHOLD")==0) {corrPeakThresh = stof(parm);}
  else if (keyword.compare("--HILB_PEAK_THRESHOLD")==0) {hilbPeakThresh = stof(parm);}
  else if (keyword.compare("--CORR_PEAK_DIAG_THRESHOLD")==0) {corrPeakDiagThresh = stof(parm);}
  else if (keyword.compare("--HILB_PEAK_DIAG_THRESHOLD")==0) {hilbPeakDiagThresh = stof(parm);}
  else if (keyword.compare("--MAX_DEV_FROM_SOUTH")==0) {southCutDist = stof(parm);}
  else if (keyword.compare("--TRIG_PHISEC_PROX")==0) {trigPhisecProx = stoi(parm);}
  else if (keyword.compare("--CUT_ORDER_LIMIT")==0) {aCutOrderStage1Limit = stoi(parm);}
  else if (keyword.compare("--POL_FRAC_METHOD")==0) {polFracMethod = stoi(parm);}
  else if (keyword.compare("--USE_ABBYS_CUT_VALUES")==0) {
    surfSatThreshold = 1000;
    dcOffsetThreshold = 50;
    minWfLength = 240;
    nadirNoiseThreshold = 0.5;
    minPeakTheta = -35;
    minLinPolFrac = 0.3;
    corrPeakThresh = 0.075;
    hilbPeakThresh = 15;
    corrPeakDiagThresh = 0.2;
    hilbPeakDiagThresh = 70;
    peakRatioThreshold = 0.9;
    trigPhisecProx = 1;
  }
  else result = -1;
  vcout << "keyword parameter " << kwp << " result code " << result << endl;
  return result;
}

int canANITASeeStripe(double lon_anita,double lat_anita,double lat_satellite,double lon_satellite,double height_anita,double height_satellite) {

  double phi_anita = M_PI/2. - lon_anita*M_PI/180.;
  double theta_anita = M_PI/2. - lat_anita*M_PI/180.;
  double theta_satellite=M_PI/2.-lat_satellite*M_PI/180.;

  double phi_satellite = M_PI/2. - lon_satellite*M_PI/180.;
  double nx=(R_EARTH+height_anita)*sin(theta_anita)*cos(phi_anita)-(R_EARTH+height_satellite)*sin(theta_satellite)*cos(phi_satellite);
  //cout << "term 1 is " << (R_EARTH+height_anita)*sin(theta_anita)*cos(phi_anita) << "\n";
  //cout << "term 2 is " << (R_EARTH+height_satellite)*sin(theta_satellite)*cos(phi_satellite) << "\n";
  double ny=(R_EARTH+height_anita)*sin(theta_anita)*sin(phi_anita)-(R_EARTH+height_satellite)*sin(theta_satellite)*sin(phi_satellite);
  double nz=(R_EARTH+height_anita)*cos(theta_anita)-(R_EARTH+height_satellite)*cos(theta_satellite);

  double norm=sqrt(nx*nx+ny*ny+nz*nz);

  nx=nx/norm;
  ny=ny/norm;
  nz=nz/norm;
  //cout << "phi_satellite, phi_anita are " << phi_satellite << "\t" << phi_anita << "\n";
  //cout << "theta_satellite, theta_anita are " << theta_satellite << "\t" << theta_anita << "\n";
  //cout << "nx, ny, nz are " << nx << "\t" << ny << "\t" << nz << "\n";


  // this is the distance along the trajectory where the distance to center of earth is a minimum
  double s_0=-1.*(R_EARTH+height_satellite)*
      (nx*sin(theta_satellite)*cos(phi_satellite) +
       ny*sin(theta_satellite)*sin(phi_satellite) +
       nz*cos(theta_satellite));

  // cout << "s_0 is " << s_0 << "\n";
  //cout << "norm is " << norm << "\n";
  if (s_0>norm) {
    return 1;
  }
  else {
    //double n_vec[3]={nx,ny,nz};
    double lengthsquared= -2.*(R_EARTH+height_satellite)
                             *(R_EARTH+height_satellite)
                             *sin(theta_satellite)
                             *( nx*ny*sin(theta_satellite)*sin(phi_satellite)*cos(phi_satellite)
                               +nx*nz*cos(theta_satellite)*cos(phi_satellite)
                               +ny*nz*cos(theta_satellite)*sin(phi_satellite) )
                          + (R_EARTH+height_satellite)
                             *(R_EARTH+height_satellite)
                             *(1.- ( nx*nx*sin(theta_satellite)*sin(theta_satellite)*cos(phi_satellite)*cos(phi_satellite)
                                    +ny*ny*sin(theta_satellite)*sin(theta_satellite)*sin(phi_satellite)*sin(phi_satellite)
                                    +nz*nz*cos(theta_satellite)*cos(theta_satellite) )
                              );

    // cout << "radius of earth squared is " << R_EARTH*R_EARTH << "\n";
    //cout << "lengthsquared is " << lengthsquared << "\n";
    //cout << "lengthsquared is " << lengthsquared << "\n";
    if (lengthsquared>R_EARTH*R_EARTH) { return 1; }
    else { return 0; }
  }
  cout << "I shouldn't be here.\n";
  return 0;
}


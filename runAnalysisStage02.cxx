// TODO get linear discriminant intercept from a table instead of using fixed value

#include "InterfUtil.h" 
#include "analysisCuts.h" 
#include "TChain.h" 
#include "TTree.h" 
#include "TFile.h" 
#include "TTreeIndex.h" 
#include "TCanvas.h" 
#include "TStyle.h" 
#include "TROOT.h"
//#include "TApplication.h"
#include "TStyle.h" 
#include "TPad.h" 
#include "TText.h" 
#include "TDirectory.h" 
#include "TMarker.h"

#include "BedmapReader.h" 
#include "AnitaEventSummary.h"

#include "Adu5Pat.h" 
#include "UsefulAdu5Pat.h"

#include <healpix_base.h>

#include <sstream>

using namespace std;


vector<float> surfaceEllipseParms1(float phi, float theta, float dist, float errPh, float errTh,
        UsefulAdu5Pat* gps, BedmapReader* bedmap);

int canANITASeeStripe(double lon_anita,double lat_anita,double lat_satellite,double lon_satellite,double height_anita,double height_satellite);

pair<double,double> rotateOnSphere(double theta, double phi, double A, double B);

int main(int argc, char* argv[]) {
  printf("Hello world %s \n ", argv[0]);
  //TApplication* app = 0; app = new TApplication("app", &argc, argv);
  
  gROOT->SetBatch(kTRUE);
  aCutOrderStage2 = { // these cuts will be calculated and applied in this program
                                    //A_CUT_SATALLITE_AREA,
                                    A_CUT_CPOL_PEAK_SEPARATION, //only applied if final cuts are enabled
                                    A_CUT_CPOL_PEAK_STRENGTH, //only applied if final cuts are enabled -Jacob!
                                    A_CUT_PEAK_RATIO,
                                    A_CUT_CORR_PEAK,
                                    A_CUT_HILBERT_PEAK
                                    //A_CUT_DEADTIME,
                                    //A_CUT_SOLAR_DIRECT,
                                    //A_CUT_GEOSTAT_DIRECT,
                                    //A_CUT_GEOSTAT_REFL,
                                    //A_CUT_HP_BIN_REJECTED,
                                    //A_CUT_SATALLITE_AREA, //moved to start
                                    //A_CUT_CPOL_PEAK_SEPARATION, //only applied if final cuts are enabled
                                    //A_CUT_CPOL_PEAK_STRENGTH //only applied if final cuts are enabled -Jacob!
                                    //A_CUT_DIAGONAL
                                };

  aCutOrderStage1 = { // just for deciding wheteher to put in the cut table; stage 1 calculated these
                                    A_CUT_SOLAR_REFL,
                                    A_CUT_CONTINENT,
                                    A_CUT_THETA,
                                    A_CUT_L3_TRIG_DIR,
                                    A_CUT_WAIS,
                                    A_CUT_LDB
                             };
                      
  //aCutOrderStage2 = {}; // for testing and comparing to output from analyzereResultsIterator06
  //int aCutOrderStage2Limit = aCutOrderStage2.size();
  int startRunNum = 333;
  int endRunNum = 352;
  int pNum = 0;
  char pChars[3] = "HV";
  char pCharsC[3] = "LR";
  char fileDir[1024];
  char interfFileDir[1024];
  char runDesc[32];
  bool simulationMode = false;
  bool saveOutput = true;
  int maxEvents = 0;
  cPolPeakSeparationThreshold = 46.0;
  minCircPeakThreshold = 0.015;
  ldSlope = -6.0;
  ldInterceptThreshold = 0;
  bool doRebinning = false;           // set to true if you want events rebinned, will need to set offsets.
  bool offsetType = false;            // changes the way the offsets shift the healpix map, false is just by adding, which creates skewed bins.  true is by rotating the sphere about the x/y axis by the number of degrees corisponding to the offset, which does not create skews.
  bool applyFinalCuts = false;        // true to do final cuts  
  bool snagInterferometry = false;   // this should never need to be on
  bool use90Data = false;             // set to true to use the 90% anlysisOutput files.
  //char fileSuffix[1024] = "_h_90"; // suffix used for passingEvents file
  char fileSuffix[1024];
  int testEvent = 0; // set to 0 to turn this off, if not zero, it only looks at this event.
  char inputFileName[1024] = "";
  char whichPol = 'V'; //this is really just used to pick the hp orientation right now.
  float hpPhiOffset = 0;
  float hpThOffset = 0;

  for (int i=1; i<argc; ++i) {
    if (argv[i][0]=='-') {
      if (strlen(argv[i]) > 1) {
        if (argv[i][1] == 'D' || argv[i][1] == 'd') {
          if (strlen(argv[i]) > 2) {
            //char thisParm[8];
            unsigned int k=2;
            for (; k<=strlen(argv[i]); k++) {
              fileDir[k-2] = argv[i][k];
            }
            fileDir[k-2] = 0;
          }
        }
        if (argv[i][1] == 'P' || argv[i][1] == 'p') {
          if (strlen(argv[i]) > 2) {
            whichPol = toupper(argv[i][2]);
          }
        }
        if (argv[i][1] == 'F' || argv[i][1] == 'f') {
          if (strlen(argv[i]) > 2) {
            //char thisParm[8];
            unsigned int k=2;
            for (; k<=strlen(argv[i]); k++) {
              inputFileName[k-2] = argv[i][k];
            }
            inputFileName[k-2] = 0;
          }
        }
        if (argv[i][1] == 'S' || argv[i][1] == 's') {
          if (strlen(argv[i]) > 2) {
            //char thisParm[8];
            unsigned int k=2;
            for (; k<=strlen(argv[i]); k++) {
              fileSuffix[k-2] = argv[i][k];
            }
            fileSuffix[k-2] = 0;
          }
        }
        if (argv[i][1] == 'I' || argv[i][1] == 'd') {
          if (strlen(argv[i]) > 2) {
            //char thisParm[8];
            unsigned int k=2;
            for (; k<=strlen(argv[i]); k++) {
              interfFileDir[k-2] = argv[i][k];
            }
            interfFileDir[k-2] = 0;
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
        }
        else if (argv[i][1] == 'M' || argv[i][1] == 'm') {
          if (strlen(argv[i]) > 2) {
            char thisParm[8];
            unsigned int k=2;
            for (; k<=strlen(argv[i]); k++) {
              thisParm[k-2] = argv[i][k];
              //vbprintf("event %i hpBin %i \n", eventNumber, hpBin);
            }
            thisParm[k-2] = 0;
            maxEvents = atoi(thisParm);
          }
        }
        else if (argv[i][1] == 'V' || argv[i][1] == 'v') {verboseMode = true;}
        else if (argv[i][1] == 'U' || argv[i][1] == 'u') {simulationMode = true;}
        else if (argv[i][1] == 'B' || argv[i][1] == 'b') {doRebinning = true;}
        else if (argv[i][1] == 'A' || argv[i][1] == 'a') {applyFinalCuts = true;}
        else if (argv[i][1] == '9') {use90Data = true;}
      }
    } else {
      ++pNum;
      if (pNum==1) {startRunNum = atoi(argv[i]); endRunNum = startRunNum;}
      else if (pNum==2) {endRunNum = atoi(argv[i]);}
    }
  }
  printf("run number %i-%i \n", startRunNum, endRunNum);
  printf("file directory is %s \n", fileDir);
  printf("interferometry results file directory is %s \n", interfFileDir);
  char filename[2048];
  if( inputFileName == '\0' ) {
    if (simulationMode) {
      //sprintf(filename, "%s/analysisOutput_sim.root", fileDir);
      sprintf(filename, "%s/analysisOutput_801_808.root", fileDir);
    } else {
      sprintf(filename, "%s/analysisOutput_%03i_%03i.root", fileDir, startRunNum, endRunNum);
    }
    if (use90Data) {
      sprintf(filename, "%s/analysisOutput90/analysisOutput_%03i_%03i_90.root", fileDir, startRunNum, endRunNum);
    }
  } else {
    sprintf(filename, "%s/%s", fileDir, inputFileName);
  }
  
  printf("input file path is %s \n", filename);
  printf("current dir is %s \n", gDirectory->GetName());
  printf("input file suffix is %s \n", fileSuffix);
  
  printf("Cuts selected: \n");
  for (int k=0; k<aCutOrderStage2.size(); ++k) {
    printf(" %i %s \n", aCutOrderStage2[k], A_CUT_DESC[aCutOrderStage2[k]]);
  }
  printf("\n");

  //Oindree's new offsets
  if ( whichPol == 'V' ) {
    hpPhiOffset = 0.56;
    hpThOffset = -5.04;
  }
  if ( whichPol == 'H' ) {
    hpPhiOffset = 3.92;
    hpThOffset = -0.00;
  }
  /*
  //Jacob's Old offsets
  if ( whichPol == 'V' ) {
    hpPhiOffset = 3.36;
    hpThOffset = -3.36;
  }
  if ( whichPol == 'H' ) {
    hpPhiOffset = 0.00;
    hpThOffset = -4.48;
  }
  */
  printf("healpix shifts: %f %f \n",hpPhiOffset,hpThOffset);


  TFile* inFile = new TFile(filename);
  TTree* tree1 = 0;
  TTree* tree0 = 0;

  tree0 = (TTree*)inFile->Get("analysisOutput0");
  printf("input tree 0 has %lli entries \n", tree0->GetEntries());
  if (!doRebinning) {
    tree1 = (TTree*)inFile->Get("analysisOutput1");
    printf("input tree 1 has %lli entries \n", tree1->GetEntries());
  } else {
    tree1 = new TTree("analysisOutput1", "analysis1");
  }
   
  BedmapReader* bedmap = BedmapReader::Instance(false);
  Healpix_Ordering_Scheme scheme = Healpix_Ordering_Scheme::RING;
  T_Healpix_Base<int>* healpix = new T_Healpix_Base<int>(healPixN, scheme);

  TFile* coastGraphFile = new TFile("antarcticCoastGraph.root");
  TGraph* coastLineGr = 0;
  coastLineGr = (TGraph*) coastGraphFile->Get("Graph");
  printf("Antarctic coastline graph contains %i points \n", coastLineGr->GetN());
  int coastStride = 40;
  for (int k=1; k<coastLineGr->GetN()-1; ++k) {for (int j=0; j<coastStride; ++j) {coastLineGr->RemovePoint(k);}}
  coastLineGr->SetPoint(coastLineGr->GetN(), coastLineGr->GetX()[0], coastLineGr->GetY()[0]);
  for (int k=0; k<coastLineGr->GetN(); ++k) {coastLineGr->SetPoint(k, coastLineGr->GetX()[k]/1000, coastLineGr->GetY()[k]/1000);}
  printf("Antarctic coastline graph contains %i points \n", coastLineGr->GetN());
  coastGraphFile->Close();

  int eventNumber = 0;
  int runNumber =0;
  char pol = 0;
  float eventWeight = 0;
  float hPeak = 0;
  float cPeak = 0;
  float cohSnr = 0;
  float cohSnr1 = 0;
  float cohSnr2 = 0;
  float cohSnr3 = 0;
  float peakRatio = 0;
  float linPolFrac = 0;
  float cPolDist = 0;
  //float cPolMinDist = 0;
  float cPolPeak[2] = {0};
  float cPolLesserPeak = 0;
  float ldCutVal = 0;
  int hpBin = 0;
  float hpWeight = 0.;
  float minCircPeakVal = 0.;
  float deadTimeFrac = 0.;

  float aCutCounts[2][NUM_A_CUTS] = {{0.}}; // number of events that fail cut ("as first")
  float aCutOrderCounts[2][NUM_A_CUTS] = {{0.}}; // number of events failing as ordered cut ("as ordered")
  float aCutLastCounts[2][NUM_A_CUTS] = {{0.}}; // number of events failing only this cut
  float hpBinRejectedCount[2] = {0.};
  float ldRejectedCount[2] = {0.};
  float lowBinWeightCount[2] = {0.};
  //int aCutCountTotal[2] = {0}; // total events that fail any cut (as-first)
  float aCutOrderCountTotal[2] = {0.}; // total events failing any ordered cut ("as ordered")
  //int aCutLastCountTotal[2] = {0}; // total events failing any as-last cut
  float eventCount[2] = {0};
  int aCuts[2][NUM_A_CUTS] = {{0}};

  float aCutCounts5[2][NUM_A_CUTS] = {{0.}}; // number of events (snr > 5) that fail cut ("as first")
  float aCutOrderCounts5[2][NUM_A_CUTS] = {{0.}}; // number of events (snr > 5) failing as ordered cut ("as ordered")
  float aCutLastCounts5[2][NUM_A_CUTS] = {{0.}}; // number of events (snr > 5) failing only this cut
  float hpBinRejectedCount5[2] = {0.};
  float ldRejectedCount5[2] = {0.};
  float lowBinWeightCount5[2] = {0.};
  float aCutOrderCountTotal5[2] = {0.}; // total events (snr > 5) failing any ordered cut ("as ordered")
  float eventCount5[2] = {0};
  int aCuts5[2][NUM_A_CUTS] = {{0}};


  //float ellParms[4];
  float theta = 0; //event theta
  float phi = 0; //event phi
  float ea = 0;
  float no = 0;
  float lat = 0;
  float lon = 0;
  float alt = 0;
  Adu5Pat* gpsRaw = 0;
  UsefulAdu5Pat* gps = 0;

  float lPhi = 0; //event phi in l-pol

  map<int, float> hpEventsAccepted[2] = {map<int,float>()};
  map<int, float> hpEventsAcceptedPreFinalCuts[2] = {map<int,float>()};
  
  // store info about accepted events
  vector<int> accEventNumber[2] = {vector<int>(0)};
  vector<int> accEventBin[2] = {vector<int>(0)};
  vector<float> accEventWeight[2] = {vector<float>(0)};
  vector<float> accEventLat[2] = {vector<float>(0)};
  vector<float> accEventLon[2] = {vector<float>(0)};
  vector<float> accEventEa[2]= {vector<float>(0)};
  vector<float> accEventNo[2]= {vector<float>(0)};


  tree0->SetBranchAddress("runNumber", &runNumber);
  tree0->SetBranchAddress("eventNumber", &eventNumber);
  tree0->SetBranchAddress("pol", &pol);

  tree0->SetBranchAddress("eventWeight", &eventWeight);
  tree0->SetBranchAddress("hPeak", &hPeak);
  tree0->SetBranchAddress("cPeak", &cPeak);
  tree0->SetBranchAddress("cohSnr", &cohSnr);
  tree0->SetBranchAddress("cohSnr1", &cohSnr1);
  tree0->SetBranchAddress("cohSnr2", &cohSnr2);
  tree0->SetBranchAddress("cohSnr3", &cohSnr3);
  tree0->SetBranchAddress("peakRatio", &peakRatio);
  tree0->SetBranchAddress("linPolFrac", &linPolFrac);
  tree0->SetBranchAddress("cPolDist", &cPolDist);
  tree0->SetBranchAddress("cPolLesserPeak", &cPolLesserPeak);
  tree0->SetBranchAddress("cPolPeak", &cPolPeak); // reactivate when cPolPeak[] is put pack into AnalyzerResultsIterator output
  //tree0->SetBranchAddress("cPolMinDist", &cPolMinDist);
  tree0->SetBranchAddress("minCircPeakVal", &minCircPeakVal);
  tree0->SetBranchAddress("deadTimeFrac", &deadTimeFrac);
  //tree0->SetBranchAddress("aCuts", &aCuts); // reactivate when aCuts is put pack into AnalyzerResultsIterator output


  //tree0->SetBranchAddress("ellParms", &ellParms); // error ellipse parameters
  tree0->SetBranchAddress("theta", &theta); // event theta
  tree0->SetBranchAddress("phi", &phi); // event phi
  tree0->SetBranchAddress("ea", &ea); // event easting
  tree0->SetBranchAddress("no", &no); // event northing
  tree0->SetBranchAddress("lat", &lat); // event latitude
  tree0->SetBranchAddress("lon", &lon); // event longitude
  tree0->SetBranchAddress("alt", &alt); // event altitude
  tree0->SetBranchAddress("gps", &gps); // payload gps info
  tree0->SetBranchAddress("gpsRaw", &gpsRaw); // raw payload gps info

  tree0->SetBranchAddress("lPhi", &lPhi);
  
  if (!doRebinning) {
    tree1->SetBranchAddress("runNumber", &runNumber);
    tree1->SetBranchAddress("eventNumber", &eventNumber);
    tree1->SetBranchAddress("pol", &pol);
    tree1->SetBranchAddress("hpBin", &hpBin);
    tree1->SetBranchAddress("hpWeight", &hpWeight);
  }
  else {
    tree1->Branch("runNumber", &runNumber, "runNumber/I");
    tree1->Branch("eventNumber", &eventNumber, "eventNumber/I");
    tree1->Branch("pol", &pol, "pol/B");
    tree1->Branch("hpBin", &hpBin, "hpBin/I");
    tree1->Branch("hpWeight", &hpWeight, "hpWeight/F");
  }



  tree0->BuildIndex("eventNumber", "pol");
  TTreeIndex* tree0Index = (TTreeIndex*)tree0->GetTreeIndex();
  //tree1->BuildIndex("eventNumber", "pol");
  //TTreeIndex* tree1Index = (TTreeIndex*)tree1->GetTreeIndex();
  
  /////////////////////////////////////////////////////////
  //
  // Rebin events for new healpix orientations
  //
  // - Only do it if its set to on with a flag and earthTheta and earthPhi offsets are given as inputs
  // - If you do do it, loop through tree0Index->GetN()
  // -- Get error ellipse for each data point
  // -- Get the weights
  // -- Remake tree1
  // -- Remake tree1 index
  //
  ///////////////////////////////////////////////////////

  // new healpix for healpix stuff
  //Healpix_Ordering_Scheme scheme = Healpix_Ordering_Scheme::RING;
  //T_Healpix_Base<int>* healpix = new T_Healpix_Base<int>(healPixN, scheme);
  //lprintf("healpix map has %i pixels \n", healpix->Npix());

  float thError = 0.25; //this might depend on SNR later?
  float phiError = 0.50;
 
  //float hpPhiOffset = 0;
  //float hpThOffset = 0; this happens up above now
  //if (whichPol = 'V') {
  // hpPhiOffset = 3.36;
  // hpThOffset = -3.36;
  //}
  //if (whichPol = 'H') {
  // hpPhiOffset = 0.00;
  // hpThOffset = -4.48;
  //}


  if ( doRebinning )
  {
    vector<float> errEllipseParms1 = {0, 0, 0, 0};
    vector<vector<double> > errHexVertices = vector<vector<double> >(0);
    vector<vector<double> > errHexPixels = vector<vector<double> >(0);
    for (int e = 0; e<tree0Index->GetN(); ++e)
    {
      int entNum = tree0Index->GetIndex()[e];
      tree0->GetEntry(entNum);
      int getRes = tree0->GetEntryWithIndex(eventNumber, pol);
      if (testEvent != 0) {
        if (eventNumber != testEvent) {continue;}
      }
      bool useThisEntry = true;
      if (getRes <= 0) {lprintf("database error event %i pol %c - skipping \n", eventNumber, pol); useThisEntry = false;}
      //if (circPeakSepThreshold>0 && cPolDist>circPeakSepThreshold) {useThisEntry = false;} // i dont think we want these in here, they are done for final cuts only later.
      //if (circPeakStrengthThreshold>0 && minCircPeakVal<circPeakStrengthThreshold) {useThisEntry = false;}
      if (!useThisEntry) { continue; }
      vbprintf("Rebinning event %i \n", eventNumber);


      // now we have the event info so lets do the rebinning.
      // first get the error ellipse.

      UsefulAdu5Pat* gps = new UsefulAdu5Pat(gpsRaw);
      BedmapReader* bedmap = BedmapReader::Instance(false);

      double plEa = 0; double plNo = 0;
      bedmap->LonLattoEaNo(gps->longitude, gps->latitude, plEa, plNo);
      float plSrcPhi = atan2(no-plNo, ea-plEa);
      float plSrcDist = sqrt((no-plNo)*(no-plNo) + (ea-plEa)*(ea-plEa));
      errEllipseParms1 = surfaceEllipseParms1(phi, theta, plSrcDist, phiError, thError, gps, bedmap);

      // now turn those into pixels

      errHexVertices = hexagonVertices(errEllipseParms1[2], errEllipseParms1[3], errEllipseParms1[0], errEllipseParms1[1], plSrcPhi);
      //vbprintf(" Error hexagon vertices :"); for (int k=0; k<6; ++k) {vbprintf(" %f, %f \n", errHexVertices[0][k], errHexVertices[1][k$
      errHexPixels = ellipsePixels(errEllipseParms1[2], errEllipseParms1[3], errEllipseParms1[0], errEllipseParms1[1], plSrcPhi);


      map<int, float> hpWeights = map<int, float>();
      for (int k=0; k<errHexPixels[0].size(); ++k) {
        double thisLat, thisLon;
        bedmap->EaNoToLonLat(errHexPixels[0][k], errHexPixels[1][k], thisLon, thisLat);
        //vbprintf(" pixel ea/no = %f,%f \n", errHexPixels[0][k], errHexPixels[1][k]);
        //vbprintf(" pixel lat/lon = %f,%f \n", thisLat, thisLon);

        double thisTheta = 0;
        double thisPhi = 0;
        if (!offsetType) {
          thisTheta = (-thisLat + 90.0 + hpThOffset) * M_PI/180.0;
          thisPhi = (thisLon + hpPhiOffset) * M_PI/180.0;
        } else {
          pair<double,double> thisPoint = rotateOnSphere((90.0-thisLat)*M_PI/180.0, thisLon*M_PI/180.0, hpPhiOffset*M_PI/180.0, hpThOffset*M_PI/180.0);
          thisTheta = thisPoint.first;
          thisPhi = thisPoint.second;
        }
        while (thisPhi < -M_PI) thisPhi += M_PI;
        //vbprintf(" pixel theta/phi = %f,%f \n", thisTheta, thisPhi);
        pointing point(thisTheta, thisPhi);
        int pixelNo = healpix->ang2pix(point);


        float thisWeight = 0;
        if (simulationMode) {
          thisWeight = errHexPixels[2][k] * eventWeight;
        } else {
          thisWeight = errHexPixels[2][k]; // * eventWeight;
        }
        // I think eventWeight is just 1.0 for data, It is different for sim

        map<int, float>::iterator thisEntry = hpWeights.find(pixelNo);
        if (thisEntry == hpWeights.end()) {
          //vbprintf(" creating bin list entry for event %i %cPol bin %i %8.3f \n", eventNumber, pol, pixelNo, thisWeight);
          hpWeights.insert(pair<int, float>(pixelNo, thisWeight));
        } else {
          //vbprintf(" updating bin list entry for event %i %cPol bin %i %8.3f \n",
          // eventNumber, pol, pixelNo, thisWeight);
          hpWeights[pixelNo] += thisWeight;
        }
      }
      for (pair<int, float> thisEntry : hpWeights) {
        //runNumber = runNumber;
        hpBin = thisEntry.first;
        hpWeight = thisEntry.second;
        //eventNumber = eventNumber;
        //pol = pol;
        vbprintf(" writing bin list entry for event %i %cPol bin %i %8.3f\n",
                eventNumber, pol, hpBin, hpWeight);
        tree1->Fill();
      }
    }
  }
  tree1->BuildIndex("eventNumber", "pol");
  TTreeIndex* tree1Index = (TTreeIndex*)tree1->GetTreeIndex();

  /////////////////////////////////////////////////
  //
  // Done rebinning
  //
  ////////////////////////////////////////////////

  // get the intercept values by healpix bin
  map<int, float> hpBinInterceptMap;
  set<int> binsToKeep;
  if (applyFinalCuts) {
    char interceptFilename[1024];
    //sprintf(interceptFilename, "%s/intercept_no_cpol.txt", fileDir);
    sprintf(interceptFilename, "%s/intercept%s.txt", fileDir, fileSuffix);
    printf("intercept filename is %s \n", interceptFilename);
    ifstream interceptFile(interceptFilename);
    string interceptLine;
    while (!interceptFile.eof()) {
    getline(interceptFile, interceptLine);
      printf("%s \n", interceptLine.c_str());
      if (interceptLine.c_str()[0] != '#') {
        float dummyVar;
        int intcptHpBin;
        int intcptBinStatus;
        int intcptBinBinStatus;
        int intcpt5Bins;
        float intcptIntercept;
        float intcptPVal;
        int intcptPValOk;
        float intcptBg;
        stringstream interceptStr(interceptLine);
        interceptStr >> intcptHpBin >> intcptBinStatus >> intcptBinBinStatus >> dummyVar >> dummyVar >> dummyVar >> intcpt5Bins >> intcptIntercept >> intcptPVal >> 
intcptPValOk >> dummyVar >> intcptBg;
        printf("bin #%i status=%i >5Bins=%i intercept=%f pVal=%f pValOk=%i bg=%f \n",
                intcptHpBin, intcptBinStatus, intcpt5Bins, intcptIntercept, intcptPVal, intcptPValOk, intcptBg);
        if (intcptBinStatus==0 && intcptBinBinStatus==0) {
          hpBinInterceptMap.insert(pair<int, float>(intcptHpBin, intcptIntercept));
        }
      }
    }
    printf("linear discriminant intercept table: \n");
    for (pair<int, float> thisEntry : hpBinInterceptMap) {
      printf("bin %i intercept %f \n", thisEntry.first, thisEntry.second);
    }
    char binKeepFilename[1024];
    //sprintf(binKeepFilename, "%s/binsOver0.01_no_cpol.txt", fileDir);
    sprintf(binKeepFilename, "%s/binsOver0.01%s.txt", fileDir, fileSuffix);
    printf("bins-to-keep filename is %s \n", binKeepFilename);
    ifstream binKeepFile(binKeepFilename);
    while (!binKeepFile.eof()) {
      int thisBin;
      binKeepFile >> thisBin;
      binsToKeep.insert(thisBin);
    }
    printf("bins over 0.01 accumulated fraction of sim events passing\n ");
    for (int thisBin : binsToKeep) {
      printf(" %4i \n", thisBin);
    }
  }
  //return 0; // for testing read of hpbin y-intercepts
  
  char outputFilename[1024];

  if (!use90Data) {
    if (simulationMode) {
      if (applyFinalCuts) {
        sprintf(outputFilename, "%s/analysisOutput_2_finalCuts_sim.root", fileDir);
      } else {
        sprintf(outputFilename, "%s/analysisOutput_2_sim.root", fileDir);
      }
    } else {
      if (applyFinalCuts) {
        sprintf(outputFilename, "%s/analysisOutput_2_finalCuts_%03i_%03i.root", fileDir, startRunNum, endRunNum);
      } else {
        sprintf(outputFilename, "%s/analysisOutput_2_%03i_%03i.root", fileDir, startRunNum, endRunNum);
      }
    }
  } else {
    if (simulationMode) {
      if (applyFinalCuts) {
        sprintf(outputFilename, "%s/analysisOutput_2_90_finalCuts_sim_%c.root", fileDir, whichPol);
      } else {
        sprintf(outputFilename, "%s/analysisOutput_2_90_sim.root", fileDir);
      }
    } else {
      if (applyFinalCuts) {
        sprintf(outputFilename, "%s/analysisOutput_2_90_finalCuts_%03i_%03i_%c.root", fileDir, startRunNum, endRunNum, whichPol);
      } else {
        sprintf(outputFilename, "%s/analysisOutput_2_90_%03i_%03i.root", fileDir, startRunNum, endRunNum);
      }
    }
  }
    
  TFile* outputFile = new TFile(outputFilename, "RECREATE");

  TTree* outputTree0 = tree0->CloneTree(0);
  TTree* outputTree1 = tree1->CloneTree(0);

  int prevEvent = 0;
  char prevPol = 0;
  printf("cPol peak separation threshold is %f \n", cPolPeakSeparationThreshold);
  
  // open file for saving acceteped events
  sprintf(filename, "%s/acceptedEventsList_%03i_%03i%s.txt", fileDir, startRunNum, endRunNum, fileSuffix);
  FILE* acceptedEventFile = fopen(filename, "w");
  
  TH2F* circPeakSepLDHist[2] = {0};
  TH2F* circPeakSepSNRHist[2] = {0};
  TH2F* corrHilbHist[2] = {0};
  TH2F* corrSnrHist[2] = {0};
  TH2F* lphiLonHist[2] = {0};
  TH1F* minCircValHist[2] = {0};
  TH2F* minCircValLDHist[2] = {0};
  TH1F* peakRatioHist[2] = {0};
  TH1F* snrHist0[2] = {0};
  TH1F* snrHist1[2] = {0};
  TH1F* snrHist2[2] = {0};
  TH1F* snrHist3[2] = {0};
  TH2F* hpCountsHist[2] = {0};
  map<int,TH1F*> diffHistMap[2];
  
  TH1F* peakSepCircHist = 0; // separation (angle) between Lpol and RPol peaks; satellite identification
  TH2F* corrPeakPolHist = 0;
  TH2F* cPolPeakDistHist = 0; // lesser of RPol, LPol peak vs RPol,LPol peak separation
  TH2F* circPeakHist = 0; // LPol, RPol correlation peak

  char name[128];
  char title[128];
  // instantiate single histograms
  sprintf(title, "CPol peak value: run %i-%i, %s; %cpol; %cpol",
          startRunNum, endRunNum, runDesc, pCharsC[0], pCharsC[1]);
  sprintf(name, "circPeakHist");
  circPeakHist = new TH2F(name, title, 100, 0, 0.3, 100, 0, 0.3);

  sprintf(title, "Correlation peak by pol: run %i-%i, %s; %cpol; %cpol",
          startRunNum, endRunNum, runDesc, pChars[0], pChars[1]);
  sprintf(name, "corrPeakPolHist");
  corrPeakPolHist = new TH2F(name, title, 100, 0, 0.6, 100, 0, 0.6);

  sprintf(title, "Lesser CPol peak value: run %i-%i, %s; lesser cpol corr peak; cpol peak separation ",
          startRunNum, endRunNum, runDesc);
  sprintf(name, "cPolPeakDistHist");
  cPolPeakDistHist = new TH2F(name, title, 100, 0, 0.6, 180, 0, 180);

  sprintf(title, "CPol peak separation: run %i-%i, %s",
          startRunNum, endRunNum, runDesc);
  sprintf(name, "peakSepCircHist");
  peakSepCircHist = new TH1F(name, title, 180, 0, 180);

  // instantiate 2-pol histograms
  for (int p=0; p<2; ++p) {
    sprintf(title, "circPol peak separation vs SNR_hilb: run %i-%i, %cPol %s; SNR_hilb; peakSep",
            startRunNum, endRunNum, pChars[p], runDesc);
    sprintf(name, "circPeakSepLDHist%1i", p);
    circPeakSepLDHist[p] = new TH2F(name, title, 90, 0, 30, 180, 0, 180);
 
    sprintf(title, "circPol peak separation vs SNR: run %i-%i, %cPol %s; SNR; peakSep",
            startRunNum, endRunNum, pChars[p], runDesc);
    sprintf(name, "circPeakSepSNRHist%1i", p);
    circPeakSepSNRHist[p] = new TH2F(name, title, 100, 0, 20, 180, 0, 180);

    sprintf(title, "Correlation peak vs Hilbert peak: run %i-%i, %cPol %s;hilb;corr",
            startRunNum, endRunNum, pChars[p], runDesc);
    sprintf(name, "corrHilbHist%1i", p);
    corrHilbHist[p] = new TH2F(name, title, 100, 0, 150, 100, 0, 0.4);

    sprintf(title, "Correlation peak vs SNR: run %i-%i, %cPol %s;SNR;corr",
            startRunNum, endRunNum, pChars[p], runDesc);
    sprintf(name, "corrSnrHist%1i", p);
    corrSnrHist[p] = new TH2F(name, title, 100, 0, 15, 100, 0, 0.4);

    sprintf(title, "LCP Phi vs Payload Longitude: run %i-%i, %cPol %s;lphi;longitude",
            startRunNum, endRunNum, pChars[p], runDesc);
    sprintf(name, "lphiLonHist%1i", p);
    lphiLonHist[p] = new TH2F(name, title, 100, 0.0, 360.0, 100, -180.0, 180.0);

    sprintf(title, "Lesser CPol correlation at %cPol peak location : run %i-%i, %s; lesser cpol value ",
            pChars[p], startRunNum, endRunNum, runDesc);
    sprintf(name, "minCircValHist%1i", p);
    minCircValHist[p] = new TH1F(name, title, 100, 0, 0.5);

    sprintf(title, "Lesser CPol correlation at %cPol peak location vs SNR_hilb: run %i-%i, %s;SNR_hilb;lesser cpol value ",
            pChars[p], startRunNum, endRunNum, runDesc);
    sprintf(name, "minCircValLDHist%1i", p);
    minCircValLDHist[p] = new TH2F(name, title, 100, 0, 30, 100, 0, 0.25);

    sprintf(title, "2nd/1st peak ratio: run %i-%i, %cPol %s",
            startRunNum, endRunNum, pChars[p], runDesc);
    sprintf(name, "peakRatioHist%1i", p);
    peakRatioHist[p] = new TH1F(name, title, 120, 0, 1.2);

    sprintf(title, "Coherent reconstruction SNR (method 0): run %i-%i, %cPol %s",
            startRunNum, endRunNum, pChars[p], runDesc);
    sprintf(name, "snrHist0%1i", p);
    snrHist0[p] = new TH1F(name, title, 175, 0, 35);

    sprintf(title, "Coherent reconstruction SNR (method 1): run %i-%i, %cPol %s",
            startRunNum, endRunNum, pChars[p], runDesc);
    sprintf(name, "snrHist1%1i", p);
    snrHist1[p] = new TH1F(name, title, 175, 0, 35);

    sprintf(title, "Coherent reconstruction SNR (method 2): run %i-%i, %cPol %s",
            startRunNum, endRunNum, pChars[p], runDesc);
    sprintf(name, "snrHist2%1i", p);
    snrHist2[p] = new TH1F(name, title, 175, 0, 35);

    sprintf(title, "Coherent reconstruction SNR (method 3): run %i-%i, %cPol %s",
            startRunNum, endRunNum, pChars[p], runDesc);
    sprintf(name, "snrHist3%1i", p);
    snrHist3[p] = new TH1F(name, title, 175, 0, 35);
        
    sprintf(title, "Binned Continent Weighted Events: run %i-%i, %cPol %s;ea(km);no(km) ",
            startRunNum, endRunNum, pChars[p], runDesc);
    sprintf(name, "hpCountHist%1i", p);
    hpCountsHist[p] = new TH2F(name, title, 600, -3000, 3000, 600, -3000, 3000);

    // special settings for setting up diff plot hists for each healpix bin in binsToKeep and hpBinInterceptMap
    for (int thisBin : binsToKeep) {
      map<int, float>::iterator thisHpBinEntry = hpBinInterceptMap.find(thisBin);
      if (thisHpBinEntry != hpBinInterceptMap.end()) {
        sprintf(title, "LD Histogram of data: bin %i, run %i-%i, %cPol %s",
                thisBin, startRunNum, endRunNum, pChars[p], runDesc);
        sprintf(name, "diffPlot%1i0%i", p, thisBin);
        TH1F* tempHist = new TH1F(name, title, 50, 0, 25);
        diffHistMap[p].insert(pair<int,TH1F*>(thisBin, tempHist));
      }
    }
        
  }
  
  bool keepIt;
  bool kept = false;
  if (maxEvents == 0) {maxEvents = tree1Index->GetN();}
  cout << "maxEvents = " << maxEvents << endl;
  cout << "maxEvents0 = " << tree0Index->GetN() << endl;
  for (int e=0; e<maxEvents; ++e) {
  //for (int e=0; e<1000; ++e) {
    int entNum = tree1Index->GetIndex()[e];
    tree1->GetEntry(entNum);
    // look up the bin in the intercept table and ignore this entry if bin not found
    map<int, float>::iterator thisHpBinEntry = hpBinInterceptMap.find(hpBin);
    ldInterceptThreshold = 0;
    if (thisHpBinEntry != hpBinInterceptMap.end()) {
      if (binsToKeep.find(hpBin) != binsToKeep.end()) {
        ldInterceptThreshold = hpBinInterceptMap[hpBin];
      }
    }
    //printf("event %i hpBin %i \n", eventNumber, hpBin);
    if (testEvent != 0) {
      if (eventNumber != testEvent) {
        continue; }
      cout << "Event number: " << testEvent << " in hpBin: " << hpBin << endl;
     
    }
    vbprintf("index number %i \n", e);
    vbprintf("entry number %i \n", entNum);
    bool newEventPol = (eventNumber != prevEvent || pol != prevPol);
    bool newEvent = (eventNumber != prevEvent);
    if (newEvent) {
      kept = false;
    }
    int p;
    if (newEventPol) {
      keepIt = true;
      p = -1;
      int getRes = tree0->GetEntryWithIndex(eventNumber, pol);
      if (pol == 'H') {p=0;} else if (pol == 'V') {p=1;}
      if (getRes <= 0) {printf("database error - skipping this Healpix bin \n"); continue;}
      //printf(" event number %i %cPol \n", eventNumber, pol) ;
      //printf(" indexed tree get result is %i \n", getRes);
      //cutVal = -cutSlope*cPeak + cohSnr2;
      //if (hpBin == ourBin) {vbprintf(" event number %i %cPol, bin=%i, cPeak=%.4f, hPeak=%.1f, cutval=%.4f \n", eventNumber, pol, hpBin, cPeak, hPeak, cutVal);}
      
      if (!simulationMode) { eventCount[p]+=eventWeight; } //should this be hp weight too? idk?
      else { eventCount[p] += eventWeight; }
      if (simulationMode && cohSnr2 > 5) { eventCount5[p] += eventWeight; }
      prevEvent = eventNumber;
      // unnormalize min cPol peak window value
      minCircPeakVal *= cPeak;
      prevPol = pol;
      vbprintf("event number %i pol=%c cPolDist=%f cPolLesserPeak=%f \n", eventNumber, pol, cPolDist, cPolLesserPeak);
      
      // apply the smart cuts: forget whatever stage 1 said
      // and reassign the flag if the cut is in the cut order.
      // otherwise, just clear the flag

      // We are goign to try using the R pol cPeak here instead of H or V
      //ldCutVal = -2*ldSlope * cPolPeak[1] + cohSnr2; //factor of two because we expect H and V to be 1/2 as strong in R
      ldCutVal = -ldSlope * cPeak + cohSnr2;
      //ldCutVal = -ldSlope * cPeak + cohSnr3;

      //cout << ldCutVal << endl;

      aCuts[p][A_CUT_PEAK_RATIO] = 0;
      if (simulationMode && cohSnr2 > 5.0) {
        aCuts5[p][A_CUT_PEAK_RATIO] = 0;
      }
      if ((find (aCutOrderStage2.begin(), aCutOrderStage2.end(), A_CUT_PEAK_RATIO) !=aCutOrderStage2.end()) && (peakRatio > peakRatioThreshold)) {
        vbprintf(" 2nd/1st peak ratio is over threshold: cutting this event\n");
        aCutCounts[p][A_CUT_PEAK_RATIO] += eventWeight;
        aCuts[p][A_CUT_PEAK_RATIO] = 1;
        if (simulationMode && cohSnr2 > 5.0) {
          aCutCounts5[p][A_CUT_PEAK_RATIO] += eventWeight;
          aCuts5[p][A_CUT_PEAK_RATIO] = 1;
        }
        keepIt = false;
      }
      aCuts[p][A_CUT_CORR_PEAK] = 0;
      if (simulationMode && cohSnr2 > 5.0) {
        aCuts5[p][A_CUT_CORR_PEAK] = 0;
      }
      if ((find (aCutOrderStage2.begin(), aCutOrderStage2.end(), A_CUT_CORR_PEAK) !=aCutOrderStage2.end()) && (cPeak < corrPeakThresh)) {
        vbprintf(" correlation peak is under threshold: cutting this event\n");
        aCutCounts[p][A_CUT_CORR_PEAK] += eventWeight;
        aCuts[p][A_CUT_CORR_PEAK] = 1;
        if (simulationMode && cohSnr2 > 5.0) {
          aCutCounts5[p][A_CUT_CORR_PEAK] += eventWeight;
          aCuts5[p][A_CUT_CORR_PEAK] = 1;
        }
        keepIt = false;
      }
      aCuts[p][A_CUT_HILBERT_PEAK] = 0;
      if (simulationMode && cohSnr2 > 5.0) {
        aCuts5[p][A_CUT_HILBERT_PEAK] = 0;
      }
      if ((find (aCutOrderStage2.begin(), aCutOrderStage2.end(), A_CUT_HILBERT_PEAK) !=aCutOrderStage2.end()) && (hPeak < hilbPeakThresh)) {
        vbprintf(" hilbert peak is under threshold: cutting this event\n");
        aCutCounts[p][A_CUT_HILBERT_PEAK] += eventWeight;
        aCuts[p][A_CUT_HILBERT_PEAK] = 1;
        if (simulationMode && cohSnr2 > 5.0) {
          aCutCounts5[p][A_CUT_HILBERT_PEAK] += eventWeight;
          aCuts5[p][A_CUT_HILBERT_PEAK] = 1;
        }
        keepIt = false;
      }
      aCuts[p][A_CUT_LINPOL_FRAC] = 0;
      if (simulationMode && cohSnr2 > 5.0) {
        aCuts5[p][A_CUT_LINPOL_FRAC] = 0;
      }
      if ((find (aCutOrderStage2.begin(), aCutOrderStage2.end(), A_CUT_LINPOL_FRAC) !=aCutOrderStage2.end()) && (linPolFrac < minLinPolFrac)) {
        vbprintf(" linear polarization fraction is under threshold: cutting this event\n");
        aCutCounts[p][A_CUT_LINPOL_FRAC] += eventWeight;
        aCuts[p][A_CUT_LINPOL_FRAC] = 1;
        if (simulationMode && cohSnr2 > 5.0) {
          aCutCounts5[p][A_CUT_LINPOL_FRAC] += eventWeight;
          aCuts5[p][A_CUT_LINPOL_FRAC] = 1;
        }
        keepIt = false;
      }
      aCuts[p][A_CUT_DEADTIME] = 0;
      if (simulationMode && cohSnr2 > 5.0) {
        aCuts5[p][A_CUT_DEADTIME] = 0;
      }
      if ((find (aCutOrderStage2.begin(), aCutOrderStage2.end(), A_CUT_DEADTIME) !=aCutOrderStage2.end()) && (deadTimeFrac > deadTimeThreshold)) {
        vbprintf(" deadtime fraction is over threshold: cutting this event\n");
        aCutCounts[p][A_CUT_DEADTIME] += eventWeight;
        aCuts[p][A_CUT_DEADTIME] = 1;
        if (simulationMode && cohSnr2 > 5.0) {
          aCutCounts5[p][A_CUT_DEADTIME] += eventWeight;
          aCuts5[p][A_CUT_DEADTIME] = 1;
        }
        keepIt = false;
      }

      //if (true) {
      if (applyFinalCuts) {
        aCuts[p][A_CUT_CPOL_PEAK_SEPARATION] = 0;
        if (simulationMode && cohSnr2 > 5.0) {
          aCuts5[p][A_CUT_CPOL_PEAK_SEPARATION] = 0;
        }
        if ((find (aCutOrderStage2.begin(), aCutOrderStage2.end(), A_CUT_CPOL_PEAK_SEPARATION) !=aCutOrderStage2.end()) && (cPolDist > cPolPeakSeparationThreshold)) {
          vbprintf(" cPol separation %f is over threshold: cutting this event\n", cPolPeakSeparationThreshold);
          aCutCounts[p][A_CUT_CPOL_PEAK_SEPARATION] += eventWeight;
          aCuts[p][A_CUT_CPOL_PEAK_SEPARATION] = 1;
          if (simulationMode && cohSnr2 > 5.0) {
            aCutCounts5[p][A_CUT_CPOL_PEAK_SEPARATION] += eventWeight;
            aCuts5[p][A_CUT_CPOL_PEAK_SEPARATION] = 1;
          }
          keepIt = false;
        }
        aCuts[p][A_CUT_CPOL_PEAK_STRENGTH] = 0;
        if (simulationMode && cohSnr2 > 5.0) {
          aCuts5[p][A_CUT_CPOL_PEAK_STRENGTH] = 0;
        }
        if ((find (aCutOrderStage2.begin(), aCutOrderStage2.end(), A_CUT_CPOL_PEAK_STRENGTH) !=aCutOrderStage2.end()) && (minCircPeakVal < minCircPeakThreshold)) {
          vbprintf(" minimum cPol value in linPol peak window %f is under threshold: cutting this event\n", minCircPeakVal);
          aCutCounts[p][A_CUT_CPOL_PEAK_STRENGTH] += eventWeight;
          aCuts[p][A_CUT_CPOL_PEAK_STRENGTH] = 1;
          if (simulationMode && cohSnr2 > 5.0) {
            aCutCounts5[p][A_CUT_CPOL_PEAK_STRENGTH] += eventWeight;
            aCuts5[p][A_CUT_CPOL_PEAK_STRENGTH] = 1;
          }
          keepIt = false;
        }
      }
      
      // calculate payload lon - event L-pol phi here before cut

      // correct for heading
      float modLPhi = fmod((lPhi - gps->heading + 360), 360);

      float mapPos = gps->longitude - modLPhi;
      while (mapPos > 160) { mapPos -= 360; } //skewed limits because of our logic checks.
      while (mapPos < -200) { mapPos += 360; }
      float plLat = gps->latitude;
      float plLon = gps->longitude;
      float plAlt = (gps->altitude); //in m?
      float sAlt = 35786000; //in m
      float sLat = -10.0;
      aCuts[p][A_CUT_SATALLITE_AREA] = 0;
      if (simulationMode && cohSnr2 > 5.0) {
        aCuts5[p][A_CUT_SATALLITE_AREA] = 0;
      }
      /*
      if ( (find (aCutOrderStage2.begin(), aCutOrderStage2.end(), A_CUT_SATALLITE_AREA) != aCutOrderStage2.end())
          && ( (-165.115 > mapPos && mapPos > -197.715 && canANITASeeStripe(plLon, plLat, sLat, -182.615, plAlt, sAlt) )
              || ( -77.035 > mapPos && mapPos > -120.535 && canANITASeeStripe(plLon, plLat, sLat, -100.435, plAlt, sAlt) ) ( 3.465 > mapPos && mapPos > -48.235 && 
              || canANITASeeStripe(plLon, plLat, sLat, -20.035, plAlt, sAlt) ) ( 142.4225 > mapPos && mapPos > 39.485 && canANITASeeStripe(plLon, plLat, sLat, 90.95, plAlt, 
              || sAlt) )
              //|| ( 112.65 > mapPos && mapPos > 60.97 && canANITASeeStripe(plLon, plLat, sLat, 86.815, plAlt, sAlt) )
             )
          && (cPolPeak[0] * 1.1 > cPolPeak[1] ) // if lPeak*1.1 > rPeak then cut the event
          )
      */
      if ( (find (aCutOrderStage2.begin(), aCutOrderStage2.end(), A_CUT_SATALLITE_AREA) != aCutOrderStage2.end())
          && ( (-193.4150 < mapPos && mapPos < -171.8150 && canANITASeeStripe(plLon, plLat, sLat, -182.6150, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 1.7 )
              || (-106.6350 < mapPos && mapPos < -89.4350 && canANITASeeStripe(plLon, plLat, sLat, -91.0350, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 2.2 ) 
              || ( -46.3350 < mapPos && mapPos < -33.4350 && canANITASeeStripe(plLon, plLat, sLat, -39.8850, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 2.2 ) 
              || (  59.3100 < mapPos && mapPos <  72.1100 && canANITASeeStripe(plLon, plLat, sLat,  65.7100, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 1.7 ) 
              || (  52.4950 < mapPos && mapPos <  58.9850 && canANITASeeStripe(plLon, plLat, sLat,  55.7400, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 2.0 ) 
              || (  69.4850 < mapPos && mapPos <  87.6850 && canANITASeeStripe(plLon, plLat, sLat,  78.5850, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 2.0 ) 
              || ( 104.3225 < mapPos && mapPos < 109.8225 && canANITASeeStripe(plLon, plLat, sLat, 107.0725, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 2.0 ) 
              || ( 124.5225 < mapPos && mapPos < 134.1225 && canANITASeeStripe(plLon, plLat, sLat, 129.3225, plAlt, sAlt) && cPolPeak[1]/cPolPeak[0] < 2.0 )
             )
          //&& (cPolPeak[0] * 1.1 > cPolPeak[1] ) // if lPeak*1.1 > rPeak then cut the event
         )
      {
        //cout << " mapPos: " << mapPos << endl;
        //cout << " gps: (" << plLon <<", " << plLat << ", " << plAlt << ")" << endl;
        //cout << " sat: (???, " << sLat << ", " << sAlt << ")" << endl;
        vbprintf(" event is within satallite stripe: cutting this event, plLon=%f; lPolPhi=%f\n",gps->longitude,modLPhi);
        //sleep(1);
        aCutCounts[p][A_CUT_SATALLITE_AREA] += eventWeight;
        aCuts[p][A_CUT_SATALLITE_AREA] = 1;
        if (simulationMode && cohSnr2 > 5.0) {
          aCutCounts5[p][A_CUT_SATALLITE_AREA] += eventWeight;
          aCuts5[p][A_CUT_SATALLITE_AREA] = 1;
        }
        keepIt = false;
      }
        

      //this is just for testing.  Its a weird version of the ldCut thats not bin dependent
      aCuts[p][A_CUT_DIAGONAL] = 0;
      //if ((find (aCutOrderStage2.begin(), aCutOrderStage2.end(), A_CUT_DIAGONAL) !=aCutOrderStage2.end()) && (ldCutVal < ldInterceptThreshold)) {
      if ((find (aCutOrderStage2.begin(), aCutOrderStage2.end(), A_CUT_DIAGONAL) !=aCutOrderStage2.end()) && (ldCutVal < 4.0)) {
        vbprintf(" linear discriminant value %f is under threshold: cutting this event\n", ldCutVal);
        aCutCounts[p][A_CUT_DIAGONAL] += eventWeight;
        aCuts[p][A_CUT_DIAGONAL] = 1;
        keepIt = false;
      }

      // count the events failing as-ordered cuts
      // iterate through both cut lists until a fail is reached, then increment the counter
      bool failed = false;
      for (int k=0; k<aCutOrderStage1.size() && !failed; ++k) {
        if (aCuts[p][aCutOrderStage1[k]] == 1) {
          aCutOrderCounts[p][aCutOrderStage1[k]]+=eventWeight;
          aCutOrderCountTotal[p]+=eventWeight;
          failed = true;
        }
      }
      for (int k=0; k<aCutOrderStage2.size() && !failed; ++k) {
        if (aCuts[p][aCutOrderStage2[k]] == 1) {
          aCutOrderCounts[p][aCutOrderStage2[k]]+=eventWeight;
          aCutOrderCountTotal[p]+=eventWeight;
          failed = true;
        }
      }
      // do the same thing for the sim events of > 5 snr
      if (simulationMode && cohSnr2 > 5.0) {
        bool failed5 = false;
        for (int k=0; k<aCutOrderStage1.size() && !failed5; ++k) {
          if (aCuts5[p][aCutOrderStage1[k]] == 1) {
            aCutOrderCounts5[p][aCutOrderStage1[k]]+=eventWeight;
            aCutOrderCountTotal5[p]+=eventWeight;
            failed5 = true;
          }
        }
        for (int k=0; k<aCutOrderStage2.size() && !failed5; ++k) {
          if (aCuts5[p][aCutOrderStage2[k]] == 1) {
            aCutOrderCounts5[p][aCutOrderStage2[k]]+=eventWeight;
            aCutOrderCountTotal5[p]+=eventWeight;
            failed5 = true;
          }
        }
      }

      // count the events failing as-last cuts
      int numFails = 0;
      // iterate through the cut list, if only one cut was failed, count it
      for (int k=0; k<aCutOrderStage1.size(); ++k) {if (aCuts[p][aCutOrderStage1[k]] == 1) {++numFails;}}
      for (int k=0; k<aCutOrderStage2.size(); ++k) {if (aCuts[p][aCutOrderStage2[k]] == 1) {++numFails;}}
      if (numFails == 1) {
        for (int k=0; k<aCutOrderStage1.size(); ++k) {if (aCuts[p][aCutOrderStage1[k]] == 1) {aCutLastCounts[p][aCutOrderStage1[k]]+=eventWeight;}}
        for (int k=0; k<aCutOrderStage2.size(); ++k) {if (aCuts[p][aCutOrderStage2[k]] == 1) {aCutLastCounts[p][aCutOrderStage2[k]]+=eventWeight;}}
      }
      // do the same thing for sim events with snr > 5
      if ( simulationMode && cohSnr2 > 5.0 ) {
        int numFails5 = 0;
        // iterate through the cut list, if only one cut was failed, count it
        for (int k=0; k<aCutOrderStage1.size(); ++k) {if (aCuts5[p][aCutOrderStage1[k]] == 1) {++numFails5;}}
        for (int k=0; k<aCutOrderStage2.size(); ++k) {if (aCuts5[p][aCutOrderStage2[k]] == 1) {++numFails5;}}
        if (numFails5 == 1) {
          for (int k=0; k<aCutOrderStage1.size(); ++k) {if (aCuts5[p][aCutOrderStage1[k]] == 1) {aCutLastCounts5[p][aCutOrderStage1[k]]+=eventWeight;}}
          for (int k=0; k<aCutOrderStage2.size(); ++k) {if (aCuts5[p][aCutOrderStage2[k]] == 1) {aCutLastCounts5[p][aCutOrderStage2[k]]+=eventWeight;}}
        }
      }

      if (keepIt) {
        vbprintf(" writing event master %i weight is %f\n", eventNumber, eventWeight);
        vbprintf(" minCircPeakVal = %f \n", minCircPeakVal);
        outputTree0->Fill();
        
        // fill the histograms
        if (!kept) {
          vbprintf("filling dual-pol histograms\n");
          vbprintf(" cPolPeak[0]=%f, cPolPeak[1]=%f eventWeight=%f \n", cPolPeak[0], cPolPeak[1], eventWeight);
          vbprintf(" cPolLesserPeak=%f, cPolDist=%f \n", cPolLesserPeak, cPolDist);
          
          circPeakHist->Fill(cPolPeak[0], cPolPeak[1], eventWeight); // the event weight should be the same for both pols
          cPolPeakDistHist->Fill(cPolLesserPeak, cPolDist, eventWeight);
          peakSepCircHist->Fill(cPolDist, eventWeight);
          corrPeakPolHist->Fill(cPolPeak[0], cPolPeak[1], eventWeight);
        }      
        circPeakSepSNRHist[p]->Fill(cohSnr2, cPolDist, eventWeight);
        circPeakSepLDHist[p]->Fill(ldCutVal, cPolDist, eventWeight);
        corrHilbHist[p]->Fill(hPeak, cPeak, eventWeight);
        corrSnrHist[p]->Fill(cohSnr2, cPeak, eventWeight);
        minCircValHist[p]->Fill(minCircPeakVal, eventWeight);
        minCircValLDHist[p]->Fill(ldCutVal, minCircPeakVal, eventWeight);
        peakRatioHist[p]->Fill(peakRatio, eventWeight);
        vbprintf("snr0=%f, weight=%f polNum=%i \n", cohSnr, eventWeight, p);
        snrHist0[p]->Fill(cohSnr, eventWeight);
        snrHist1[p]->Fill(cohSnr1, eventWeight);
        snrHist2[p]->Fill(cohSnr2, eventWeight);
        snrHist3[p]->Fill(cohSnr3, eventWeight);
        lphiLonHist[p]->Fill( fmod((lPhi-gps->heading + 360), 360) , gps->longitude, eventWeight);
        kept = true;
      }      
    }
    if (keepIt) {

      if (hpEventsAcceptedPreFinalCuts[p].find(hpBin) == hpEventsAcceptedPreFinalCuts[p].end()) {
        hpEventsAcceptedPreFinalCuts[p].insert(pair<int, float>(hpBin, hpWeight));
      } else {
        hpEventsAcceptedPreFinalCuts[p][hpBin] += hpWeight;
      }


      // do the final cuts
      bool keepDetail = true;

      if (applyFinalCuts) {
        //if (binsToKeep.find(hpBin) != binsToKeep.end()) { // healpix bin was rejected
        if (ldInterceptThreshold == 0) {
          keepDetail = false;
          hpBinRejectedCount[p] += hpWeight;
          aCutOrderCountTotal[p] += hpWeight;
          if (simulationMode && cohSnr2 > 5.0 ) {
            hpBinRejectedCount5[p] += hpWeight;
            aCutOrderCountTotal5[p] += hpWeight;
          }
        } else if (ldCutVal < ldInterceptThreshold) {
          keepDetail = false;
          ldRejectedCount[p] += hpWeight;
          aCutOrderCountTotal[p] += hpWeight;
          if (simulationMode && cohSnr2 > 5.0 ) {
            ldRejectedCount5[p] += hpWeight;
            aCutOrderCountTotal5[p] += hpWeight;
          }
        } else if (  (!simulationMode && hpWeight < 0.5)
                   ||( simulationMode && hpWeight/eventWeight < 0.5) ) {
          keepDetail = false;
          lowBinWeightCount[p] += hpWeight;
          aCutOrderCountTotal[p] += hpWeight;
          if (simulationMode && cohSnr2 > 5.0 ) {
            lowBinWeightCount5[p] += hpWeight;
            aCutOrderCountTotal5[p] += hpWeight;
          }
        }
        if (ldInterceptThreshold != 0) {
          //cout << "Bin # " << hpBin << ", ldThres " << ldInterceptThreshold << endl;
          diffHistMap[p][hpBin]->Fill(ldCutVal, hpWeight);
        }
      }
      if (keepDetail) {
        vbprintf(" writing healpix detail event %i bin %i \n", eventNumber, hpBin);
        outputTree1->Fill();
        // update the bin count map
        if (hpEventsAccepted[p].find(hpBin) == hpEventsAccepted[p].end()) {
          hpEventsAccepted[p].insert(pair<int, float>(hpBin, hpWeight));
        } else {
          hpEventsAccepted[p][hpBin] += hpWeight;
        }
        printf(" accepting event %8i %cPol run %3i bin %5i weight %5.4f \n", eventNumber, pChars[p], runNumber, hpBin, hpWeight);
        //printf(" ldCutVal %5.4f ldThresh %5.4f \n", ldCutVal, ldInterceptThreshold);
        //sleep(1);
        fprintf(acceptedEventFile, " %8i %c %3i %5i %5.4f \n", eventNumber, pChars[p], runNumber, hpBin, hpWeight);
        
        // criminally inelegant hack job to snag localization of accepted event because some unnamed idiot didn't do it in the first stage of analysis
        // better hope we don't accept very many events

        // So this part of the code slows shit down like crazy.  heres the deal.
        // Not sure we need it.  also not sure we need it all the time.
        // If we let more than a few events get through the code runs like... something slow.
        // Added flag around it so we can turn it off.  -Jacob!
        // -I added all this info to tree0, so this just needs to be re written to use it.
        // if you need to use this.  It shouldnt be that hard... lets... leave
        // leave it as an excersize for the reader.  -Jacob

        if (snagInterferometry) {
          sprintf(filename, "%s/analyzerResults_run%03i.root", interfFileDir, runNumber);
          printf(" interferometry filename is %s \n", filename);
          TFile* interfResultsFile = new TFile(filename);
          interfResultsFile->cd();
          TTree* interfResultsTree = (TTree*) interfResultsFile->Get("resultTree");
          printf(" interferometry result tree has %lli entries \n", interfResultsTree->GetEntries());
          interfResultsTree->BuildIndex("eventNumber");
          AnitaEventSummary* eventSummary = 0;
          interfResultsTree->SetBranchAddress("eventSummary", &eventSummary);
          int getRes = interfResultsTree->GetEntryWithIndex(eventNumber);
          if (getRes > 0) {
            double lat = eventSummary->peak[p][0].latitude;
            double lon = eventSummary->peak[p][0].longitude;
            double ea, no;
            bedmap->LonLattoEaNo(lon, lat, ea, no);
            printf(" localization of event %i is lat/lon (%f,%f) ea/no (%f,%f) \n", eventNumber, lat, lon, ea, no);
            accEventNumber[p].push_back(eventNumber);
            accEventBin[p].push_back(hpBin);
            accEventWeight[p].push_back(hpWeight);
            accEventLat[p].push_back(lat);
            accEventLon[p].push_back(lon);
            accEventEa[p].push_back(ea);
            accEventNo[p].push_back(no);
          } else {
            printf(" localization of event %i not found \n", eventNumber);
          }
        }
      }
    }
    if ((e+1)%1000==0) {vbprintf("%i entries read \n", e+1);}
  }
  fclose(acceptedEventFile);
  for (int p=0; p<2; ++p) {
    //printf("snrHist0[%i] has %f entries \n", p, snrHist0[p]->GetSumOfWeights());
    // populate the binned continent map
    for (int gBin = 0; gBin <hpCountsHist[p]->GetNcells()-1; ++gBin) {
      int xBin, yBin, zBin;
      hpCountsHist[p]->GetBinXYZ(gBin, xBin, yBin, zBin);
      float ea = hpCountsHist[p]->GetXaxis()->GetBinCenter(xBin) * 1000;
      float no = hpCountsHist[p]->GetYaxis()->GetBinCenter(yBin) * 1000;
      double thisLat, thisLon;
      bedmap->EaNoToLonLat(ea, no, thisLon, thisLat);

      double thisTheta = 0;
      double thisPhi = 0;
      if (!offsetType) {
        thisTheta = (-thisLat+90.0+hpThOffset) * M_PI/180.0;
        thisPhi = (thisLon+hpPhiOffset) * M_PI/180.0;
      } else {
        pair<double,double> thisPoint = rotateOnSphere((90.0-thisLat)*M_PI/180.0, thisLon*M_PI/180.0, hpPhiOffset*M_PI/180.0, hpThOffset*M_PI/180.0);
        thisTheta = thisPoint.first;
        thisPhi = thisPoint.second;
      }
      //while (thisPhi < -M_PI) thisPhi += 2*M_PI;
      pointing point(thisTheta, thisPhi);
      int binNo = healpix->ang2pix(point);
      map<int, float>::iterator thisEntry = hpEventsAcceptedPreFinalCuts[p].find(binNo);
      float eventCounts = 0;
      // change to make the heat map show the number of events in a bin BEFORE the linear disciminate cut! -Jacob!
      if (thisEntry != hpEventsAcceptedPreFinalCuts[p].end()) {
        //
        eventCounts = hpEventsAccepted[p][binNo];
        eventCounts = hpEventsAcceptedPreFinalCuts[p][binNo];
      }
      hpCountsHist[p]->Fill(ea/1000, no/1000, eventCounts);
      //printf(" ea=%8.0f, no=%8.0f, lat=%6.3f, lon=%6.3f, theta=%4.3f, phi=%4.3f, pixelNo=%i, count=%8.1f \n",
      // ea, no, thisLat, thisLon, thisTheta, thisPhi, pixelNo, eventCounts);
      //if (count%100 == 0) {printf(" %8i entries processed", count);}
    }
    
  }
  
  TCanvas* circPeakSepLDCanv = 0;
  TCanvas* circPeakSepSNRCanv = 0;
  TCanvas* corrHilbCanv = 0;
  TCanvas* corrSnrCanv = 0;
  TCanvas* lphiLonCanv = 0;
  TCanvas* minCircValCanv = 0;
  TCanvas* minCircValLDCanv = 0;
  TCanvas* peakRatioCanv = 0;
  TCanvas* snrHist0Canv = 0;
  TCanvas* snrHist1Canv = 0;
  TCanvas* snrHist2Canv = 0;
  TCanvas* snrHist3Canv = 0;

  TCanvas* peakSepCircCanv = 0;
  TCanvas* corrPeakPolCanv = 0;
  TCanvas* cPolPeakDistCanv = 0;
  TCanvas* circPeakCanv = 0;
  TCanvas* hpCountsCanv = 0;

  TCanvas* diffCanv = 0;
  
  circPeakSepSNRCanv = new TCanvas("circPeakSepSNRCanv", "circPeakSepSNRCanv", 400, 800);
  circPeakSepSNRCanv->Divide(1, 2);
  circPeakSepLDCanv = new TCanvas("circPeakSepLDCanv", "circPeakSepLDCanv", 400, 800);
  circPeakSepLDCanv->Divide(1, 2);
  corrHilbCanv = new TCanvas("corrHilbCanv", "Correlation peak vs. Hilbert peak value", 450, 900);
  corrHilbCanv->Divide(1, 2);
  corrSnrCanv = new TCanvas("corrSnrCanv", "Correlation peak vs. SNR", 450, 900);
  corrSnrCanv->Divide(1, 2);
  lphiLonCanv = new TCanvas("lphiLonCanv", "LCP Phi vs Payload Longitude", 450, 900);
  lphiLonCanv->Divide(1, 2);
  minCircValCanv = new TCanvas("minCircValCanv", "minCircValCanv", 450, 900);
  minCircValCanv->Divide(1, 2);
  minCircValLDCanv = new TCanvas("minCircValLDCanv", "minCircValLDCanv", 450, 900);
  minCircValLDCanv->Divide(1, 2);
  peakRatioCanv = new TCanvas("peakRatioCanv", "Coherent sum peak ratio", 450, 900);
  peakRatioCanv->Divide(1, 2);
  snrHist0Canv = new TCanvas("snrHist0Canv", "SNR", 450, 900);
  snrHist0Canv->Divide(1, 2);
  snrHist1Canv = new TCanvas("snrHist1Canv", "SNR", 450, 900);
  snrHist1Canv->Divide(1, 2);
  snrHist2Canv = new TCanvas("snrHist2Canv", "SNR", 450, 900);
  snrHist2Canv->Divide(1, 2);
  snrHist3Canv = new TCanvas("snrHist3Canv", "SNR", 450, 900);
  snrHist3Canv->Divide(1, 2);
  hpCountsCanv = new TCanvas("hpCountsCanv", "hpCounts", 450, 900);
  hpCountsCanv->Divide(1, 2);
  
  peakSepCircCanv = new TCanvas("peakSepCircCanv", "cPol peak separation", 500, 500);
  corrPeakPolCanv = new TCanvas("corrPeakPolCanv", "corr peak by pol", 500, 500);
  cPolPeakDistCanv = new TCanvas("cPolPeakDistCanv", "cPol peak vs distance", 500, 500);
  circPeakCanv = new TCanvas("circPeakCanv", "cPol peak", 500, 500);

  diffCanv = new TCanvas("diffCanv", "LD Hist", 450, 900);
  diffCanv->Divide(1,2);

  for (int p=0; p<2; ++p) {
    if (circPeakSepLDCanv) {
      circPeakSepLDCanv->cd(p+1);
      circPeakSepLDHist[p]->Draw("COLZ");
      gStyle->SetOptStat("eou");
      circPeakSepLDCanv->Draw();
    }
    if (circPeakSepSNRCanv) {
      circPeakSepSNRCanv->cd(p+1);
      circPeakSepSNRHist[p]->Draw("COLZ");
      gStyle->SetOptStat("eou");
      circPeakSepSNRCanv->Draw();
    }
    if (corrHilbCanv) {
      corrHilbCanv->cd(p+1);
      gPad->SetLogz();
      corrHilbHist[p]->Draw("COLZ");
      gStyle->SetOptStat("eou");
      corrHilbCanv->Draw();
    }
    if (corrSnrCanv) {
      corrSnrCanv->cd(p+1);
      gPad->SetLogz();
      corrSnrHist[p]->Draw("COLZ");
      gStyle->SetOptStat("eou");
      corrSnrCanv->Draw();
    }
    if (lphiLonCanv) {
      lphiLonCanv->cd(p+1);
      gPad->SetLogz();
      lphiLonHist[p]->Draw("COLZ");
      gStyle->SetOptStat("eou");
      lphiLonCanv->Draw();
    }
    if (minCircValCanv) {
      minCircValCanv->cd(p+1);
      minCircValHist[p]->Draw("");
      gStyle->SetOptStat("eou");
      minCircValCanv->Draw();
    }
    if (minCircValLDCanv) {
      minCircValLDCanv->cd(p+1);
      minCircValLDHist[p]->Draw("COLZ");
      gStyle->SetOptStat("eou");
      minCircValLDCanv->Draw();
    }
    if (peakRatioCanv) {
      peakRatioCanv->cd(p+1);
      peakRatioHist[p]->Draw();
      gStyle->SetOptStat("eou");
      peakRatioCanv->Draw();
    }
    if (snrHist0Canv) {
      snrHist0Canv->cd(p+1);
      snrHist0[p]->Draw();
      gStyle->SetOptStat("eou");
      snrHist0Canv->Update();
    }
    if (snrHist1Canv) {
      snrHist1Canv->cd(p+1);
      snrHist1[p]->Draw();
      gStyle->SetOptStat("eou");
      snrHist1Canv->Draw();
    }
    if (snrHist2Canv) {
      snrHist2Canv->cd(p+1);
      snrHist2[p]->Draw();
      gStyle->SetOptStat("eou");
      snrHist2Canv->Draw();
    }
    if (snrHist3Canv) {
      snrHist3Canv->cd(p+1);
      snrHist3[p]->Draw();
      gStyle->SetOptStat("eou");
      snrHist3Canv->Draw();
    }
    if (hpCountsCanv) {
      hpCountsCanv->cd(p+1);
      //gStyle->SetOptStat(0);
      coastLineGr->SetMarkerColor(14);
      coastLineGr->SetFillColor(18);
      coastLineGr->SetLineWidth(1);
      coastLineGr->SetFillStyle(1001);
      hpCountsHist[p]->SetStats(kFALSE);
      hpCountsHist[p]->Draw("COLZ");
      gPad->SetLogz();
      //hpCountsCanv->Draw();
      coastLineGr->Draw("C");
      double ea, no; // for bedmap conversion
      for (pair<int, float> thisEntry : hpEventsAccepted[p]) {
        pointing thisPoint = healpix->pix2ang(thisEntry.first);

        double thisLat = 0;
        double thisLon = 0;
        if (offsetType) {
          pair<double,double> thisPointNew = rotateOnSphere(thisPoint.theta, thisPoint.phi, -hpPhiOffset*M_PI/180.0, -hpThOffset*M_PI/180.0);  //new way
          thisLat = 90.0 - (thisPointNew.first * 180.0/M_PI);
          thisLon = thisPointNew.second * 180.0/M_PI;
        } else {
          thisLat = 90.0 - (thisPoint.theta * 180.0/M_PI) + hpThOffset;    //old way of doing rotation
          thisLon = thisPoint.phi * 180.0/M_PI - hpPhiOffset;
        }

        bedmap->LonLattoEaNo(thisLon, thisLat,	ea, no);
        ea /= 1000;
        no /= 1000;
        printf("hp bin %i lat,lon is %f,%f ea,no is %f,%f \n", thisEntry.first, thisLat, thisLon, ea, no);
        char hpBinText[8]; sprintf(hpBinText, "%i", thisEntry.first);
        printf("drawing hp bin number at %f, %f, %f, %f \n", ea-125, no, ea+125, no+100);
        bool posQuad = (ea*no > 0);
        double textEa = posQuad ? ea : ea-75;
        double textNo = posQuad ? no-125 : no+100;
        TText* thisText = new TText(textEa, textNo, hpBinText);
        thisText->SetTextSize(0.025);
        thisText->SetTextAngle(posQuad ? 55 : -55);
        thisText->SetTextFont(42);
        thisText->Draw();
      }
      
      for (int k=0; k<accEventNumber[p].size(); ++k) {
        printf(" marking (%f,%f) \n", accEventEa[p][k]/1000, accEventNo[p][k]/1000);
        TMarker* thisMarker = new TMarker(accEventEa[p][k]/1000, accEventNo[p][k]/1000, kPlus);
        thisMarker->SetMarkerColor(kRed);
        thisMarker->Draw();
        thisMarker = 0;
      }
      bedmap->LonLattoEaNo(MCMURDO[1], MCMURDO[0],	ea, no);
      TMarker* thisMarker = new TMarker(ea/1000, no/1000, kPlus);
      thisMarker->SetMarkerColor(kBlack);
      thisMarker->Draw();
      thisMarker = 0;
      
      //gStyle->SetOptStat(0);
    }
  }
  
  if (diffCanv) {
    for (int thisBin : binsToKeep) {
      map<int, float>::iterator thisHpBinEntry = hpBinInterceptMap.find(thisBin);
      if (thisHpBinEntry != hpBinInterceptMap.end()) {
        for (int p = 0; p<2; p++) {
          diffCanv->cd(p+1);
          diffHistMap[p][thisBin]->Draw();
          gStyle->SetOptStat("eou");
          gPad->SetLogy();
          diffCanv->Draw();
        }
        sprintf(filename, "%s/diffHist0%i_%03i_%03i_%c.png", fileDir, thisBin, startRunNum, endRunNum, whichPol);
        diffCanv->SaveAs(filename);
        sprintf(filename, "%s/diffHist0%i_%03i_%03i_%c.root", fileDir, thisBin, startRunNum, endRunNum, whichPol);
        diffCanv->SaveAs(filename);
        //diffCanv = 0;
      }
    }
  }


  if (peakSepCircCanv) {
    peakSepCircCanv->cd(0);
    gPad->SetLogy();
    peakSepCircHist->Draw("COLZ");
    gStyle->SetOptStat("eou");
    peakSepCircCanv->Draw(0);
  }
  if (corrPeakPolCanv) {
    corrPeakPolCanv->cd(0);
    corrPeakPolHist->Draw("COLZ");
    gStyle->SetOptStat("eou");
    corrPeakPolCanv->Draw();
  }
  if (cPolPeakDistCanv) {
    cPolPeakDistCanv->cd(0);
    cPolPeakDistHist->Draw("COLZ");
    gStyle->SetOptStat("eou");
    cPolPeakDistCanv->Draw();
  }
  if (circPeakCanv) {
    circPeakCanv->cd(0);
    circPeakHist->Draw("COLZ");
    gStyle->SetOptStat("eou");
    circPeakCanv->Draw();
  }

  printf ("result tree 0 has %llu entries \n", outputTree0->GetEntries());
  printf ("result tree 1 has %llu entries \n", outputTree1->GetEntries());
  //printf ("current directory is %s \n", gDirectory->GetName());
  //printf ("parent directory is %s \n", gDirectory->GetMotherDir()->GetName());
  
  // print out the cut table
  char cutTableStr[16384] = "";
  sprintf(cutTableStr+strlen(cutTableStr), "\nAnalysis cuts: \n");
  for (int p=0; p<2; ++p) {
    char pol = (p==0) ? 'H' : 'V';
    sprintf(cutTableStr+strlen(cutTableStr), " polarization %c\n", pol);
    sprintf(cutTableStr+strlen(cutTableStr), " total events processed %8.1f \n", eventCount[p]);
    sprintf(cutTableStr+strlen(cutTableStr), " description as first cut as ordered cut as last cut \n");
    sprintf(cutTableStr+strlen(cutTableStr), " number fraction number fraction number fraction \n");
    if (false) {
      for (int k=0; k<aCutOrderStage1.size(); ++k) {
        sprintf(cutTableStr+strlen(cutTableStr), "%28s %8.1f %8.5f %8.1f %8.5f %8.1f %8.5f \n",
                A_CUT_DESC[aCutOrderStage1[k]] ,aCutCounts[p][aCutOrderStage1[k]], (float)aCutCounts[p][aCutOrderStage1[k]]/eventCount[p],
                aCutOrderCounts[p][aCutOrderStage1[k]], (float)aCutOrderCounts[p][aCutOrderStage1[k]]/eventCount[p],
                aCutLastCounts[p][aCutOrderStage1[k]], (float)aCutLastCounts[p][aCutOrderStage1[k]]/eventCount[p]);
      }
      sprintf(cutTableStr+strlen(cutTableStr), "\n");
    }
    for (int k=0; k<aCutOrderStage2.size(); ++k) {
      sprintf(cutTableStr+strlen(cutTableStr), "%28s %8.1f %8.5f %8.1f %8.5f %8.1f %8.5f \n",
              A_CUT_DESC[aCutOrderStage2[k]] ,aCutCounts[p][aCutOrderStage2[k]], (float)aCutCounts[p][aCutOrderStage2[k]]/eventCount[p],
              aCutOrderCounts[p][aCutOrderStage2[k]], (float)aCutOrderCounts[p][aCutOrderStage2[k]]/eventCount[p],
              aCutLastCounts[p][aCutOrderStage2[k]], (float)aCutLastCounts[p][aCutOrderStage2[k]]/eventCount[p]);
    }
    sprintf(cutTableStr+strlen(cutTableStr), "%28s %8.1f %8.5f \n", "healpix bin rejected",
            hpBinRejectedCount[p], hpBinRejectedCount[p]/eventCount[p]);
    sprintf(cutTableStr+strlen(cutTableStr), "%28s %8.1f %8.5f \n", "linear discriminant cut",
            ldRejectedCount[p], ldRejectedCount[p]/eventCount[p]);
    sprintf(cutTableStr+strlen(cutTableStr), "%28s %8.1f %8.5f \n", "low event bin weight",
            lowBinWeightCount[p], lowBinWeightCount[p]/eventCount[p]);
    sprintf(cutTableStr+strlen(cutTableStr), "\n");
    sprintf(cutTableStr+strlen(cutTableStr), " total events cut: %8.1f %8.5f \n",
            aCutOrderCountTotal[p], (float)aCutOrderCountTotal[p]/eventCount[p]);
    sprintf(cutTableStr+strlen(cutTableStr), " surviving events: %8.1f %8.5f \n\n",
            eventCount[p]-aCutOrderCountTotal[p], (float)(eventCount[p]-aCutOrderCountTotal[p])/eventCount[p]);
  }
  printf(cutTableStr);

  // print out cuttable for sim events with snr > 5, if relevent.
  if (simulationMode) {
    // print out the cut table
    char cutTableStr5[16384] = "";
    sprintf(cutTableStr5+strlen(cutTableStr5), "\nAnalysis cuts: \n");
    for (int p=0; p<2; ++p) {
      char pol = (p==0) ? 'H' : 'V';
      sprintf(cutTableStr5+strlen(cutTableStr5), " polarization %c\n", pol);
      sprintf(cutTableStr5+strlen(cutTableStr5), " total events processed %8.1f \n", eventCount5[p]);
      sprintf(cutTableStr5+strlen(cutTableStr5), " description as first cut as ordered cut as last cut \n");
      sprintf(cutTableStr5+strlen(cutTableStr5), " number fraction number fraction number fraction \n");
      if (false) {
        for (int k=0; k<aCutOrderStage1.size(); ++k) {
          sprintf(cutTableStr5+strlen(cutTableStr5), "%28s %8.1f %8.5f %8.1f %8.5f %8.1f %8.5f \n",
                  A_CUT_DESC[aCutOrderStage1[k]] ,aCutCounts5[p][aCutOrderStage1[k]], (float)aCutCounts5[p][aCutOrderStage1[k]]/eventCount5[p],
                  aCutOrderCounts5[p][aCutOrderStage1[k]], (float)aCutOrderCounts5[p][aCutOrderStage1[k]]/eventCount5[p],
                  aCutLastCounts5[p][aCutOrderStage1[k]], (float)aCutLastCounts5[p][aCutOrderStage1[k]]/eventCount5[p]);
        }
        sprintf(cutTableStr5+strlen(cutTableStr5), "\n");
      }
      for (int k=0; k<aCutOrderStage2.size(); ++k) {
        sprintf(cutTableStr5+strlen(cutTableStr5), "%28s %8.1f %8.5f %8.1f %8.5f %8.1f %8.5f \n",
                A_CUT_DESC[aCutOrderStage2[k]] ,aCutCounts5[p][aCutOrderStage2[k]], (float)aCutCounts5[p][aCutOrderStage2[k]]/eventCount5[p],
                aCutOrderCounts5[p][aCutOrderStage2[k]], (float)aCutOrderCounts5[p][aCutOrderStage2[k]]/eventCount5[p],
                aCutLastCounts5[p][aCutOrderStage2[k]], (float)aCutLastCounts5[p][aCutOrderStage2[k]]/eventCount5[p]);
      }
      sprintf(cutTableStr5+strlen(cutTableStr5), "%28s %8.1f %8.5f \n", "healpix bin rejected",
              hpBinRejectedCount5[p], hpBinRejectedCount5[p]/eventCount5[p]);
      sprintf(cutTableStr5+strlen(cutTableStr5), "%28s %8.1f %8.5f \n", "linear discriminant cut",
              ldRejectedCount5[p], ldRejectedCount5[p]/eventCount5[p]);
      sprintf(cutTableStr5+strlen(cutTableStr5), "%28s %8.1f %8.5f \n", "low event bin weight",
              lowBinWeightCount5[p], lowBinWeightCount5[p]/eventCount5[p]);
      sprintf(cutTableStr5+strlen(cutTableStr5), "\n");
      sprintf(cutTableStr5+strlen(cutTableStr5), " total events cut: %8.1f %8.5f \n",
              aCutOrderCountTotal5[p], (float)aCutOrderCountTotal5[p]/eventCount5[p]);
      sprintf(cutTableStr5+strlen(cutTableStr5), " surviving events: %8.1f %8.5f \n\n",
              eventCount5[p]-aCutOrderCountTotal5[p], (float)(eventCount5[p]-aCutOrderCountTotal5[p])/eventCount5[p]);
    }
    printf(cutTableStr5);
  }
  
  if (saveOutput) {
    outputFile->cd();
    outputTree0->SetDirectory(outputFile);
    outputTree1->SetDirectory(outputFile);
    outputTree0->Write();
    outputTree1->Write();
    outputFile->Save();

    if (circPeakSepLDCanv) {
      sprintf(filename, "%s/circPeakSepLDCanv_%03i_%03i_%s.png", fileDir, startRunNum, endRunNum, pChars);
      circPeakSepLDCanv->SaveAs(filename);
      sprintf(filename, "%s/circPeakSepLDCanv_%03i_%03i_%s.root", fileDir, startRunNum, endRunNum, pChars);
      circPeakSepLDCanv->SaveAs(filename);
    }
    if (circPeakSepSNRCanv) {
      sprintf(filename, "%s/circPeakSepSNRCanv_%03i_%03i_%s.png", fileDir, startRunNum, endRunNum, pChars);
      circPeakSepSNRCanv->SaveAs(filename);
      sprintf(filename, "%s/circPeakSepSNRCanv_%03i_%03i_%s.root", fileDir, startRunNum, endRunNum, pChars);
      circPeakSepSNRCanv->SaveAs(filename);
    }
    if (corrHilbCanv) {
      sprintf(filename, "%s/corrHilbCanv_%03i_%03i_%s.png", fileDir, startRunNum, endRunNum, pChars);
      corrHilbCanv->SaveAs(filename);
      sprintf(filename, "%s/corrHilbCanv_%03i_%03i_%s.root", fileDir, startRunNum, endRunNum, pChars);
      corrHilbCanv->SaveAs(filename);
    }
    if (corrSnrCanv) {
      sprintf(filename, "%s/corrSnrCanv_%03i_%03i_%s.png", fileDir, startRunNum, endRunNum, pChars);
      corrSnrCanv->SaveAs(filename);
      sprintf(filename, "%s/corrSnrCanv_%03i_%03i_%s.root", fileDir, startRunNum, endRunNum, pChars);
      corrSnrCanv->SaveAs(filename);
    }
    if (lphiLonCanv) {
      sprintf(filename, "%s/lphiLonCanv_%03i_%03i_%s.png", fileDir, startRunNum, endRunNum, pChars);
      lphiLonCanv->SaveAs(filename);
      sprintf(filename, "%s/lphiLonCanv_%03i_%03i_%s.root", fileDir, startRunNum, endRunNum, pChars);
      lphiLonCanv->SaveAs(filename);
    }
    if (minCircValCanv) {
      sprintf(filename, "%s/minCircValCanv_%03i_%03i_%s.png", fileDir, startRunNum, endRunNum, pChars);
      minCircValCanv->SaveAs(filename);
      sprintf(filename, "%s/minCircValCanv_%03i_%03i_%s.root", fileDir, startRunNum, endRunNum, pChars);
      minCircValCanv->SaveAs(filename);
    }
    if (minCircValLDCanv) {
      sprintf(filename, "%s/minCircValLDCanv_%03i_%03i_%s.png", fileDir, startRunNum, endRunNum, pChars);
      minCircValLDCanv->SaveAs(filename);
      sprintf(filename, "%s/minCircValLDCanv_%03i_%03i_%s.root", fileDir, startRunNum, endRunNum, pChars);
      minCircValLDCanv->SaveAs(filename);
    }
    if (peakRatioCanv) {
      sprintf(filename, "%s/peakRatioCanv_%03i_%03i_%s.png", fileDir, startRunNum, endRunNum, pChars);
      peakRatioCanv->SaveAs(filename);
      sprintf(filename, "%s/peakRatioCanv_%03i_%03i_%s.root", fileDir, startRunNum, endRunNum, pChars);
      peakRatioCanv->SaveAs(filename);
    }
    if (peakSepCircCanv) {
      sprintf(filename, "%s/peakSepCircCanv_%03i_%03i_%s.png", fileDir, startRunNum, endRunNum, pChars);
      peakSepCircCanv->SaveAs(filename);
      sprintf(filename, "%s/peakSepCircCanv_%03i_%03i_%s.root", fileDir, startRunNum, endRunNum, pChars);
      peakSepCircCanv->SaveAs(filename);
    }
    if (corrPeakPolCanv) {
      sprintf(filename, "%s/corrPeakPolCanv_%03i_%03i_%s.png", fileDir, startRunNum, endRunNum, pChars);
      corrPeakPolCanv->SaveAs(filename);
      sprintf(filename, "%s/corrPeakPolCanv_%03i_%03i_%s.root", fileDir, startRunNum, endRunNum, pChars);
      corrPeakPolCanv->SaveAs(filename);
    }
    if (cPolPeakDistCanv) {
      sprintf(filename, "%s/cPolPeakDistCanv_%03i_%03i_%s.png", fileDir, startRunNum, endRunNum, pChars);
      cPolPeakDistCanv->SaveAs(filename);
      sprintf(filename, "%s/cPolPeakDistCanv_%03i_%03i_%s.root", fileDir, startRunNum, endRunNum, pChars);
      cPolPeakDistCanv->SaveAs(filename);
    }
    if (circPeakCanv) {
      sprintf(filename, "%s/circPeakCanv_%03i_%03i_%s.png", fileDir, startRunNum, endRunNum, pChars);
      circPeakCanv->SaveAs(filename);
      sprintf(filename, "%s/circPeakCanv_%03i_%03i_%s.root", fileDir, startRunNum, endRunNum, pChars);
      circPeakCanv->SaveAs(filename);
    }

    if (snrHist0Canv) {
      sprintf(filename, "%s/snr0_%03i_%03i_%s.png", fileDir, startRunNum, endRunNum, pChars);
      snrHist0Canv->SaveAs(filename);
      sprintf(filename, "%s/snr0_%03i_%03i_%s.root", fileDir, startRunNum, endRunNum, pChars);
      snrHist0Canv->SaveAs(filename);
    }
    if (snrHist1Canv) {
      sprintf(filename, "%s/snr1_%03i_%03i_%s.png", fileDir, startRunNum, endRunNum, pChars);
      snrHist1Canv->SaveAs(filename);
      sprintf(filename, "%s/snr1_%03i_%03i_%s.root", fileDir, startRunNum, endRunNum, pChars);
      snrHist1Canv->SaveAs(filename);
    }
    if (snrHist2Canv) {
      sprintf(filename, "%s/snr2_%03i_%03i_%s.png", fileDir, startRunNum, endRunNum, pChars);
      snrHist2Canv->SaveAs(filename);
      sprintf(filename, "%s/snr2_%03i_%03i_%s.root", fileDir, startRunNum, endRunNum, pChars);
      snrHist2Canv->SaveAs(filename);
    }
    if (snrHist3Canv) {
      sprintf(filename, "%s/snr3_%03i_%03i_%s.png", fileDir, startRunNum, endRunNum, pChars);
      snrHist3Canv->SaveAs(filename);
      sprintf(filename, "%s/snr3_%03i_%03i_%s.root", fileDir, startRunNum, endRunNum, pChars);
      snrHist3Canv->SaveAs(filename);
    }
    if (hpCountsCanv) {
      sprintf(filename, "%s/hpCounts_%03i_%03i_%s.png", fileDir, startRunNum, endRunNum, pChars);
      hpCountsCanv->SaveAs(filename);
      sprintf(filename, "%s/hpCounts_%03i_%03i_%s.root", fileDir, startRunNum, endRunNum, pChars);
      hpCountsCanv->SaveAs(filename);
    }

   outputFile->Close();
   
  }  
  printf("goodbye world \n");
  //app->Run();
  return 0;
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
    if (lengthsquared>R_EARTH*R_EARTH) { return 1; }
    else { return 0; }
  }
  cout << "I shouldn't be here.\n";
  return 0;
}


// return center (Ea, No), radial and lateral semi-axes
vector<float> surfaceEllipseParms1(float phi, float theta, float plSrcDist, float errPh, float errTh, UsefulAdu5Pat* gps, BedmapReader* bedmap) {
  // use bedmap to obtain ellipse left, right, near and far vertices, easting/northing
  vbprintf("phi=%f, theta=%f, dist=%f, errPhi=%f, errTheta=%f \n", phi, theta, plSrcDist, errPh, errTh);
  double thisLat, thisLon, thisAlt, thisAdj;
  double far[2], near[2];
  vector<float> fail = {-9999, -9999, -9999, -9999};
  // pointing peak (just for test)
  gps->traceBackToContinent(phi *M_PI/180.0, theta *M_PI/180.0, &thisLon, &thisLat, &thisAlt, &thisAdj);
  if (thisLat == -9999 || thisLon == -9999 || thisAlt == -9999) {
    vbprintf("trace back failed for peak location \n");
    return fail; // structure violation yeah yeah
  }
  // far vertex
  gps->traceBackToContinent(phi *M_PI/180.0, (theta-errTh) *M_PI/180.0, &thisLon, &thisLat, &thisAlt, &thisAdj);
  if (thisLat == -9999 || thisLon == -9999 || thisAlt == -9999) {
    vbprintf("trace back failed (%.2f,%.2f) for far ellipse vertex theta-far=%f, theta-near=%f \n", thisLat, thisLon, theta-errTh, theta+errTh);
    return fail; // structure violation yeah yeah
  }
  vbprintf(" far vertex lat=%.2f, lon=%.2f, alt=%.2f, adj=%.2f \n", thisLat, thisLon, thisAlt, thisAdj);
  bedmap->LonLattoEaNo(thisLon, thisLat, far[0], far[1]);
  // near vertex
  gps->traceBackToContinent(phi *M_PI/180.0, (theta+errTh) *M_PI/180.0, &thisLon, &thisLat, &thisAlt, &thisAdj);
  if (thisLat == -9999 || thisLon == -9999 || thisAlt == -9999) {
    vbprintf("trace back failed for near ellipse vertex theta-far=%f, theta-near=%f \n", theta-errTh, theta+errTh);
    return fail; // structure violation
  }
  vbprintf(" near vertex lat=%.2f, lon=%.2f, alt=%.2f, adj=%.2f \n", thisLat, thisLon, thisAlt, thisAdj);
  bedmap->LonLattoEaNo(thisLon, thisLat, near[0], near[1]);

  float rAxis = 0.5*sqrt((far[0]-near[0])*(far[0]-near[0]) + (far[1]-near[1])*(far[1]-near[1]));
  float centerX = (far[0]+near[0])*0.5;
  float centerY = (far[1]+near[1])*0.5;
  float thAxis = plSrcDist * errPh * M_PI/180.0;; // errPh assumed small
  return vector<float>({rAxis, thAxis, centerX, centerY});
}

// Function to rotate a point by an angle A around the x-axis and an angle B
// around the y-axis, staying on the surface of the sphere.
//  Takes the initial data point, (theta,phi), and the rotation angles A and B.
pair<double,double> rotateOnSphere(double theta, double phi, double A, double B) 
{
  // rotate by B first.  (note: atan(y/x) = atan2(y,x))
  double thetaPrime = atan2( 
          pow( pow(sin(theta)*cos(phi),2) + pow(cos(A)*sin(theta)*sin(phi)-sin(A)*cos(theta),2),0.5 ),
          cos(A)*cos(theta)+sin(A)*sin(theta)*sin(phi)
         );
  double phiPrime = atan2(
          cos(A)*sin(theta)*sin(phi) - sin(A)*cos(theta),
          sin(theta)*cos(phi)
         );

  // now rotate by A.
  double thetaPrimePrime = atan2(
          pow( pow(sin(thetaPrime)*sin(phiPrime),2) + pow(cos(B)*sin(thetaPrime)*cos(phiPrime)+sin(B)*cos(thetaPrime),2),0.5 ),
          cos(B)*cos(thetaPrime)-sin(B)*sin(thetaPrime)*cos(phiPrime)
         );
  double phiPrimePrime = atan2(
          sin(thetaPrime)*sin(phiPrime),
          cos(B)*sin(thetaPrime)*cos(phiPrime)+sin(B)*cos(thetaPrime)
         );

  return pair<double,double>(thetaPrimePrime,phiPrimePrime);
}



// diffCutPlots.cxx 
//  read results from Analysis cuts pass 1
//  calculate the diagonal cut for each event
//  accumulate it into histograms by cut value, one for each pixel
// Edits by Jacob Gordon:
//  - Changed fitting to be done by minuit
//  - Changed pval pseudo experiment data to be found using root TF1 ditribution utility
//  - Added treatment of event weight to pseudo experiment
//  - Changed ll_val calc to use minuit (just like the fitting) 
//  - Added print out of additional files for the case of pval = 1
//  - Added error calc to backgrounds from parameter uncertainties
//  - Added calculation of leakage between bins
//  - Added uncertanty due to leakage between bins to S_up calculation
//  - Added options to recalculate healpix bins here


#include "TFile.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "Math/SpecFuncMathCore.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFitResult.h"
#include "TVirtualFitter.h"
#include "InterfUtil.h"
#include "TRandom.h"
#include "TMinuit.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TText.h"
#include "TFrame.h"
#include <sstream>


#include "RawAnitaHeader.h"
#include "AnitaEventSummary.h"
#include "Adu5Pat.h"
#include "UsefulAdu5Pat.h"
//#include "WaveformInfo.h"
#include "AnitaConventions.h"
#include "BedmapReader.h"
#include "AnitaGeomTool.h"
#include "UCUtil.h"



#include <healpix_base.h>

#include "fitFunc.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnContours.h"

#include "gsl/gsl_roots.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_math.h"


int processKeywordParm(char* kwp);

struct dataRow_t {
  float numEvents = 0;
  TH1F* cutHist = 0;
  TH1F* cutHistPse = 0;
  TH1F* cutHistTotal = 0;
  TH1F* cutHistWeighted = 0;
  TH1F* simCutHist = 0;
  TH1F* sHist = 0;
  TH1F* sUpHist = 0;
  TH1F* s_sUpHist = 0;
  TH2F* corrPeakSnrHist = 0;
  TH2F* corrPeakSnrSimHist = 0;
  TH1F* eventSNRSimHist = 0;
  TH1F* passingEventSNRSimHist = 0;
  TH1F* effSNRSimHist = 0;

  int   status = -99;
  int   binStatus = 0;

  //params for minuit2 fit
  float fitSlope = 0;
  float fitSlopeError = 0;
  float fitIntercept = 0;
  float fitInterceptError = 0;
  //params for minuit2 powerlaw fit
  float fitPowScale = 0;
  float fitPowOrder = 0;
  //params for cov matrix
  float fitSlopeM = 0;
  float fitInterceptM = 0;
  float fitCrossM = 0;
  //params for error contour info
  float fitErrorSlopeAvg = 0;
  float fitErrorInterceptAvg = 0;
  float fitErrorSlopeSigma = 0;
  float fitErrorInterceptSigma = 0;
  float fitErrorRow = 0;
  pair <float,float> fitContourMax; 
  pair <float,float> fitContourMin;

  int   startBin = 0;
  int   endBin = 0;
  float xMin = 0;
  float xMax = 0;
  float xMaxPse = 0;
  float xMaxPowPse = 0;
  float yMin = 0;
  float yMax = 0;

  float logL = 0;
  TH1F* logLikeHist = 0;
  float logLikePVal = 0;
  float logLPow = 0;
  TH1F* logLikeHistPow = 0;
  float logLikePValPow = 0;

  float optCutVal = 0;
  // background stuff
  float expBackground = 0;
  float expBGLowError = 0;           //low background bound
  float expBGHighError = 0;          //high background bound
  float fitErrorHigh = 0;            //high error bar from fit uncertainty
  float fitErrorLow = 0;             //low error bar from fit uncertainty
  float outflow = 0;                 //low error bar from spillover
  float inflow = 0;                  //high error bar from spillover
  float samplePassing = 0;
  float sampleInFitRange = 0;
  float sampleInFitRangeTotal = 0;
  float sampleInFitRangeWeighted = 0;
  float simPassing = 0;
  float simEventsBeforeLDCut = 0;
};

struct dataRow_final {
  float simPassing = 0;
  float optCutVal = 0;
  float expBackground = 0;
  float expBGHighError = 0;
  float expBGLowError = 0;
  TH1F* bgHist = 0;
};
  
struct f_parms { double b; double alpha; };   // parm structure for likelihood function
double likelihood0W(double arg, void* parms);
bool sortFunc(pair<int, dataRow_t> r1, pair<int, dataRow_t> r2);

vector<float> surfaceEllipseParms1(float phi, float theta, float dist, float errPh, float errTh,
        UsefulAdu5Pat* gps, BedmapReader* bedmap);

float getRandLogNormal(float mean, float var);

pair<float,float> getRandParamBackground ( float sAvg, float iAvg, float sSigma, float iSigma, float row );

pair<float,float> getSpilloverOutIn( map<int, map<int,pair<double,int> > > spillover,  map<int, dataRow_t> cutValHist, int thisHpBin, float cutVal);


float cutSlope = -6.0;
float circPeakSepThreshold = 46.0;
float circPeakStrengthThreshold = 0.015;
float sampleFrac = 0.1;
float minPVal = 0.05;
float hpThOffset = -0.0;
float hpPhiOffset = 0.0;


int main(int argc, char* argv[]) {
  printf("Hello World:  %s \n", argv[0]);
  cout << "This is Oindree.." << endl; 
  char filename[1024];  
  int startRun = 175;
  int endRun = 439;
  float startScaleFac = 0.1;                 //note this get changed down below if pol=H
  float endScaleFac = 0.001;
  float scaleFacInc = pow(10.0, 0.025);
  float scaleFacIncFine = pow(10.0,0.0025);
  
  //int num_bg_exp = 100;
  int num_bg_exp = 1; // oindree trying to make it run fast during test
  int trailingZeros = 21;                     // number of trailing zero data points to include in histograms

  int pNum = 0;
  bool doPValFlatness = false;
  bool doRebinning = false;
  bool loadSpilloverFromFile = true;
  bool drawIt = false;
  int fitErrorFlag = 1;                        // fitErrorFlag defines how the fit error is calculated
  int ourBin = 3032;                           //  0 = based on uncertainty in minuit exp fit parameters
  int pseudoExpCount = 10000;                  //  1 = based on difference between a fit using exponential vs power law
  char whichPol = 'V';
  bool noisyMode = false;
  char inputDir[1024] = "/users/PAS0174/osu8620/anita/SamsAnalysis/results";
  //char inputDir[1024] = "../results";
  //char inputDir[1024] = "../analysisResults";
  for (int i=1; i<argc; ++i) {
    if (argv[i][0] == '-') {
      if (strlen(argv[i]) > 1) {
        if (argv[i][1] == 'P' || argv[i][1] == 'p') {
          if (strlen(argv[i]) > 2) {
            whichPol = toupper(argv[i][2]);
          }
        }
        else if (argv[i][1] == 'r' || argv[i][1] == 'R') {doRebinning = true;} 
        else if (argv[i][1] == 'v') {verboseMode = true;}
        else if (argv[i][1] == 'V') {verboseMode = true; noisyMode = true;}
        else if (argv[i][1] == 'g' || argv[i][1] == 'G') {drawIt = true;}
        else if (argv[i][1] == 'f' || argv[i][1] == 'F') {doPValFlatness = true;}
        else if (argv[i][1] == '-') {
          processKeywordParm(argv[i]); }
        else if (argv[i][1] == 'D' || argv[i][1] == 'd') {
          if (strlen(argv[i]) > 2) {
            //char thisParm[8];
            unsigned int k=2;
            for (; k<=strlen(argv[i]); k++) {
              inputDir[k-2] = argv[i][k];
            }
            inputDir[k-2] = 0;
          }
        }                    
        //else if (argv[i][1] == 'O' || argv[i][1] == 'o') {
        //  if (strlen(argv[i]) > 2) {
        //    //char thisParm[8];
        //    unsigned int k=2;
        //    for (; k<=strlen(argv[i]); k++) {
        //      outputDir[k-2] = argv[i][k];
        //    }
        //    outputDir[k-2] = 0;
        //  }
        //}                    
      }
    } else {
      if (++pNum==1) {startRun = atoi(argv[i]);}
      else if (pNum==2) {endRun = atoi(argv[i]);}
      else if (pNum==3) {ourBin = atoi(argv[i]);}
    }
  }

  if(whichPol == 'H') {
    startScaleFac = 0.1;       
    endScaleFac = 0.01;
  }

  char outputDir[1024];
  //sprintf(outputDir, "%s/optResults_no_cpol_cut/cpp%02.1f/cps%06.4f/slope%02.0f", inputDir, circPeakSepThreshold, circPeakStrengthThreshold, -cutSlope);
  //sprintf(outputDir, "%s/optResults/cpp%02.1f/cps%06.4f/slope%02.0f", inputDir, circPeakSepThreshold, circPeakStrengthThreshold, -cutSlope);
  
  //oindree-- going to comment out below line and add a different line
  sprintf(outputDir, "%s/../resultsHighBg/optResults_%c/phiOffset%02.1f/thOffset%06.4f", inputDir, whichPol, hpPhiOffset, hpThOffset);
  //sprintf(outputDir, "/users/PAS0654/osu0426/Diffuse/results/optResults_%c/phiOffset%02.1f/thOffset%06.4f", whichPol, hpPhiOffset, hpThOffset);

  //char sumOutputDir[1024];
  //sprintf(sumOutputDir, "%s/%2.1f/%1.3f", inputDir, circPeakSepThreshold, circPeakStrengthThreshold);
  printf("cut slope = %f \n", cutSlope);    
  printf("cpol peak separation threshold = %f \n", circPeakSepThreshold);
  printf("cpol peak strength threshold = %f \n", circPeakStrengthThreshold);
  printf("theta healpix offset = %f \n", hpThOffset);
  printf("phi healix offset = %f \n", hpPhiOffset);
  char outputDirCmd[1024];
  sprintf(outputDirCmd, "mkdir -p %s", outputDir);
  printf("output directory command is %s \n", outputDirCmd);
  system(outputDirCmd);
  char logFilename[1024]; sprintf(logFilename, "%s/diffCutPlots_%i_%i.log", outputDir, startRun, endRun);
  printf("log file is %s \n", logFilename);
  logFile = fopen(logFilename, "w");
  //fprintf(logFile, "test");
  if (!drawIt) {
    gROOT->SetBatch(kTRUE);
  }
  char fitOpt[16] = "WL S";  //Jacob! Prob need to change this up/it might not be needed now?
  //if (!verboseMode) strcat(fitOpt, " Q");
  if (!noisyMode) strcat(fitOpt, " Q");
  if (noisyMode) strcat(fitOpt, " V");
  printf("rotated cut slope value = %f \n", cutSlope);
  
  gRandom->SetSeed(0);
  TApplication* app = 0;
  app = new TApplication("app", &argc, argv);  

  //char inFilename[1024]; sprintf(inFilename, "%s/analysisOutput_2_%03i_%03i.root", inputDir, startRun, endRun);  //oindree changing to below
  char inFilename[1024]; sprintf(inFilename, "/users/PAS0174/osu8620/anita/SamsAnalysis/resultsNewSim/analysisOutput_2_%03i_%03i.root", startRun, endRun);
  //char inFilename[1024]; sprintf(inFilename, "%s/analysisOutput_2_%03i_%03i_CpolInfo.root", inputDir, startRun, endRun);
  lprintf("input filename is %s \n", inFilename);

  char dummyFilename[1024]; sprintf(dummyFilename, "%s/dummy/dummy_%c_%f_%f.root", inputDir, whichPol, hpThOffset, hpPhiOffset);

  TFile* inFile = new TFile(inFilename);
  TTree* tree0 = (TTree*)inFile->Get("analysisOutput0");
  lprintf("input tree 0 has %lli entries \n", tree0->GetEntries());
  TTree* tree1 = 0;
  if (!doRebinning) {
    tree1 = (TTree*)inFile->Get("analysisOutput1");
    lprintf("input tree 1 has %lli entries \n", tree1->GetEntries());
  }

  TFile* dummyFile = new TFile(dummyFilename, "RECREATE");
  if (doRebinning) {
    tree1 = new TTree("analysisOutput1", "analysis1");
  }

  //int oRunNumber = 0;
  //int oEventNumber = 0;
  //char oPol = 0;
  //int oHpBin = 0;
  //float oHpWeight = 0;


  int eventNumber = 0;
  int runNumber = 0;
  char pol = '0';
  float hPeak = 0;
  float cPeak = 0;
  float cohSnr = 0;
  float cohSnr1 = 0;
  float cohSnr2 = 0;
  float cohSnr3 = 0;
  float linPolFrac = 0;
  float minCircPeakVal = 0;
  float cPolDist = 0;
  float cPolPeak[2] = {0};
  int hpBin = 0;
  float hpWeight = 0;

  float theta = 0;
  float phi = 0;
  float ea = 0;
  float no = 0;
  float lat = 0;
  float lon = 0;
  float alt = 0;
  //float ellParms[4] = {0,0,0,0};
  //UsefulAdu5Pat* gps = 0;
  Adu5Pat* gpsRaw = 0;
    



  tree0->SetBranchAddress("runNumber", &runNumber);
  tree0->SetBranchAddress("eventNumber", &eventNumber);
  tree0->SetBranchAddress("pol", &pol);
  tree0->SetBranchAddress("hPeak", &hPeak);
  tree0->SetBranchAddress("cPeak", &cPeak);
  tree0->SetBranchAddress("cohSnr", &cohSnr);
  tree0->SetBranchAddress("cohSnr1", &cohSnr1);
  tree0->SetBranchAddress("cohSnr2", &cohSnr2);
  tree0->SetBranchAddress("cohSnr3", &cohSnr3);
  tree0->SetBranchAddress("linPolFrac", &linPolFrac);
  tree0->SetBranchAddress("minCircPeakVal", &minCircPeakVal);
  tree0->SetBranchAddress("cPolDist", &cPolDist);
  tree0->SetBranchAddress("cPolPeak", &cPolPeak);    

  if ( !doRebinning ) {
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
  
  //tree0->SetBranchAddress("ellParms", &ellParms);                                   // error ellipse parameters
  tree0->SetBranchAddress("gpsRaw", &gpsRaw);                                         // raw payload gps info
  //tree0->SetBranchAddress("gps", &gps);                                             // payload gps info
  tree0->SetBranchAddress("theta", &theta);                                           // event theta
  tree0->SetBranchAddress("phi", &phi);                                               // event phi
  tree0->SetBranchAddress("ea", &ea);                                                 // event easting
  tree0->SetBranchAddress("no", &no);                                                 // event northing
  tree0->SetBranchAddress("lat", &lat);                                               // event latitude
  tree0->SetBranchAddress("lon", &lon);                                               // event longitude
  tree0->SetBranchAddress("alt", &alt);                                               // event altitude



  tree0->BuildIndex("eventNumber", "pol");
  //tree1->BuildIndex("eventNumber", "pol");

  //TTreeIndex* tree1Index = (TTreeIndex*)tree1->GetTreeIndex();
  TTreeIndex* tree0Index = (TTreeIndex*)tree0->GetTreeIndex();

  char histTitle[128];
  float cutValMin = 0;
  //float cutValMax = 25;
  float cutValMax = 100;
  float cutValBins = 2*cutValMax;
  float cutValSimBins = 10*cutValMax;

  map<int, dataRow_t> cutValHist;
  int prevEvent = 0;
  char prevPol = 0;
  float cutVal = 0;
  //  float slope = -corrPeakDiagThresh/hilbPeakDiagThresh;
  lprintf(" discriminant slope = %.5f \n", cutSlope);
  
  sprintf(histTitle, "SNR vs Corr Peak run %i-%i;corr peak;SNR", startRun, endRun);
  TH2F* corrSnrHist = new TH2F("corrSnrHist", histTitle, 50, 0, 0.5, 50, 0, 15);
  bool useThisEntry = false;

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
  Healpix_Ordering_Scheme scheme = Healpix_Ordering_Scheme::RING;
  T_Healpix_Base<int>* healpix = new T_Healpix_Base<int>(healPixN, scheme);
  //lprintf("healpix map has %i pixels \n", healpix->Npix());

  float thError = 0.25;     //this might depend on SNR later?
  float phiError = 0.50;

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
      useThisEntry = true;
      if (getRes <= 0) {lprintf("database error event %i pol %c - skipping \n", eventNumber, pol); useThisEntry = false;}
      if (circPeakSepThreshold>0 && cPolDist>circPeakSepThreshold) {useThisEntry = false;}
      if (circPeakStrengthThreshold>0 && minCircPeakVal<circPeakStrengthThreshold) {useThisEntry = false;}
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
      //vbprintf("     Error hexagon vertices :"); for (int k=0; k<6; ++k) {vbprintf("       %f, %f \n", errHexVertices[0][k], errHexVertices[1][k$
      errHexPixels = ellipsePixels(errEllipseParms1[2], errEllipseParms1[3], errEllipseParms1[0], errEllipseParms1[1], plSrcPhi);


      map<int, float> hpWeights = map<int, float>();
      for (int k=0; k<errHexPixels[0].size(); ++k) {
        double thisLat, thisLon;
        bedmap->EaNoToLonLat(errHexPixels[0][k], errHexPixels[1][k], thisLon, thisLat);
        //vbprintf(" pixel ea/no = %f,%f \n", errHexPixels[0][k], errHexPixels[1][k]);
        //vbprintf(" pixel lat/lon = %f,%f \n", thisLat, thisLon);
        double thisTheta = (-thisLat + 90.0 + hpThOffset) * M_PI/180.0;
        double thisPhi = (thisLon + hpPhiOffset) * M_PI/180.0;
        while (thisPhi < -M_PI) thisPhi += M_PI;
        //vbprintf("  pixel theta/phi = %f,%f \n", thisTheta, thisPhi);
        pointing point(thisTheta, thisPhi);
        int pixelNo = healpix->ang2pix(point);

        float thisWeight = errHexPixels[2][k]; // * eventWeight;  // I think eventWeight is just 1.0 for data, It is different for sim events

        map<int, float>::iterator thisEntry = hpWeights.find(pixelNo);
        if (thisEntry == hpWeights.end()) {
          vbprintf("  creating bin list entry for event %i %cPol bin %i  %8.3f \n", eventNumber, pol, pixelNo, thisWeight);
          hpWeights.insert(pair<int, float>(pixelNo, thisWeight));
        } else {
          vbprintf("  updating bin list entry for event %i %cPol bin %i  %8.3f \n",
                   eventNumber, pol, pixelNo, thisWeight);
          hpWeights[pixelNo] += thisWeight;
        }
      }
      for (pair<int, float> thisEntry : hpWeights) {
        //runNumber = runNumber;
        hpBin = thisEntry.first;
        hpWeight = thisEntry.second;
        //eventNumber = eventNumber;
        //pol = pol;
        vbprintf(" writing bin list entry for event %i %cPol bin %i  %8.3f \n",
                eventNumber, pol, hpBin, hpWeight);
        tree1->Fill();
      }
    }
  }


  tree1->BuildIndex("eventNumber", "pol");

  TTreeIndex* tree1Index = (TTreeIndex*)tree1->GetTreeIndex();

      




  // ------- Build the histograms by cut value for every continent bin put data structure into the map ------------------------------------
  for (int e=0; e<tree1Index->GetN(); ++e) {
  //for (int e=0; e<1; ++e) {
    int entNum = tree1Index->GetIndex()[e];
    //lprintf("entry number %i \n", entNum); 
    tree1->GetEntry(entNum);
    if (pol != whichPol) {continue;}
    //lprintf(" event number %i %cPol bin %i \n", eventNumber, pol, hpBin);
    //int p = (pol=='H' || pol=='L') ? 0 : 1;
    vbprintf("Creating Hists for event %i \n", eventNumber);
    if (eventNumber != prevEvent || pol != prevPol) {
      int getRes = tree0->GetEntryWithIndex(eventNumber, pol);
      useThisEntry = true;
      if (getRes <= 0) {lprintf("database error event %i pol %c bin %i - skipping \n", eventNumber, pol, hpBin); useThisEntry = false;}
      if (circPeakSepThreshold>0 && cPolDist>circPeakSepThreshold) {useThisEntry = false;}
      if (circPeakStrengthThreshold>0 && minCircPeakVal<circPeakStrengthThreshold) {useThisEntry = false;}
      //lprintf(" event number %i %cPol \n", eventNumber, pol) ; 
      //lprintf(" indexed tree get result is %i \n", getRes);

      // we trying something diff -Jacob
      //cutVal = -cutSlope*2*cPolPeak[1] + cohSnr2;
      cutVal = -cutSlope*cPeak + cohSnr2;
      //cutVal = -cutSlope*cPeak + cohSnr3;

      if (hpBin == ourBin) {vbprintf(" event number %i %cPol, bin=%i, cPeak=%.4f, hPeak=%.1f, cutval=%.4f \n", eventNumber, pol, hpBin, cPeak, hPeak, cutVal);}
      prevEvent = eventNumber;
      prevPol = pol;
    }
    if (!useThisEntry) {continue;}
    if (cutValHist.find(hpBin) == cutValHist.end()) {
      // instantiate new data row and histograms for this bin and put to the map
      dataRow_t thisRow; 
      //lprintf("  creating new histogram for %cPol bin %i \n", pol, hpBin);
      char title[128]; sprintf(title, "SNR-correlation cut value bin %i, run %i-%i, %cPol;y-intercept;count", hpBin, startRun, endRun, pol);
      char name[128]; sprintf(name, "cutValHist%c%05i", pol, hpBin);
      TH1F* thisHist = new TH1F(name, title, cutValBins, cutValMin, cutValMax);
      thisRow.cutHist = thisHist;

      sprintf(title, "SNR-correlation cut value bin %i with all events weight=1, run %i-%i, %cPol;y-intercept;count", hpBin, startRun, endRun, pol);
      sprintf(name, "cutValHistTotal%c%05i", pol, hpBin);
      TH1F* thisHistTotal = new TH1F(name, title, cutValBins, cutValMin, cutValMax);
      thisRow.cutHistTotal = thisHistTotal;

      sprintf(title, "SNR-correlation cut value bin %i with only weighted events, weight=1, run %i-%i, %cPol;y-intercept;count", hpBin, startRun, endRun, pol);
      sprintf(name, "cutValHistWeighted%c%05i", pol, hpBin);
      TH1F* thisHistWeighted = new TH1F(name, title, cutValBins, cutValMin, cutValMax);
      thisRow.cutHistWeighted = thisHistWeighted;

      sprintf(title, "SNR-correlation cut value bin %i for pseudo experiment, run %i-%i, %cPol;y-intercept;count", hpBin, startRun, endRun, pol);
      sprintf(name, "cutValHistPse%c%05i", pol, hpBin);
      TH1F* thisHistPseu = new TH1F(name, title, cutValBins, cutValMin, cutValMax);
      thisRow.cutHistPse = thisHistPseu;

      sprintf(title, "SNR vs correlation peak bin %i, run %i-%i, %cPol;corr Peak;snr", hpBin, startRun, endRun, pol);
      sprintf(name, "corrSnrHist%c%05i", pol, hpBin);
      TH2F* thisCorrPeakSnrHist = new TH2F(name, title, 50, 0, 0.5, 50, 0, 20); 
      thisRow.corrPeakSnrHist = thisCorrPeakSnrHist;

      sprintf(title, "SNR vs correlation peak bin %i, run %i-%i, %cPol (sim);corr Peak;snr", hpBin, startRun, endRun, pol);
      sprintf(name, "corrSnrSimHist%c%05i", pol, hpBin);
      TH2F* thisCorrPeakSnrSimHist = new TH2F(name, title, 50, 0, 0.5, 50, 0, 20); 
      thisRow.corrPeakSnrSimHist = thisCorrPeakSnrSimHist;

      sprintf(title, "SNR-correlation cut value bin %i, run %i-%i, %cPol (sim);y-intercept;count", hpBin, startRun, endRun, pol);
      sprintf(name, "simCutHist%c%05i", pol, hpBin);
      TH1F* thisSimHist = new TH1F(name, title, cutValSimBins, cutValMin, cutValMax);    //we want more presision in cutVal than 0.5, so increase num bins here -Jacob! 
      thisRow.simCutHist  = thisSimHist; 

      sprintf(title, "S and S_up, bin %i, %cPol (sim);y-intercept", hpBin, pol);
      sprintf(name, "simSigHist%c%05i", pol, hpBin);
      thisRow.sHist = new TH1F(name, title, cutValSimBins, cutValMin, cutValMax); 
      thisRow.sHist = new TH1F(name, title, cutValSimBins, cutValMin, cutValMax); 

      sprintf(title, "S_up bin %i, %cPol (sim);y-intercept", hpBin, pol);
      sprintf(name, "sUpHist%c%05i", pol, hpBin);
      thisRow.sUpHist = new TH1F(name, title, cutValSimBins, cutValMin, cutValMax); 
      
      sprintf(title, "S/S_up bin %i, %cPol (sim);y-intercept;S/S_up", hpBin, pol);
      sprintf(name, "s_sUpHist%c%05i", pol, hpBin);
      thisRow.s_sUpHist = new TH1F(name, title, cutValSimBins, cutValMin, cutValMax); 
            
      cutValHist.insert(pair<int, dataRow_t>(hpBin, thisRow));
      thisSimHist = 0;      
      thisHist = 0;

      //oindree deleting hists to make memory leak warnings go away
      delete thisHist, thisHistTotal, thisHistWeighted, thisHistPseu, thisCorrPeakSnrHist, thisCorrPeakSnrSimHist, thisSimHist; 
    } 
    //lprintf(" bin number %i, %.3f, %.4f \n", hpBin, hpWeight, cutVal);
    
    // If all of these are set up to cut any events above a yint of 25, 
    // If we want to include those we need to set up a quick redefinition if 
    // statement to catch the extra events.... not really sure how though.
    // come back to this if Amy thinks it matters alot.

    cutValHist[hpBin].numEvents += hpWeight;
    cutValHist[hpBin].cutHist->Fill(cutVal, hpWeight);
    cutValHist[hpBin].cutHistTotal->Fill(cutVal, 1.0);                                 // Added by Jacob!
    if (fabs(hpWeight - 1.0) > 1e-3)                                                   // so we know how many total events and weighted 
    {                                                                                  // events we have
      cutValHist[hpBin].cutHistWeighted->Fill(cutVal, 1.0);
      //printf("%0.16f\n", hpWeight);
      //cout << "WHHHYYYYYY!!!!!!! " << hpWeight << endl; 
    }                           
    cutValHist[hpBin].corrPeakSnrHist->Fill(cPeak, cohSnr2, hpWeight);
    corrSnrHist->Fill(cPeak, cohSnr2, hpWeight);
    //break;
  }  


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //  Calculate leakage between bins for each passing event
  //
  //   -Loop over Events
  //   --Get hpBin, payload position, and event positon for each event
  //   --Loop over grid in phi and theta
  //   ---Grid is 100 by 100 points (for now)
  //   ---Grid coveres (-5sigma) through (5sigma) in phi and theta
  //   ---Find Prob of event at grid point G reconstructing to event point E
  //   ---Find hpBin for G
  //   
  //   -Info is saved to spillover.
  //   --for each Bin, spillover saves its event contribution to other bins
  //   --For each bin, total prob is normalized to 1
  //   --When including probibility contributions that do not reconstruct to the earth
  //
  //   -If loadSpilloverFromFile is set to true, then instead of calculating the 
  //    spillover it pulls it from a default file location based on the input dir and 
  //    the orientation parameters
  //   --This is waaaay faster
  //   --I was lazy, so there is no protection logic agaisnt it not finding the filess.... sooo...
  //     you know... be careful.
  //
  ///////////////////////////////////////////////////////////////////////////////////////////////////


  // spillover explained!
  // (BinFrom, (BinTo, (Contribution, Count)))
  map<int,map<int,pair<double,int> > > spillover;

  if(!loadSpilloverFromFile){

    int nSteps = 100;

    // Pre calculate a vector of vectors of probablilities to save time later

    vector<vector<double> > probArray(100,vector<double>(100,0));

    for (int xi = 0; xi < nSteps; xi++)
      {
      double thisTh = 0.1*(50-(float)xi)*thError ;
      for(int yi = 0; yi < nSteps; yi++)
      {
        double thisPhi = 0.1*(50-(float)yi)*phiError;
        // Find P for this grid point
        double thisP = exp( -1*(thisPhi)*(thisPhi)/(2*phiError*phiError) - (thisTh)*(thisTh)/(2*thError*thError) );
        thisP *= 1/(2*M_PI*phiError*thError);       // normailize to 1
        probArray[xi][yi] = thisP;
      }
    }
  
    // define the minimum probablility a grid point can have for us to care about it.
    // this value is the hieght of a normalized gaussian at the 4 sigma point.
    double minProb = 0.00013383;
  
    map<int,map<int,int> > lastevent;


    // bedmap for bedmap stuff
    //BedmapReader* bedmap = BedmapReader::Instance(false);

    // for geom stuff? not sure we need this
    //AnitaGeomTool* geom = AnitaGeomTool::Instance();


    //loop over events
    //for (int e=0; e<tree1Index->GetN(); ++e) {
    for (int e=0; e<1; ++e) {
      int entNum = tree1Index->GetIndex()[e];
      tree1->GetEntry(entNum);
      if (pol != whichPol) {continue;}
      vbprintf("event %i in hpBin %i \n", eventNumber, hpBin);
      // do circPol cuts and get info if new event
      if (eventNumber != prevEvent || pol != prevPol) {
        int getRes = tree0->GetEntryWithIndex(eventNumber, pol);
        useThisEntry = true;
        if (getRes <= 0) {lprintf("database error event %i pol %c bin %i - skipping \n", eventNumber, pol, hpBin); useThisEntry = false;}
        if (circPeakSepThreshold>0 && cPolDist>circPeakSepThreshold) {useThisEntry = false;}
        if (circPeakStrengthThreshold>0 && minCircPeakVal<circPeakStrengthThreshold) {useThisEntry = false;}
        //lprintf(" event number %i %cPol \n", eventNumber, pol) ;
        //lprintf(" indexed tree get result is %i \n", getRes);
        prevEvent = eventNumber;
        prevPol = pol;
      }
    
      if (spillover.find(hpBin) == spillover.end())
      {
        map<int,pair<double,int> > thisMap;
        spillover.insert(pair<int,map<int,pair<double,int> > >(hpBin, thisMap) );
        map<int,int> thisMap2;
        lastevent.insert(pair<int,map<int,int> >(hpBin, thisMap2) );
      }
      //cout << " spillover[" << hpBin << "] is " << spillover[hpBin] << endl;


      // At this point we have the event and payload positions and phi and theta.
      //
      // ea = easting
      // no = northing
      // lat = latitude
      // lon = longitude 
      // gps = gps payload info
      //


      //UsefulAdu5Pat* gps = new UsefulAdu5Pat(gpsRaw);

      double eventTh = 0; 
      double eventPhi = 0;
    
      eventTh = theta;
      eventPhi = phi;

      // loop over x and y to make Grid in phi and theta
      //   let x and y be iterators for the angle from the payload
      //   x is the theta itterator
      //   y is the phi itterator
      int nSteps = 100;
      for (int xi = 0; xi < nSteps; xi++) 
      {
        //double thisTh = 0.1*(50-(float)xi)*thError + eventTh; 
        for(int yi = 0; yi < nSteps; yi++)
        {
          //double thisPhi = 0.1*(50-(float)yi)*phiError + eventPhi;
  
          // Find P for this grid point
        
          double thisP = probArray[xi][yi];
          if (thisP < minProb) {continue;}
          thisP *= hpWeight;

          // Find which bin this grid point falls in.
          int thisHpBin = 0;
        
          //  Trace back to continent to find lon and lat
          double thisLat = 0; double thisLon= 0; double thisAlt = 0; //double thisAdj = 0;
          //int traceStatus = gps->traceBackToContinent(thisPhi*M_PI/180.0, thisTh*M_PI/180.0, &thisLon, &thisLat, &thisAlt, &thisAdj);

          if (xi == 50 && yi == 50){
            //cout << "(thisTh,thisPhi): (" << thisTh << "," << thisPhi << ") (eventTh,eventPhi): (" << eventTh << "," << eventPhi << ")" << endl;
            //cout << "(thisLat,thisLon): (" << thisLat << "," << thisLon << ") (eventLat,eventLon): (" << lat << "," << lon << ")" << endl;
          }
          //cout << "traceback status " << traceStatus << endl;
          //cout << "M_PI = " << M_PI << endl;

          //  Note: If traceBack does not reconstruct to the earth then
          //          lon lat and/or alt will be set to -9999.
          //        If that happens set bin = 0 and say that prob 
          //          contribution is lost (event doesnt reconstruct to the earth)

          if (thisLat == -9999 || thisLon == -9999 || thisAlt == -9999) { thisHpBin = 0; }
          else if (thisLat == 0 && thisLon == 0 && thisAlt == 0) {thisHpBin = 0; }
          else {

            //  Turn lon and lat into earth polar coordinates 
            //   (Note: these are different theta and phi)

            double thisETheta = (-thisLat + 90.0 + hpThOffset) * M_PI/180.0;
            double thisEPhi = (thisLon + hpPhiOffset) * M_PI/180.0;
            while (thisEPhi < -M_PI) thisEPhi += M_PI;

            //  Turn Ephi and Etheta into a healpix bin number

            pointing thisPoint(thisETheta, thisEPhi);
            thisHpBin = healpix->ang2pix(thisPoint);
          }
          //cout << "Bin " << hpBin << " contributes " << thisP << " to Bin " << thisHpBin << endl;

          // Now for this gridpoint we have the probability P and its healpix bin.
          if (spillover[hpBin].find(thisHpBin) == spillover[hpBin].end())
          {
            pair<double,int> temppair = pair<double,int>(thisP,1);
            spillover[hpBin].insert(pair<int,pair<double,int> >(thisHpBin,temppair));
            lastevent[hpBin].insert(pair<int,int>(thisHpBin,eventNumber));
          } else {
            spillover[hpBin][thisHpBin].first += thisP;
            if (lastevent[hpBin][thisHpBin] != eventNumber)
            {
              spillover[hpBin][thisHpBin].second += 1;
              lastevent[hpBin][thisHpBin] = eventNumber;
            }
          }
        }
      }
    }

    // We still need to normalize the spillover per bin.
    //  -Each bin should have a total of 1 in all the bins in it's map.
    //  -While we are here, lets also make an output file.


    char spillFilename[1024];
    sprintf(spillFilename, "%s/spillover%c.txt", outputDir, whichPol);
    FILE*  spillFileOut = fopen(spillFilename, "w");
    ofstream spillFile(spillFilename, ios_base::app);
    fprintf(spillFileOut, "Giving Bin #, Recieving Bin #, percent spillover, number events spilling over \n");

    // Loop over bins
    for (map<int, map<int,pair<double,int> > >::iterator thisPair=spillover.begin(); thisPair!=spillover.end(); ++thisPair) {
      int hpBin = thisPair->first;
      map<int,pair<double,int> > thisMap = thisPair->second;
      float thisNorm = 0;
      // Loop over bins hpBin contributes to, to find normalization
      for ( map<int,pair<double,int> >::iterator thisCont=thisMap.begin(); thisCont!=thisMap.end(); ++thisCont) {
        double contribution = thisCont->second.first;
        thisNorm += contribution;
      }
      // Loop over bins hpBin contribues background to, to normalize contribution.
      for ( map<int,pair<double,int> >::iterator thisCont=thisMap.begin(); thisCont!=thisMap.end(); ++thisCont) {
        int inHpBin = thisCont->first;
        thisCont->second.first /= thisNorm;
        cout << "Bin " << hpBin << " contributes " << thisCont->second.first << " percent of its background to bin ";
        cout << inHpBin << " from " << thisCont->second.second << " events." << endl;
        fprintf(spillFileOut, "%5i, %5i, %6.4f, %5i \n", hpBin, inHpBin, thisCont->second.first, thisCont->second.second);
      }
    }  
    fclose(spillFileOut);
    //cin.ignore();
  }
  else // if we are loading the spillover data:
  {
    //this just a dumn directory thing
    if ( hpPhiOffset < 0.01 && hpThOffset > -0.01 )  // if both are zero
    {
      hpThOffset = -0.0;
    }
    else if (hpThOffset > -0.01)  // if just thoffset = 0
    {
      hpThOffset = 0.0;
    }     

    char spillDir[1024];
    sprintf(spillDir, "%s/spillover_skew/phiOff%03.4f/thOff%03.4f", inputDir, hpPhiOffset, hpThOffset);

    char spillFilename[1024];
    sprintf(spillFilename, "%s/spillover%c.txt", spillDir, whichPol);
    printf("spillover filename is %s \n", spillFilename);
    ifstream spillFile(spillFilename);
    string spillLine;
    while (!spillFile.eof()) {
      getline(spillFile, spillLine);
      //printf("%s \n", spillLine.c_str());
      if (spillLine.c_str()[0] != 'G') {
        int hpBin;
        int inHpBin;
        float cont;
        int numEvents;
        stringstream spillStr(spillLine);

        string field;
        getline( spillStr, field, ',' );
        stringstream fs(field);
        fs >> hpBin;
        getline( spillStr, field, ',' );
        stringstream fs1(field);
        fs1 >> inHpBin;
        getline( spillStr, field, ',' );
        stringstream fs2(field);
        fs2 >> cont;
        getline( spillStr, field, ',' );
        stringstream fs3(field);
        fs3 >> numEvents;

        printf("bin #%i  inBin=%i  cont=%f  numEvents=%i \n", hpBin, inHpBin, cont, numEvents);

        // put into map

        if (spillover.find(hpBin) == spillover.end())
        {
          map<int,pair<double,int> > thisMap;
          spillover.insert(pair<int,map<int,pair<double,int> > >(hpBin, thisMap) );
        }
        if (spillover[hpBin].find(inHpBin) == spillover[hpBin].end())
        {
          pair<double,int> temppair = pair<double,int>(cont,numEvents);
          spillover[hpBin].insert(pair<int,pair<double,int> >(inHpBin,temppair));
        } else {
          spillover[hpBin][inHpBin].first += cont;
          spillover[hpBin][inHpBin].second += numEvents;
        }
      }
    }
    spillFile.close();
  }

  

  // draw the diagonal cut plot
  TCanvas* corrSnrCanv = new TCanvas("corrSnrCanv", "corrSnrCanv", 600, 600);
  corrSnrCanv->cd(0);
  corrSnrHist->Draw("COLZ");
  corrSnrCanv->Draw();
  // draw the cutlines
  float xMin = gPad->GetUxmin();
  float xMax = gPad->GetUxmax();
  float yMin = gPad->GetUymin();
  float yMax = gPad->GetUymax();  
  float deltaB = (float)(cutValMax-cutValMin)/(float)cutValBins;
  //lprintf("xMin=%f, yMin=%f, xMax=%f, yMax=%f, dB=%f \n", xMin, yMin, xMax, yMax, deltaB);
  TLine* cutLine = 0;
  for (float thisB = yMin; cutSlope*xMax + thisB < yMax; thisB += deltaB) {
    float x1 = max(xMin, (yMax-thisB) / cutSlope);
    float x2 = min(-thisB/cutSlope, xMax);
    float y1 = cutSlope * x1 + thisB;
    float y2 = cutSlope * x2 + thisB;
    cutLine = new TLine(x1, y1, x2, y2);
    cutLine->SetLineColor(kRed);
    cutLine->Draw();
    cutLine = 0;   // just forget about it and let our grandchildren clean it up.
  }
  //app->Run();
  //return 0;
  //gMinuit->SetPrintLevel(-1);

  TCanvas* dummyCanv = new TCanvas("dummy", "dummy", 0, 0);
  dummyCanv->cd();
  
  // ---------  iterate through each continent bin in the map and do log likelihood fit on its histogram --------------------------------------
  for (map<int, dataRow_t>::iterator thisPair=cutValHist.begin(); thisPair!=cutValHist.end(); ++thisPair) {      
    // find the bins containing the maximum value and last nonzero value
    //if (thisPair->first != ourBin) continue;
    printf("healpix bin %i \n", thisPair->first);

    /*
    //in this code we only want to look at bins Oindree set the background to 0.1 for.  keep only those bins
    if (   thisPair->first != 3012 
        && thisPair->first != 3018
        && thisPair->first != 2968
        && thisPair->first != 2997
        && thisPair->first != 3019
        && thisPair->first != 3017
        && thisPair->first != 2935
        && thisPair->first != 3004
        && thisPair->first != 2998
        && thisPair->first != 2996
        && thisPair->first != 3007
        && thisPair->first != 3021
        && thisPair->first != 2979
        && thisPair->first != 2999
        && thisPair->first != 2951
        && thisPair->first != 2973
        && thisPair->first != 2952 ) 
    {
      thisPair->second.status = -7;    // we dont care right now...
      thisPair->second.binStatus = 1;
      continue;
    }
    */

    int thisHpBin = thisPair->first; 

    float thisCutVal = 0; float thisIntercept = 0; float thisSlope = 0; 

    // these are bins with naturally calculated cutvals.
    if ( thisHpBin == 3014 ) { thisCutVal =  11.0; thisIntercept = 7.710; thisSlope = -1.119; }   //
    if ( thisHpBin == 3013 ) { thisCutVal =  10.0; thisIntercept = 7.668; thisSlope = -1.185; }   //
    if ( thisHpBin == 3037 ) { thisCutVal =  12.3; thisIntercept = 4.894; thisSlope = -0.7610; }   //
    if ( thisHpBin == 3015 ) { thisCutVal =  10.9; thisIntercept = 9.854; thisSlope = -1.389; }   //
    if ( thisHpBin == 3016 ) { thisCutVal =   9.4; thisIntercept = 9.685; thisSlope = -1.546; }   //
    if ( thisHpBin == 2971 ) { thisCutVal =  11.1; thisIntercept = 8.001; thisSlope = -1.2100; }   //
    if ( thisHpBin == 2970 ) { thisCutVal =  12.0; thisIntercept = 7.759; thisSlope = -1.113; }   //
    if ( thisHpBin == 2936 ) { thisCutVal =   9.3; thisIntercept = 9.715; thisSlope = -1.461; }   //
    if ( thisHpBin == 2989 ) { thisCutVal =   9.9; thisIntercept = 9.955; thisSlope = -1.515; }   //
    if ( thisHpBin == 2937 ) { thisCutVal =  10.2; thisIntercept = 7.238; thisSlope = -1.119; }   //
    if ( thisHpBin == 2990 ) { thisCutVal =  10.0; thisIntercept = 8.117; thisSlope = -1.212; }   //
    if ( thisHpBin == 2988 ) { thisCutVal =  10.1; thisIntercept = 8.826; thisSlope = -1.439; }   //
    if ( thisHpBin == 2991 ) { thisCutVal =   8.7; thisIntercept = 9.920; thisSlope = -1.725; }   //
    if ( thisHpBin == 3029 ) { thisCutVal =  10.9; thisIntercept = 8.694; thisSlope = -1.289; }   //
    if ( thisHpBin == 3003 ) { thisCutVal =  14.4; thisIntercept = 6.901; thisSlope = -0.8453; }   //
    if ( thisHpBin == 3030 ) { thisCutVal =   9.9; thisIntercept = 8.427; thisSlope = -1.301; }   //
    if ( thisHpBin == 2938 ) { thisCutVal =  10.3; thisIntercept = 5.842; thisSlope = -1.099; }   //
    if ( thisHpBin == 3011 ) { thisCutVal =   9.9; thisIntercept = 9.734; thisSlope = -1.561; }   //
    if ( thisHpBin == 3010 ) { thisCutVal =   9.4; thisIntercept = 10.44; thisSlope = -1.727; }   //
    if ( thisHpBin == 3008 ) { thisCutVal =  10.4; thisIntercept = 8.767; thisSlope = -1.331; }   //
    if ( thisHpBin == 2939 ) { thisCutVal =  12.1; thisIntercept = 6.908; thisSlope = -1.107; }   //
    if ( thisHpBin == 2901 ) { thisCutVal =  10.8; thisIntercept = 5.470; thisSlope = -0.9245; }   //
    if ( thisHpBin == 2977 ) { thisCutVal =  11.2; thisIntercept = 7.163; thisSlope = -1.041; }   //
    //if ( thisHpBin == 3001 ) { thisCutVal =  11.4; }   //  sideband
    //if ( thisHpBin == 2978 ) { thisCutVal =  10.0; }   //  sideband
    //if ( thisHpBin == 2900 ) { thisCutVal =  13.7; }   //  sideband


    // these are bins with adjusted cutVals
    if ( thisHpBin == 3012 ) { thisCutVal =  12.5; thisIntercept = 7.439; thisSlope = -1.011; }   //12.5
    if ( thisHpBin == 3018 ) { thisCutVal =  10.1; thisIntercept = 7.292; thisSlope = -1.226; }   //10.1
    if ( thisHpBin == 2968 ) { thisCutVal =  16.2; thisIntercept = 3.418; thisSlope = -0.5673; }   //16.2
    if ( thisHpBin == 2997 ) { thisCutVal =   9.3; thisIntercept = 7.455; thisSlope = -1.339; }   //9.3
    if ( thisHpBin == 3019 ) { thisCutVal =  13.5; thisIntercept = 4.971; thisSlope = -0.7723; }   //13.5
    if ( thisHpBin == 3017 ) { thisCutVal =  14.6; thisIntercept = 4.248; thisSlope = -0.6751; }   //14.6
    if ( thisHpBin == 2935 ) { thisCutVal =  14.3; thisIntercept = 4.000; thisSlope = -0.6749; }   //14.3
    if ( thisHpBin == 3004 ) { thisCutVal =  11.7; thisIntercept = 8.691; thisSlope = -1.174; }   //11.7
    if ( thisHpBin == 2998 ) { thisCutVal =  12.2; thisIntercept = 4.977; thisSlope = -0.8474; }   //12.2
    if ( thisHpBin == 2996 ) { thisCutVal =  57.4; thisIntercept = -0.7038; thisSlope = -0.116; }   //57.4
    if ( thisHpBin == 3007 ) { thisCutVal =  11.2; thisIntercept = 8.921; thisSlope = -1.248; }   //11.2
    if ( thisHpBin == 3021 ) { thisCutVal =  26.3; thisIntercept = 5.424; thisSlope = -0.4362; }   //26.3
    if ( thisHpBin == 2979 ) { thisCutVal =  13.1; thisIntercept = 4.557; thisSlope = -0.7685; }   //13.1
    if ( thisHpBin == 2999 ) { thisCutVal =  30.8; thisIntercept = 4.871; thisSlope = -0.3605; }   //30.8
    if ( thisHpBin == 2951 ) { thisCutVal = 168.6; thisIntercept = -0.6156; thisSlope = -0.05241; }   //168.6
    if ( thisHpBin == 2973 ) { thisCutVal =  11.5; thisIntercept = 7.846; thisSlope = -1.131; }   //11.5
    if ( thisHpBin == 2952 ) { thisCutVal = 114.0; thisIntercept = -0.7287; thisSlope = -0.06341; }   //114.0

    if (thisCutVal == 0) { 
      thisPair->second.status = -7;
      thisPair->second.binStatus = 1;
    }
    
    thisPair->second.optCutVal = thisCutVal;


    printf("healpix bin %i \n", thisPair->first);

    dataRow_t& thisRow = thisPair->second;
    TH1F*& thisHist = thisRow.cutHist;
    if (thisHist->GetSumOfWeights() > 0) {
      int maxBin = thisHist->GetMaximumBin();
      int startBin = maxBin+2;
      int lastNzBin = -1;
      int nzBinCount = 0;
      for (int k=thisHist->GetNbinsX(); k>=0; --k) {
        if (thisHist->GetBinContent(k) > 0) {
          if (lastNzBin < 0) {lastNzBin = k;}
          ++nzBinCount;
        }
      }

      ////////////////// OINDREE WORKING ON CHANGING THIS ////////////////////////////////////////////
      //int lastDescBin = -1; bool done = false;
      //for (int k=startBin; k<=thisHist->GetNbinsX() && !done; ++k) {
        //if (thisHist->GetBinContent(k) < thisHist->GetBinContent(k-1) && thisHist->GetBinContent(k)>0) {
        //if (thisHist->GetBinContent(k)>0.0001) {
          //lastDescBin = k;
        //} else {
          //done = true;
        //}
      //}
      //////////////////////////////////////////////////////////////////////////////////////////////////

      double hist_sum_content = 0.0; 
      int descBinCount = 0; 
    
      for (int k = startBin; k <= thisHist->GetNbinsX(); ++k) {
        hist_sum_content = hist_sum_content + thisHist->GetBinContent(k);
        if (thisHist->GetBinContent(k) > 0.0001) {
          descBinCount = descBinCount + 1; 
        }
      }

      //int descBinCount = lastDescBin - maxBin;
      //int descBinCount = lastDescBin - (startBin-1);
      //vbprintf("  hp-bin %i:  maxBin=%i,   lastDecBin=%i, descBinCount=%i   \n", thisPair->first, maxBin, lastDescBin, descBinCount);
      vbprintf("  hp-bin %i:  maxBin=%i,   GetNbinsX()=%i, descBinCount=%i   \n", thisPair->first, maxBin, thisHist->GetNbinsX(), descBinCount);
      thisHist->Sumw2();
      if (descBinCount > 1) { 
      //if (true) { //this is oindree trying to fit every bin 

        // how we were fitting:  (Jacob!)               
        //int fitStatus = thisHist->Fit("expo", fitOpt, "", thisHist->GetBinLowEdge(startBin), thisHist->GetBinLowEdge(lastNzBin+1));
        // How we are fitting now using Minuit:
        //  We fit to both an exponetial and a power law
        //  The exponential is used for our actual fit, and our 
        //   pval calc and background and most stuff
        //  The power law is used for estimating our fit uncertainty.
        // Make out data vectors
        vector<double> pos;
        vector<double> meas;
        vector<double> var;

        //float thisHistSpacing = (cutValMax-cutValMin)/cutValBins;
        for (int iii=1; iii<thisHist->GetSize()-1; iii++)
        {
          double pos_here = thisHist->GetBinCenter(iii);
          if ( thisHist->GetBinLowEdge(startBin) <= pos_here && pos_here <= thisHist->GetBinLowEdge(lastNzBin+1)+trailingZeros*cutValBins )
          {
            pos.push_back(pos_here);
            meas.push_back(thisHist->GetBinContent(iii));
            var.push_back(thisHist->GetBinError(iii));
          }
        }
        // define our FCN for minuit (see fitFunc.h/.cc)
        FitFCN theFCN(meas, pos, var);
        PowerFCN powerFCN(meas, pos, var);
        
        // define fit parameters
        //  name of paramter, starting value, initial error
        // fit(x_i) = exp(a+b*x_i)
        // powFit(x_i) = c*x_i^d

        MnUserParameters upar;
        MnUserParameters uparPow;

        upar.Add("a", thisIntercept, 0.1 );
        upar.Add("b", thisSlope, 0.1 );
        //upar.Add("a", 10.0, 0.1 );  
        //upar.Add("b", -1.0, 0.1 );
        uparPow.Add("c", 10000,0.1 );
        uparPow.Add("d", -5.0,0.1 );

        MnMigrad migrad(theFCN,upar);
        MnMigrad migradPow(powerFCN,uparPow);

        cout << "Fitting HpBin " << thisHpBin << endl;
        // Lets make sure we get good fits by fitting 
        // the parameters individually then fitting them
        // together
        // expo fit first
        /*
        migrad.Fix("a");
        FunctionMinimum mmin = migrad();
        migrad.Release("a");
        migrad.Fix("b");
        mmin = migrad();
        migrad.Release("b");
        mmin = migrad();
        */
        migrad.Fix("a");
        migrad.Fix("b");
        FunctionMinimum mmin = migrad();

        // now the pow fit
        migradPow.Fix("c");
        FunctionMinimum mminPow = migradPow();
        migradPow.Release("c");
        migradPow.Fix("d");
        mminPow = migradPow();
        migradPow.Release("d");
        mminPow = migradPow();

        cout << " Retriving Contours" << endl;
        MnContours mcont(theFCN,mmin);

        // for root, fitstatus = 0 means its all good.
        // for minuit, fitstatus = 1 means its all good, 
        // so we add this little if after the fitStatus gets set. 
        int fitStatus;
        if (mmin.IsValid() == 1) { fitStatus = 0; }
        if (mmin.IsValid() == 0) { fitStatus = 1; }  //?? 

        //lprintf("fit status = %i \n", fitStatus);
        float minX = thisHist->GetBinLowEdge(startBin);
        float maxX = thisHist->GetBinLowEdge(lastNzBin+1);
        //float minY = thisHist->GetBinContent(lastNzBin);
        //float maxY = thisHist->GetBinContent(maxBin+1);
        double l0; double  edm, errdef, l0Pow;
        //int nvpar, nparx;

        // Change in how to get data from fit (using minuit now) Jacob!
        //TVirtualFitter::GetFitter()->GetStats(l0, edm, errdef, nvpar, nparx);
        l0 = mmin.Fval();
        l0Pow = mminPow.Fval();
        edm = mmin.Edm();
        errdef = mmin.Up();
        //nvpar = 2 //number of function variables. its 2. I dont know how to call it
        //nparx = > 7?  //total number of variables/ its 2+num bins, I dont know how to call it, but its not used anywhere. 

        vbprintf("    l0 =  %f \n", l0);    
        //thisPair->second.logL = l0/thisHist->GetSumOfWeights();
        thisPair->second.logLPow = l0Pow;
        thisPair->second.logL = l0;
        thisRow.endBin = lastNzBin;

        // this is also different now with minuit Jacob!
        //float thisSlope = TVirtualFitter::GetFitter()->GetParameter(1);       
        //float thisIntercept = TVirtualFitter::GetFitter()->GetParameter(0);
        //float thisSlope = mmin.UserParameters().Value(1);
        //float thisIntercept = mmin.UserParameters().Value(0);
      
        thisSlope = mmin.UserParameters().Value(1);
        thisIntercept = mmin.UserParameters().Value(0);
        float thisSlopeError = mmin.UserParameters().Error(1);
        float thisInterceptError = mmin.UserParameters().Error(0);
        // also get the cov matrix parameters
        float thisCovM00 = 0; //mmin.UserCovariance()(0,0);
        float thisCovM01 = 0; //mmin.UserCovariance()(0,1);
        float thisCovM11 = 0; //mmin.UserCovariance()(1,1);

        // get the pow fit
        float thisPowC = mminPow.UserParameters().Value(0);
        float thisPowD = mminPow.UserParameters().Value(1);

        //cout << " the cov matrix is:" << endl;
        //cout << "  | " << thisCovM00 << " " << thisCovM01 << " |" << endl;
        //cout << "  | " << thisCovM01 << " " << thisCovM11 << " |" << endl;
        //sleep(1);

        // use function i built to find the cov matrix?
        //FitFunc ff;
        //double  M11,M12,M21,M22;
        //ff.findCovM(M11,M12,M21,M22,pos,meas,thisIntercept,thisSlope);

        //cout << " the hesse matrix is (wrong):" << endl;
        //cout << "  | " << M11 << " " << M21 << " |" << endl;
        //cout << "  | " << M12 << " " << M22 << " |" << endl;
        //sleep(1);

        vector<pair<double,double> > oval;
        oval = mcont(0,1,40);
        /*
        for (int ic=0; ic<oval.size(); ic++)
        {
          cout << "{" << oval[ic].first << ", " << oval[ic].second << "}," << endl;
        }
        */
        // from our contour we need to calc the following
        //  - slopeAvg and interAvg
        //  - slopeSigma and interSigma
        //  - innerSigma (it goes with the intercept)
        //  - row (the corrolation between the parameters)
        double slopeAvg = 0;
        double interAvg = 0;        
        for (int ic=0; ic<oval.size(); ic++)
        {
          slopeAvg += oval[ic].second;
          interAvg += oval[ic].first;
        }
        slopeAvg /= oval.size();
        interAvg /= oval.size(); 
        float slopeSigma = 0;
        float interSigma = 0;
        float innerSigma = 0;
        float mindiff = 10000;
        for (int ic=0; ic<oval.size(); ic++)
        {
          if (abs(oval[ic].first-interAvg) >= interSigma)
          {
            interSigma = abs(oval[ic].first-interAvg);
          }
          if (abs(oval[ic].second-slopeAvg) >= slopeSigma)
          {
            slopeSigma = abs(oval[ic].second-slopeAvg);
          }
          if (abs(oval[ic].second-slopeAvg) <= mindiff)
          {
            mindiff = abs(oval[ic].second-slopeAvg);
            innerSigma = abs(oval[ic].first-interAvg);
          }
        }
        float row = pow(1-innerSigma*innerSigma/interSigma/interSigma,0.5);

        // also save the limiting cases of the contour, like, basically the 4 corners, for a later plotting thing

        float contourMinX = 0; float contourMinY = 0; float contourMaxX = 0; float contourMaxY = 0;
        float maxRad = 0;
        float minRad = 1000;
        for (int ic=0; ic<oval.size();ic++)
        {
          float thisRad = (oval[ic].first-interAvg)*(oval[ic].first-interAvg)+(oval[ic].second-slopeAvg)*(oval[ic].second-slopeAvg);
          if ( maxRad <= thisRad )
          {
            maxRad = thisRad;
            contourMaxX = abs(oval[ic].first-interAvg);
            contourMaxY = abs(oval[ic].second-slopeAvg);
          }
          if ( minRad >= thisRad )
          {
            minRad = thisRad;
            contourMinX = abs(oval[ic].first-interAvg);
            contourMinY = abs(oval[ic].second-slopeAvg);
          }          
        }
        // if the contour is empty or failed or whatever we need a backup. 
        if (innerSigma == 0 || interSigma == 0 || isnan(row))
        {
          interSigma = thisInterceptError ;
          slopeSigma = thisSlopeError;
          interAvg = thisIntercept;
          slopeAvg = thisSlope;
          row = pow(1-innerSigma*innerSigma/interSigma/interSigma,0.5);
        }
        thisRow.fitContourMax.first = contourMaxX;
        thisRow.fitContourMax.second = contourMaxY;
        thisRow.fitContourMin.first = contourMinX;
        thisRow.fitContourMin.second = contourMinY;
        // now save the parameters
        thisRow.fitErrorSlopeAvg = slopeAvg;
        thisRow.fitErrorInterceptAvg = interAvg;
        thisRow.fitErrorSlopeSigma = slopeSigma;
        thisRow.fitErrorInterceptSigma = interSigma;
        thisRow.fitErrorRow = row;
        thisRow.fitSlope = thisSlope;
        thisRow.fitIntercept = thisIntercept;
        thisRow.fitSlopeError = thisSlopeError;
        thisRow.fitInterceptError = thisInterceptError;
        thisRow.fitPowScale = thisPowC;
        thisRow.fitPowOrder = thisPowD;
        thisRow.fitSlopeM = thisCovM11;
        thisRow.fitInterceptM = thisCovM00;
        thisRow.fitCrossM = thisCovM01;
        thisRow.status = fitStatus;
        thisRow.xMin = minX;
        thisRow.xMax = maxX;
        thisRow.yMin = exp(minX*thisSlope + thisIntercept);
        thisRow.yMax = exp(maxX*thisSlope + thisIntercept);
        thisRow.startBin = startBin;
        if (thisRow.status == 1) {
          thisRow.binStatus = 2;
        }
        if (thisRow.status == 0 && thisRow.fitSlope >= 0) {
          thisRow.status = -2;
          thisRow.binStatus = 2;
        }
        vbprintf("     slope=%f, intercept=%f \n", thisRow.fitSlope, thisRow.fitIntercept);
        vbprintf("     logL=%f \n", thisPair->second.logL);

        TF1* myexpo = new TF1("myexpo","exp([0]+[1]*x)");
        myexpo->FixParameter(0,thisIntercept);
        myexpo->FixParameter(1,thisSlope);
        int fitStatus_fix = thisHist->Fit("myexpo", fitOpt, "", thisRow.xMin, thisRow.xMax);
      } // if  
      if (descBinCount < 5 || hist_sum_content < 5.0) {  
      //if (descBinCount < 4 || hist_sum_content < 4.0) { // tried loosen for hpol -- oindree 
        thisRow.status = -1;    // short on fit bins
        thisRow.binStatus = 1;
      }
    }   // if this hist get sum of weights > 0  
    else {
      thisRow.status = -2;     // empty histogram   
      thisRow.binStatus = 1;
    }
    //delete myexpo; 
  } //for 

  //for (int p=0; p<2; ++p) {for (pair<int, dataRow_t> thisPair : cutValHist) {thisPair.second.cutHist->Delete();}}
  //return 0;

    //for (pair<int, TH1F*> thisPair : cutValHist) {
  //   
  
  //  -----  draw and save histograms only: don't put any analysis processing in this loop  --------------------------------------------------
  for (pair<int, dataRow_t> thisPair : cutValHist) {
    dataRow_t& thisRow = thisPair.second;
    //if (thisRow.status == 0) {  //oindree commenting out, were we able to fit this? 
    if (true) { //oindree -- do this anyway
      if (thisPair.first == ourBin || !drawIt) {
        vbprintf(" bin %i, status=%i, slope=%f, intercept=%f, logL=%f, xMin=%f, xMax=%f, yMin=%f, yMax=%f \n", 
                thisPair.first, thisRow.status, thisRow.fitSlope, thisRow.fitIntercept, thisRow.logL,
                thisRow.xMin, thisRow.xMax, thisRow.yMin, thisRow.yMax);
        TH1F*& thisHist = thisPair.second.cutHist;
        //theHist = thisHist;
        vbprintf("hist pointer = %p  n=%f \n", thisPair.second.cutHist, thisPair.second.cutHist->GetSumOfWeights());
        vbprintf("hist pointer = %p  n=%f \n", thisHist, thisHist->GetSumOfWeights());
        vbprintf("  drawing diff plot for bin %i, \n", thisPair.first);
        TCanvas* thisCanv =  new TCanvas(thisHist->GetName(), thisHist->GetName(), 450, 450);  
        //canvi.push_back(thisCanv);
        TLegend oindreeLegend(0.6, 0.65, 0.4, 0.5); 
        oindreeLegend.SetTextSize(0.04);
        oindreeLegend.AddEntry(thisHist,Form("fitIntercept = %f",thisRow.fitIntercept),"p");
        oindreeLegend.AddEntry(thisHist,Form("fitSlope = %f",thisRow.fitSlope),"p");
        oindreeLegend.AddEntry(thisHist,Form("fitInterceptError = %f",thisRow.fitInterceptError),"p");
        oindreeLegend.AddEntry(thisHist,Form("fitSlopeError = %f",thisRow.fitSlopeError),"p"); 
        thisCanv->cd(0);
        gPad->SetLogy();
        gStyle->SetOptStat("");
        gStyle->SetOptFit(2);
        gStyle->SetLegendBorderSize(0); 
        thisHist->Draw();
        oindreeLegend.Draw(); 
        //singleCutValCanv->Draw();
        vbprintf("hist pointer = %p  n=%f \n", thisHist, thisHist->GetSumOfWeights());
        char outFilename[1024]; sprintf(outFilename, "%s/%s.png", outputDir, thisCanv->GetName());
        thisCanv->SaveAs(outFilename);
        sprintf(outFilename, "%s/%s.root", outputDir, thisCanv->GetName());
        thisCanv->SaveAs(outFilename);
        thisCanv = 0;
      }
    }
  }

  //if (drawIt) app->Run();
  //return 0;
  //if (doPValFlatness) {

  // ------------ do the pseudo-experiments on each continent bin and record the log likelihoods  ------------------------------------------

  for (map<int, dataRow_t>::iterator thisPair=cutValHist.begin(); thisPair!=cutValHist.end(); ++thisPair) {
    //lprintf("bin %i \n", thisPair->first);
    //if (thisPair.first == ourBin) {      
    dataRow_t& thisRow = thisPair->second;
    if (thisRow.status == 0) {
      vbprintf("Doing pseudo-experiments for bin %i \n", thisPair->first);
      vbprintf(" bin logL=%f \n", thisRow.logL);
      // do the pseudo-experiments
      // TODO instead of instantiating the histogram, set up an array
      vector<double> logLikeHistData(0);
      vector<double> logLikeHistPowData(0);
      //TH1F* thisLogLikeHist = new TH1F("logLikeHist", "pseudo-exp log likelihood", 1000, 0, 1);            
      float uMin = exp(thisRow.fitSlope * thisRow.xMin);
      float uMax = exp(thisRow.fitSlope * thisRow.xMax);
      lprintf(" fitSlope=%f, fitIntercept=%f, xMin=%f, xMax=%f, uMin=%f, uMax=%f \n", thisRow.fitSlope, thisRow.fitIntercept, thisRow.xMin, thisRow.xMax, uMin, uMax);

      //fixed for Brian's weighted method (Jacob!)

      // this is the number of events including weighted events
      float numEventsPse = thisRow.cutHist->Integral(thisRow.startBin, thisRow.endBin);
      float numEventsPseTotal = thisRow.cutHistTotal->Integral(thisRow.startBin, thisRow.endBin);
      float numEventsPseWeighted = thisRow.cutHistWeighted->Integral(thisRow.startBin, thisRow.endBin);
      float percentWeight = numEventsPseWeighted/numEventsPseTotal;
      vbprintf("numEvent: %f, numEventsTotal: %f, numEventsWeighted: %f \n", numEventsPse,numEventsPseTotal,numEventsPseWeighted);      

      thisRow.sampleInFitRange = numEventsPse;
      thisRow.sampleInFitRangeTotal = numEventsPseTotal;
      thisRow.sampleInFitRangeWeighted = numEventsPseWeighted;
      

      // make a new TF1 object with our fit.

      char functionChar [1024];
      sprintf(functionChar,"exp(%1.9f + %1.9f * x)",thisRow.fitIntercept,thisRow.fitSlope);
      cout << "function: " << functionChar << endl;
      TF1 *f1 = new TF1("f1",functionChar,thisRow.xMin,25);  // 1000);
      sprintf(functionChar,"%1.9f*pow(x,%1.9f)",thisRow.fitPowScale,thisRow.fitPowOrder);
      cout << "pow function: " << functionChar << endl;
      TF1 *f2 = new TF1("f2",functionChar,thisRow.xMax,25); 
      // if we want to increase this back up, we need to make all the hists bigger.  
      // 
      f1->SetNpx(10000);
      f2->SetNpx(10000);
      //f1->SetNormalized(1); //do we need to turn normailization on?...


      int pseSize = numEventsPseTotal;
      vbprintf(" [           pseudoExp percent completion          ] \n");
      for (int pe = 0; pe<pseudoExpCount; pe++) {
        // Jacob!
        if (pe == 0) { vbprintf(" ["); fflush(stdout); }
        else if ( floor(pe/(0.02*pseudoExpCount)) == pe/(0.02*pseudoExpCount) ) { vbprintf("="); fflush(stdout); }

        char histName[128]; sprintf(histName, "cutValHistPse%c%05i", pol, thisPair->first);
        char histTitle[128];
        sprintf(histTitle, "SNR-correlation cut value bin %i for pseudo experiment, run %i-%i, %cPol;y-intercept;count", thisPair->first, startRun, endRun, pol);
        TH1F* thisPseHist = new TH1F(histName, histTitle, cutValBins, cutValMin, cutValMax);

        sprintf(histName, "cutValHistPsePow%c%05i", pol, thisPair->first);
        sprintf(histTitle, "SNR-correlation cut value bin %i for pseudo experiment with Pow fit, run %i-%i, %cPol;y-intercept;count", thisPair->first, startRun, endRun, pol);
        TH1F* thisPsePowHist = new TH1F(histName, histTitle, cutValBins, cutValMin, cutValMax);

        //TH1F*& thisPseHist = thisRow.cutHistPse;
        //lprintf("  pseudo-exp %i, size=%i \n", pe, pseSize);         
        thisPseHist->Sumw2();
        thisPsePowHist->Sumw2();

        vector<float> xR;
        vector<float> powR;
        vector<float> wR;
        for (int e=0; e<pseSize; ++e) {

          /* 
          // This method works pretty spicifically just for 
          // an exonential fit, so if you change the function! dont use this!
          // also, we aren't even using it... so. just use roots built in distrobution stuff.

          // throw a random number in [uMin,uMax]
          float yR = gRandom->Rndm();
          yR = uMin + (uMax- uMin)*yR;
          // look it up on the inverse CDF
          float xR = (yR - thisRow.fitIntercept)/ thisRow.fitSlope;
          //xR = floor(xR);
          */

          // use TF1 to get psudo data
          xR.push_back( f1->GetRandom() );
          powR.push_back( f2->GetRandom() );

          //cout << xR[e] << " " << powR[e] << endl;
          //sleep(0.4);

          // take into account the weighted events:
          float isWeight = gRandom->Rndm();
          if (isWeight < percentWeight)
          { wR.push_back( gRandom->Rndm() ) ; }
          else
          { wR.push_back( 1.0 ); }

          // in addition to random x and w, we want random slope and intercept

          //lprintf("       filling pseudo-exp hist with value %f, weight %f \n",xR,wR);
        }

        float xRMax = *max_element(xR.begin(),xR.end());
        if ( xRMax > 25 )
        {
          float newNBins = (floor(xRMax)+1)*2;
          float newBinMin = 0;
          float newBinMax = floor(xRMax)+1;
          cout << "new bin layout: " << newNBins << ", " << newBinMin << ", " << newBinMax << endl;
          // if true, we need to resize our hists... which apperently means just making a new hist.
          thisPseHist = new TH1F(thisPseHist->GetName(),thisPseHist->GetTitle(),newNBins,newBinMin,newBinMax);
          if(pe == 100) 
          {
            thisRow.cutHistPse = new TH1F(thisRow.cutHistPse->GetName(),thisRow.cutHistPse->GetTitle(),newNBins,newBinMin,newBinMax);  
          }
        }
        float powRMax = *max_element(powR.begin(),powR.end());
        if ( powRMax > 25 )
        {
          float newNBins = (floor(powRMax)+1)*2;
          float newBinMin = 0;
          float newBinMax = floor(powRMax)+1;
          cout << "new bin layout: " << newNBins << ", " << newBinMin << ", " << newBinMax << endl;
          // if true, we need to resize our hists... which apperently means just making a new hist.
          thisPsePowHist = new TH1F(thisPsePowHist->GetName(),thisPsePowHist->GetTitle(),newNBins,newBinMin,newBinMax);
        }

        for ( int fill_i=0; fill_i < xR.size(); fill_i++ )
        { 
          thisPseHist->Fill(xR[fill_i],wR[fill_i]);
          if (pe == 100) { thisRow.cutHistPse->Fill(xR[fill_i],wR[fill_i]); }
        }
        for ( int fill_i=0; fill_i < powR.size(); fill_i++ )
        {
          thisPsePowHist->Fill(powR[fill_i],wR[fill_i]);
          //if (pe == 100) { thisRow.cutHistPse->Fill(powR[fill_i],wR[fill_i]); }
        }


        dummyCanv->cd();
        int fitStatus = 0;

        // we use minuit2 now! changes!! (Jacob!)

        //TF1* myexpo = new TF1("myexpo","exp([0]+[1]*x)");
        //myexpo->FixParameter(0,thisRow.fitIntercept);
        //myexpo->FixParameter(1,thisRow.fitSlope);
        //fitStatus = thisPseHist->Fit("myexpo", fitOpt, "", thisRow.xMin, thisRow.xMax);

        //thisRow.xMaxPse = thisRow.xMax;                                                        //original fit range?
        thisRow.xMaxPse = thisPseHist->GetBinLowEdge( 1+thisPseHist->FindLastBinAbove(0,1) );  //or range for new pse data?
        thisRow.xMaxPowPse = thisPsePowHist->GetBinLowEdge( 1+thisPsePowHist->FindLastBinAbove(0,1) );

        // Make out data vectors
        vector<double> pos_pse, meas_pse, var_pse, posPow_pse, measPow_pse, varPow_pse;

        float thisPseHistSpacing = thisPseHist->GetBinCenter(2) - thisPseHist->GetBinCenter(1);
        for (int iii=1; iii<thisPseHist->GetSize()-1; iii++)
        {
          double pos_here = thisPseHist->GetBinCenter(iii);
          if ( thisRow.xMin <= pos_here && pos_here <= thisRow.xMaxPse + thisPseHistSpacing*trailingZeros )     
          {
            pos_pse.push_back(pos_here);
            meas_pse.push_back(thisPseHist->GetBinContent(iii));
            var_pse.push_back(thisPseHist->GetBinError(iii));
          }
        }
        float thisPsePowHistSpacing = thisPsePowHist->GetBinCenter(2) - thisPsePowHist->GetBinCenter(1);
        for (int iii=1; iii<thisPsePowHist->GetSize()-1; iii++)
        {
          double pos_here = thisPsePowHist->GetBinCenter(iii);
          if ( thisRow.xMin <= pos_here && pos_here <= thisRow.xMaxPowPse + thisPsePowHistSpacing*trailingZeros )
          {
            posPow_pse.push_back(pos_here);
            measPow_pse.push_back(thisPsePowHist->GetBinContent(iii));
            varPow_pse.push_back(thisPsePowHist->GetBinError(iii));
          }
        }



        // define our FCN for minuit (see fitFunc.h/.cc)
        FitFCN thePseFCN(meas_pse, pos_pse, var_pse);
        PowerFCN psePowerFCN(measPow_pse, posPow_pse, varPow_pse);

        // define fit parameters
        // fit(x_i) = exp(a+b*x_i)
        // powFit(x_i) = c*x_i^d

        MnUserParameters upar_pse;
        MnUserParameters uparPow_pse;
        upar_pse.Add("a", thisRow.fitIntercept, 0.1 );
        upar_pse.Add("b", thisRow.fitSlope, 0.1 );
        uparPow_pse.Add("c", thisRow.fitPowScale ,0.1);
        uparPow_pse.Add("d", thisRow.fitPowOrder ,0.1);

        MnMigrad migrad_pse(thePseFCN,upar_pse);
        MnMigrad migradPow_pse(psePowerFCN, uparPow_pse);
        migrad_pse.Fix("a");
        migrad_pse.Fix("b");
        migradPow_pse.Fix("c");
        migradPow_pse.Fix("d");

        FunctionMinimum mmin_pse = migrad_pse();
        FunctionMinimum mminPow_pse = migradPow_pse();

        //fitStatus = mmin_pse.IsValid();
        if (mmin_pse.IsValid() == 1) { fitStatus = 0; }
        if (mmin_pse.IsValid() == 0) { fitStatus = 1; }  //??

        if (fitStatus > 0) {
          vbprintf("   fit error on pseudo-exp %i bin %i \n", pe, thisPair->first);
        }

        //lprintf("hist pointer = %p  n=%f \n", theHist, theHist->GetSumOfWeights());
        double l0; double edm, errdef, l0Pow;
        l0 = mmin_pse.Fval();
        edm = mmin_pse.Edm();
        errdef = mmin_pse.Up();

        l0Pow = mminPow_pse.Fval();

        //thisLogLikeHist->Fill(l0/thisPseHist->GetSumOfWeights());
        //logLikeHistData.push_back(l0/thisPseHist->GetSumOfWeights());
        logLikeHistData.push_back(l0);
        logLikeHistPowData.push_back(l0Pow);

        //cout << "L0S: " << thisRow.logL << " " << l0 << " " << thisRow.logLPow << " " << l0Pow << endl;

        thisPseHist->Delete();
        thisPseHist = 0;
        delete thisPseHist;
        
        thisPsePowHist->Delete();
        thisPsePowHist = 0;
        delete thisPsePowHist;
        if (pe == pseudoExpCount-1) { vbprintf("]\n"); }  //end of pse-exp loop
      }

      // get the maximum value from the array
      float maxLogLike = *max_element(logLikeHistData.begin(), logLikeHistData.end());
      float maxLogLikePow = *max_element(logLikeHistPowData.begin(), logLikeHistPowData.end());

      char histname[128]; sprintf(histname, "logLikeHist%05i", thisPair->first);
      TH1F* thisLogLikeHist = new TH1F(histname, "pseudo-exp log likelihood", 100, 0, max(maxLogLike, thisRow.logL));

      sprintf(histname, "logLikeHistPow%05i", thisPair->first);
      TH1F* thisLogLikeHistPow = new TH1F(histname, "pseudo-exp log likelihood for Powerlaw Fit", 100, 0, max(maxLogLikePow, thisRow.logLPow));


      vector<double> logLikeHistWeights(logLikeHistData.size(), 1.0);
      thisLogLikeHist->FillN(logLikeHistData.size(), &logLikeHistData[0], &logLikeHistWeights[0]);

      vector<double> logLikeHistPowWeights(logLikeHistPowData.size(), 1.0);
      thisLogLikeHistPow->FillN(logLikeHistPowData.size(), &logLikeHistPowData[0], &logLikeHistPowWeights[0]);

      //lprintf("");
      // calculate p-value
      float pVal = 0;
      int sBin = thisLogLikeHist->FindBin(thisRow.logL);
      if (sBin > -1) {
        for (int b=sBin; b<=thisLogLikeHist->GetNbinsX(); ++b) {
          pVal += thisLogLikeHist->GetBinContent(b);
        }
        //cout << "pval pre norm and -1 : " << pVal << endl;
        pVal /= thisLogLikeHist->GetSumOfWeights();
        //pVal = 1.0-pVal;   // TODO reckon on this
      }
      float pValPow = 0;
      int sBinPow = thisLogLikeHistPow->FindBin(thisRow.logLPow);
      if (sBinPow > -1) {
        for (int b=sBinPow; b<=thisLogLikeHistPow->GetNbinsX(); ++b) {
          pValPow += thisLogLikeHistPow->GetBinContent(b);
        }
        //cout << "pval pre norm and -1 : " << pVal << endl;
        pValPow /= thisLogLikeHistPow->GetSumOfWeights();
      }

      cout << "Bin: " << thisPair->first << " Pval: " << pVal << " PvalPow: " << pValPow << endl; 

      vbprintf(" p-value = %f \n", pVal);
      thisRow.logLikePVal = pVal;
      thisRow.logLikeHist = thisLogLikeHist;
      thisRow.logLikePValPow = pValPow;
      thisRow.logLikeHistPow = thisLogLikeHistPow;

      if (true) {
      //if (thisPair->first == ourBin || !drawIt) {
        char canvName[128]; sprintf(canvName, "logLikelihood_%04i", thisPair->first);
        TCanvas* logLikeCanv = new TCanvas("logLikeCanv", canvName, 500, 500);
        thisLogLikeHist->SetTitle(canvName);
        logLikeCanv->cd(0);
        thisLogLikeHist->Draw();
        gStyle->SetOptStat("e");
        logLikeCanv->Draw();
        //vbprintf("cut line: %f,%f,%f,%f \n", thisRow.logL, gPad->GetUymin(), thisRow.logL, gPad->GetUymax());
        TLine* cutLine = new TLine(thisRow.logL, gPad->GetUymin(), thisRow.logL, gPad->GetUymax());
        cutLine->SetLineColor(kRed);
        cutLine->Draw();
        char filename[1024]; sprintf(filename, "%s/%s.png", outputDir, canvName);
        logLikeCanv->SaveAs(filename);
        
      } 
      //else {
      //  thisLogLikeHist->Delete();
      //}
      thisLogLikeHist = 0;
      delete thisLogLikeHist;

      thisLogLikeHistPow = 0;
      delete thisLogLikeHistPow;

      // if we get a pval of one we want to output some info about it.
      // spicifically 
      //   - add one hist, an example pse data hist (with fit)
      //   - add text file output for table of log likelihoods by bin.
      // -Jacob!
      if ( thisRow.logLikePVal == 1.0 ) 
      {
        // make or get the data for the table.  also open the file and ouput the stuff. Junkal
        char outPseFilename[1024];
        sprintf(outPseFilename, "%s/Pval_1_info_%i.txt", outputDir, thisPair->first);
        FILE* pValStatsOutFile = fopen(outPseFilename, "w");
        ofstream pValStatsOut(outPseFilename, ios_base::app);
        fprintf(pValStatsOutFile, "Bin number: %i Background Fit: f(x)=exp(%f + %f * x )\n \n", thisPair->first, thisRow.fitIntercept, thisRow.fitSlope);
        fprintf(pValStatsOutFile, "Real Data: \n \n");
        fprintf(pValStatsOutFile, "Bin center, Bin content, Bin error, f(x), log-likelyhood \n");
        for (int bin_i = 1; bin_i < thisRow.cutHist->GetSize()-1; bin_i++) 
        {
          float pos_i = thisRow.cutHist->GetBinCenter(bin_i);
          float val_i = thisRow.cutHist->GetBinContent(bin_i);
          float err_i = thisRow.cutHist->GetBinError(bin_i);
          float fval_i = exp(thisRow.fitIntercept + thisRow.fitSlope*pos_i);
          float llval_i;
          if (val_i == 0) 
          {
            llval_i = fval_i - val_i;
          } else {
            llval_i = fval_i - val_i + val_i*log(val_i/fval_i);
          }
          if ( thisRow.xMin <= pos_i && pos_i <= thisRow.xMax )
          {
            fprintf(pValStatsOutFile, "%f, %f, %f, %f, %f \n", pos_i, val_i, err_i, fval_i, llval_i );
          }
        }
        fprintf(pValStatsOutFile, "\nPseudo Data: \n\n");
        fprintf(pValStatsOutFile, "Bin center, Bin content, Bin error, f(x), log-likelyhood\n");
        for (int bin_i = 1; bin_i < thisRow.cutHistPse->GetSize()-1; bin_i++)
        {
          float pos_i = thisRow.cutHistPse->GetBinCenter(bin_i);
          float val_i = thisRow.cutHistPse->GetBinContent(bin_i);
          float err_i = thisRow.cutHistPse->GetBinError(bin_i);
          float fval_i = exp(thisRow.fitIntercept + thisRow.fitSlope*pos_i);
          float llval_i;
          if (val_i == 0)
          {
            llval_i = fval_i - val_i;
          } else {
            llval_i = fval_i - val_i + val_i*log(val_i/fval_i);
          }
          if ( thisRow.xMin <= pos_i && pos_i <= thisRow.xMaxPse )
          { 
            fprintf(pValStatsOutFile, "%f, %f, %f, %f, %f\n", pos_i, val_i, err_i, fval_i, llval_i );
          }
        }
        fclose(pValStatsOutFile);

        // plot the pse hist.

        TH1F*& thisHistPse = thisRow.cutHistPse;
        TCanvas* thisCanvPse =  new TCanvas(thisHistPse->GetName(), thisHistPse->GetName(), 450, 450);
        thisCanvPse->cd(0);
        gPad->SetLogy();
        gStyle->SetOptStat("");
        gStyle->SetOptFit(2);
        thisHistPse->Draw();
        f1->Draw("same");
        thisCanvPse->Update();
        char outFilename[1024]; sprintf(outFilename, "%s/%s.png", outputDir, thisCanvPse->GetName());
        thisCanvPse->SaveAs(outFilename);
        sprintf(outFilename, "%s/%s.root", outputDir, thisCanvPse->GetName());
        thisCanvPse->SaveAs(outFilename);
        thisCanvPse = 0;
      }

      delete f1;
      delete f2;
 
    }
  }
  
  sprintf(histTitle, "p-value: cut slope=%3.1f", cutSlope);
  TH1I* pValHist    = new TH1I("pValHist", histTitle, 20, 0, 1);
  TH2F* pValHist2   = new TH2F("pValHist2", histTitle, 50, 0, 1, 50, 0,   0);
  TH1I* pValHistPow = new TH1I("pValHistPow", histTitle, 20, 0, 1);
  int binsOverLim = 0;
  //float pValScore = 0;
  for (map<int, dataRow_t>::iterator thisPair=cutValHist.begin(); thisPair!=cutValHist.end(); ++thisPair) {
    int thisBin = thisPair->first;
    dataRow_t thisRow = thisPair->second;
    vbprintf("Bin number %4i:  events=%8.3f, logL=%6.5f, pValue = %5.4f \n", thisBin, thisRow.cutHist->GetSumOfWeights(), thisRow.logL, thisRow.logLikePVal);
    if (thisRow.cutHist->GetSumOfWeights() > 30.0) {
      pValHist->Fill(thisRow.logLikePVal);
      pValHistPow->Fill(thisRow.logLikePValPow);
      pValHist2->Fill(thisRow.logLikePVal, thisRow.cutHist->GetSumOfWeights());
      //pValScore += log10(thisRow.cutHist->GetSumOfWeights()) * pow(thisRow.logLikePVal - 0.5 ,2);
      if (thisRow.logLikePVal > 0.5) ++binsOverLim;
    } else {
      vbprintf("   ignoring \n");
    }
  }

  //int fitStatus = pValHist->Fit("pol1", fitOpt, "", 0.05, 1.0);
  //double l0, edm, errdef;
  //int nvpar, nparx;        
  //int fitStatus = pValHist->Fit("pol0", fitOpt, "", 0.05, 1.0);

  //TVirtualFitter::GetFitter()->GetStats(l0, edm, errdef, nvpar, nparx);
  //float pValConst0 = TVirtualFitter::GetFitter()->GetParameter(0);
  //float pValLogL0 = l0;

  //fitStatus = pValHist->Fit("pol1", fitOpt, "", 0.05, 1.0);
  //TVirtualFitter::GetFitter()->GetStats(l0, edm, errdef, nvpar, nparx);
  //float pValSlope1 = TVirtualFitter::GetFitter()->GetParameter(1);
  //float pValConst1 = TVirtualFitter::GetFitter()->GetParameter(0);
  //float pValLogL1 = l0;

  //lprintf("final p-value distribution linear log likelihood fit\n");
  //lprintf("  const=%f, logL=%f \n", pValConst0, pValLogL0);
  
  lprintf("  bins over limit %i \n", binsOverLim);
  char outFilename[1024];
  TCanvas* pValCanv = new TCanvas("pValCanv", "pValCanv", 500, 800);
  pValCanv->Divide(1, 2);
  pValCanv->cd(1);
  pValHist->Draw();
  //gPad->SetLogy();
  //pValHist2->SetMarkerSize(6);
  //pValHist2->SetMarkerStyle(kFullDotMedium);
  //pValHist2->Draw();
  //pValHist2->GetYaxis()->SetRangeUser(10, 10000);
  pValCanv->cd(2);
  //pValHist->Draw();
  pValHistPow->Draw();
  gStyle->SetOptStat("");
  gStyle->SetOptFit(2);
  pValCanv->Update();
  sprintf(outFilename, "%s/pValHist2_slope_%03.1f.root", outputDir, -cutSlope);
  pValCanv->SaveAs(outFilename);
  sprintf(outFilename, "%s/pValHist2_slope_%03.1f.png", outputDir, -cutSlope);
  pValCanv->SaveAs(outFilename);
  //printf("pVal score for slope %f is %f \n", cutSlope, pValScore);
  
  float sumPVal = 0;
  int numBins = pValHist->GetNbinsX()-1;
  for (int k=2; k<=pValHist->GetNbinsX(); ++k) {
    sumPVal += pValHist->GetBinContent(k);
  }
  float lambdaP = sumPVal / numBins;
  lprintf("pValue sum=%f, %i bins avg=%f \n", sumPVal, numBins, lambdaP);

  // calculate the poisson log-likelihood L0 value 
  float logLSum = 0;
  int numRetainedBins = 0;
  //float logLambdaP = log(lambdaP);
  for (int k=2; k<=pValHist->GetNbinsX(); ++k) {
    int thisVal = pValHist->GetBinContent(k);
    numRetainedBins += thisVal;
    //float thisLogL = thisVal * logLambdaP;
    float thisLogL = log(TMath::Poisson(thisVal, lambdaP));
    lprintf("bin %i  value=%i  logL=%f \n", k, thisVal, thisLogL);
    logLSum += thisLogL;
  }
  //logLSum -= numRetainedBins*lambdaP;
  lprintf("number of Healpix bins retained = %i \n", numRetainedBins);
  lprintf("log likelihood total = %f \n", logLSum);
  lprintf(" mean p-value bin count %f \n", lambdaP);

  if (doPValFlatness) {
    // throw a zillion histograms, each with N (number of healpix bins)  uniform random numbers
    int aZillion = 1000;
    char title[128];
    sprintf(title, "PseudoExp Log Lilkelihoods - slope -%2.0f", cutSlope);
    TH1F* grandPoobahHist = new TH1F("slope_pseudoexp", title, 50, -40, -15);
    vector<TH1F*> metaHists(aZillion, 0);
    for (int k=0; k<aZillion; ++k) {
      //lprintf("building metaHist %i \n", k);
      char histName[64];
      sprintf(histName, "metaHist%06i", k);
      metaHists[k] = new TH1F(histName, histName, 20, 0, 1);
      for (int j=0; j<numRetainedBins; ++j) {
        float rand = gRandom->Rndm();
        //lprintf("filling event number %i value=%f \n", j, rand);
        metaHists[k]->Fill(rand, 1.0);
      }
      // calculate the Poisson log-likelihood of this histogram and store in in the grand poobah histogram
      float grandPoobahLogL = 0;
      for (int j=2; j<=metaHists[k]->GetNbinsX(); ++j) {
        float thisVal = metaHists[k]->GetBinContent(j);
        //float thisLogL = thisVal * logLambdaP;
        float thisLogL = log(TMath::Poisson(thisVal, lambdaP));
        //lprintf("bin %i (%f) value=%f  logL=%f \n", j, metaHists[k]->GetBinLowEdge(k), thisVal, thisLogL);
        grandPoobahLogL += thisLogL;
      }
      //grandPoobahLogL -= numRetainedBins*lambdaP;
      //lprintf(" for pseudo-exp %i grand poobah log likelihood is %f \n", k, grandPoobahLogL);
      grandPoobahHist->Fill(grandPoobahLogL);
    }
    TCanvas* gpCanv = new TCanvas("gpCanv", "gpCanv", 400, 400);
    gpCanv->cd(0);
    grandPoobahHist->Draw();
    gStyle->SetOptStat("emrou");
    gpCanv->Update();
    //lprintf("p-val line limits are %f,%f,   %f,%f \n", logLSum, gPad->GetUymin(), logLSum,  gPad->GetUymax());
    TLine* l0Line = new TLine(logLSum, gPad->GetUymin(), logLSum, gPad->GetUymax());
    l0Line->SetLineColor(kRed);
    l0Line->Draw();
    // finally get at the p-value of L0 
       // get the bin containing L0
    int l0Bin =  grandPoobahHist->FindFixBin(logLSum);
    float gpPCount = 0;
    for (int k=l0Bin; k<=grandPoobahHist->GetNbinsX(); ++k) {
      gpPCount += grandPoobahHist->GetBinContent(k);
    }
    float gpPVal = gpPCount / (float)(grandPoobahHist->GetSumOfWeights());
    //lprintf("grand poobah total count is %f, count to right of line is %f \n", grandPoobahHist->GetSumOfWeights(), gpPCount);
    lprintf("grand poobah p-value %f \n", gpPVal);

    sprintf(outFilename, "%s/pValStatsBySlope.txt", outputDir);
    FILE* pValStatsOutFile = fopen(outFilename, "a");
    ofstream pValStatsOut(outFilename, ios_base::app);
    fprintf(pValStatsOutFile, "%f,%i,%f\n", cutSlope, binsOverLim, gpPVal); 
    fclose(pValStatsOutFile);  

  }

  // now optimize the cut intercept for each healpix bin
  
  // get the background expectation and keep it in dataRow.cutHist
  for (map<int, dataRow_t>::iterator thisPair=cutValHist.begin(); thisPair!=cutValHist.end(); ++thisPair) {
    int thisHpBin = thisPair->first;
    dataRow_t thisRow = thisPair->second;
    // to get a background, integrate the 10% differential cut fit from y-int cut to infinity (Jacob!)  
    // butttt we dont have the opt y-int value here. what do we do?
    // real bakground est is calcualted later and saved. this one isnt saved? just printed to screen...
    float binWidth = thisRow.cutHist->GetBinWidth(1);    
    float bgNormFac = binWidth;
    cout << "bin " << thisHpBin << " bgNormFac " << bgNormFac << endl; 
    float thisMaxBin = thisRow.cutHist->GetMaximumBin();
    float thisArg = thisRow.cutHist->GetBinLowEdge(thisMaxBin+2);  //place holder
    float thisBackground = -1.0/thisRow.fitSlope/bgNormFac * exp(thisRow.fitSlope * thisArg + thisRow.fitIntercept);
    //float thisBackground = -1.0/thisRow.fitSlope/bgNormFac * exp(thisRow.fitSlope * thisOptCutVal + thisRow.fitIntercept);
    lprintf("healpix bin %i slope=%f, int=%f, xMin=%f logL=%f \n", thisHpBin, thisRow.fitSlope, thisRow.fitIntercept, thisRow.xMin, thisRow.logL );
    lprintf("  background=%f  \n", thisBackground);
    // set up a histogram with the same x-binning as corrSnrHist
    if (thisHpBin == ourBin || !drawIt) {
      char name[64]; sprintf(name, "survivingBackground_%i", thisHpBin);
      //lprintf(" bins=%i min=%f max=%f \n", corrSnrHist->GetNbinsX(), corrSnrHist->GetXaxis()->GetBinLowEdge(1), 
      //        corrSnrHist->GetXaxis()->GetBinLowEdge(corrSnrHist->GetNbinsX()+1));
      TH1F* thisBgHist = new TH1F(name, name, thisRow.cutHist->GetNbinsX(), thisRow.cutHist->GetXaxis()->GetBinLowEdge(1), 
              thisRow.cutHist->GetXaxis()->GetBinLowEdge(thisRow.cutHist->GetNbinsX()+1));
      lprintf(" surviving background hist has %i bins \n", thisBgHist->GetNbinsX());
      for (int k=0; k<=thisBgHist->GetNbinsX(); ++k) {
        float thisBinLowEdge = thisBgHist->GetXaxis()->GetBinLowEdge(k);
        float thisArg = max(thisBinLowEdge, thisRow.xMin);
        //lprintf("integrating background fit function: arg=%f \n", thisArg);
        float thisBinVal = -1.0/thisRow.fitSlope/bgNormFac * exp(thisRow.fitSlope * thisArg + thisRow.fitIntercept); 
        //lprintf("  k=%i, ctr=%f, arg= %f, val=%f \n", k, thisBinCtr, thisArg, thisBinVal);     
        thisBgHist->SetBinContent(k, thisBinVal);
      }
      sprintf(name, "survivingBackground_%i", thisHpBin);
      TCanvas* sbgCanv = new TCanvas(name, name, 300, 300);
      sbgCanv->cd(0);
      //gPad->SetLogy();
      thisBgHist->Draw();
    }
  }
  
  // make hist of the accepted simulated events by cut value
  //sprintf(inFilename, "%s/analysisOutput_2_sim.root", inputDir); //oindree changing to below
  sprintf(inFilename, "/users/PAS0174/osu8620/anita/SamsAnalysis/resultsNewSim/analysisOutput_2_sim.root");
  lprintf("simulation input filename is %s \n", inFilename);
  TFile* simFile = new TFile(inFilename);
  TTree* simTree0 = (TTree*)simFile->Get("analysisOutput0");
  lprintf("input tree 0 has %lli entries \n", simTree0->GetEntries());

  TTree* simTree1 = 0;
  if (!doRebinning) {
    simTree1 = (TTree*)simFile->Get("analysisOutput1");
    lprintf("input tree 1 has %lli entries \n", simTree1->GetEntries());
  }

  sprintf(dummyFilename, "%s/dummy/dummy2_%c_%f_%f.root", inputDir, whichPol, hpThOffset, hpPhiOffset);  
  TFile* dummyFile2 = new TFile(dummyFilename, "RECREATE");
  if (doRebinning) {
    simTree1 = new TTree("analysisOutput1", "analysis1");
  }

  float eventWeight = 0;

  simTree0->SetBranchAddress("runNumber", &runNumber);
  simTree0->SetBranchAddress("eventNumber", &eventNumber);
  simTree0->SetBranchAddress("pol", &pol);
  simTree0->SetBranchAddress("hPeak", &hPeak);
  simTree0->SetBranchAddress("cPeak", &cPeak);
  simTree0->SetBranchAddress("cohSnr", &cohSnr);
  simTree0->SetBranchAddress("cohSnr1", &cohSnr1);
  simTree0->SetBranchAddress("cohSnr2", &cohSnr2);
  simTree0->SetBranchAddress("linPolFrac", &linPolFrac);
  simTree0->SetBranchAddress("minCircPeakVal", &minCircPeakVal);
  simTree0->SetBranchAddress("cPolDist", &cPolDist);
  simTree0->SetBranchAddress("cPolPeak", &cPolPeak);

  simTree0->SetBranchAddress("eventWeight", &eventWeight);

  //simTree0->SetBranchAddress("ellParms", &ellParms);                                   // error ellipse parameters
  simTree0->SetBranchAddress("gpsRaw", &gpsRaw);                                         // raw payload gps info
  //simTree0->SetBranchAddress("gps", &gps);                                               // payload gps info
  simTree0->SetBranchAddress("theta", &theta);                                           // event theta
  simTree0->SetBranchAddress("phi", &phi);                                               // event phi
  simTree0->SetBranchAddress("ea", &ea);                                                 // event easting
  simTree0->SetBranchAddress("no", &no);                                                 // event northing
  simTree0->SetBranchAddress("lat", &lat);                                               // event latitude
  simTree0->SetBranchAddress("lon", &lon);                                               // event longitude
  simTree0->SetBranchAddress("alt", &alt);                                               // event altitude

  if ( !doRebinning ) {
    simTree1->SetBranchAddress("runNumber", &runNumber);
    simTree1->SetBranchAddress("eventNumber", &eventNumber);
    simTree1->SetBranchAddress("pol", &pol);
    simTree1->SetBranchAddress("hpBin", &hpBin);
    simTree1->SetBranchAddress("hpWeight", &hpWeight);
  } else {
    simTree1->Branch("runNumber", &runNumber, "runNumber/I");
    simTree1->Branch("eventNumber", &eventNumber, "eventNumber/I");
    simTree1->Branch("pol", &pol, "pol/B");
    simTree1->Branch("hpBin", &hpBin, "hpBin/I");
    simTree1->Branch("hpWeight", &hpWeight, "hpWeight/F");
  }


  simTree0->BuildIndex("eventNumber", "pol");
  TTreeIndex* simTree0Index = (TTreeIndex*)simTree0->GetTreeIndex();

  // do rebinning for sim events
  if ( doRebinning )
  {
    vector<float> errEllipseParms1 = {0, 0, 0, 0};
    vector<vector<double> > errHexVertices = vector<vector<double> >(0);
    vector<vector<double> > errHexPixels = vector<vector<double> >(0);

    for (int e = 0; e<simTree0Index->GetN(); ++e)
    {
      int entNum = simTree0Index->GetIndex()[e];
      simTree0->GetEntry(entNum);
      int getRes = simTree0->GetEntryWithIndex(eventNumber, pol);
      useThisEntry = true;
      if (getRes <= 0) {lprintf("database error event %i pol %c - skipping \n", eventNumber, pol); useThisEntry = false;}
      if (circPeakSepThreshold>0 && cPolDist>circPeakSepThreshold) {useThisEntry = false;}
      if (circPeakStrengthThreshold>0 && minCircPeakVal<circPeakStrengthThreshold) {useThisEntry = false;}
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
      //vbprintf("     Error hexagon vertices :"); for (int k=0; k<6; ++k) {vbprintf("       %f, %f \n", errHexVertices[0][k], errHexVertices[1][k$
      errHexPixels = ellipsePixels(errEllipseParms1[2], errEllipseParms1[3], errEllipseParms1[0], errEllipseParms1[1], plSrcPhi);


      map<int, float> hpWeights = map<int, float>();
      for (int k=0; k<errHexPixels[0].size(); ++k) {
        double thisLat, thisLon;
        bedmap->EaNoToLonLat(errHexPixels[0][k], errHexPixels[1][k], thisLon, thisLat);
        //vbprintf(" pixel ea/no = %f,%f \n", errHexPixels[0][k], errHexPixels[1][k]);
        //vbprintf(" pixel lat/lon = %f,%f \n", thisLat, thisLon);
        double thisTheta = (-thisLat + 90.0 + hpThOffset) * M_PI/180.0;
        double thisPhi = (thisLon + hpPhiOffset) * M_PI/180.0;
        while (thisPhi < -M_PI) thisPhi += M_PI;
        //vbprintf("  pixel theta/phi = %f,%f \n", thisTheta, thisPhi);
        pointing point(thisTheta, thisPhi);
        int pixelNo = healpix->ang2pix(point);

        float thisWeight = errHexPixels[2][k] * eventWeight;  // I think eventWeight is just 1.0 for data, It is different for sim events

        map<int, float>::iterator thisEntry = hpWeights.find(pixelNo);
        if (thisEntry == hpWeights.end()) {
          vbprintf("  creating bin list entry for event %i %cPol bin %i  %8.3f \n", eventNumber, pol, pixelNo, thisWeight);
          hpWeights.insert(pair<int, float>(pixelNo, thisWeight));
        } else {
          vbprintf("  updating bin list entry for event %i %cPol bin %i  %8.3f \n",
                   eventNumber, pol, pixelNo, thisWeight);
          hpWeights[pixelNo] += thisWeight;
        }
      }
      for (pair<int, float> thisEntry : hpWeights) {
        //runNumber = runNumber;
        hpBin = thisEntry.first;
        hpWeight = thisEntry.second;
        //eventNumber = eventNumber;
        //pol = pol;
        vbprintf(" writing bin list entry for event %i %cPol bin %i  %8.3f \n",
                eventNumber, pol, hpBin, hpWeight);
        simTree1->Fill();
      }
    }
  }

  simTree1->BuildIndex("eventNumber", "pol");
  TTreeIndex* simTree1Index = (TTreeIndex*)simTree1->GetTreeIndex();

    // ------- loop through the simulation results tree and populate the simulation count histograms ------------------------------------
  prevEvent = 0;
  prevPol = 0;
  for (int e=0; e<simTree1Index->GetN(); ++e) {
    int entNum = simTree1Index->GetIndex()[e];
    //lprintf("entry number %i \n", entNum); 
    simTree1->GetEntry(entNum);
    if (pol != whichPol) {continue;}
    vbprintf(" event number %i %cPol bin %i \n", eventNumber, pol, hpBin);
    //int p = (pol=='H' || pol=='L') ? 0 : 1;
    if (eventNumber != prevEvent || pol != prevPol) {
      int getRes = simTree0->GetEntryWithIndex(eventNumber, pol);
      //if (getRes <= 0) {lprintf("database error event %i pol %c bin %i - skipping \n", eventNumber, pol, hpBin); continue;}
      useThisEntry = true;
      if (getRes <= 0) {lprintf("database error event %i pol %c bin %i - skipping \n", eventNumber, pol, hpBin); useThisEntry = false;}
      //if (circPeakSepThreshold>0 && cPolDist>circPeakSepThreshold) {useThisEntry = false;}
      //if (circPeakStrengthThreshold>0 && minCircPeakVal<circPeakStrengthThreshold) {useThisEntry = false;}
      //lprintf(" event number %i %cPol \n", eventNumber, pol) ; 
      //lprintf(" indexed tree get result is %i \n", getRes);

      //trying something different in sim too -Jacob
     
      //cutVal = -cutSlope*2*cPolPeak[1] + cohSnr2;
      cutVal = -cutSlope*cPeak + cohSnr2;
      //cutVal = -cutSlope*cPeak + cohSnr3;

      if (hpBin == ourBin) {vbprintf(" event number %i %cPol, bin=%i, cPeak=%.4f, hPeak=%.1f, cutval=%.4f \n", eventNumber, pol, hpBin, cPeak, hPeak, cutVal);}
      prevEvent = eventNumber;
      prevPol = pol;
    }

    if (!useThisEntry) {continue;}
    map<int, dataRow_t>::iterator thisPair = cutValHist.find(hpBin);
    if (thisPair != cutValHist.end()) {
      thisPair->second.simCutHist->Fill(cutVal, hpWeight); 
      thisPair->second.corrPeakSnrSimHist->Fill(cPeak, cohSnr2, hpWeight); 
      //thisPair->second.corrPeakSnrSimHist->Fill(cPeak, cohSnr3, hpWeight);

    } 
    //lprintf(" bin number %i, %.3f, %.4f \n", hpBin, hpWeight, cutVal); 
    //break; 
  }  
  
  // draw and save only - no analysis processing in this loop please
  for (pair<int, dataRow_t> thisPair : cutValHist) {
    dataRow_t& thisRow = thisPair.second;
    if (thisRow.status == 0) {
      if (thisPair.first == ourBin || !drawIt) {
        //vbprintf(" bin %i, status=%i, slope=%f, intercept=%f, logL=%f, xMin=%f, xMax=%f, yMin=%f, yMax=%f \n", 
        //        thisPair.first, thisRow.status, thisRow.fitSlope, thisRow.fitIntercept, thisRow.logL,
        //        thisRow.xMin, thisRow.xMax, thisRow.yMin, thisRow.yMax);
        TH1F*& thisHist = thisPair.second.simCutHist;
        //theHist = thisHist;
        //vbprintf("hist pointer = %p  n=%f \n", thisPair.second.cutHist, thisPair.second.cutHist->GetSumOfWeights());
        //vbprintf("hist pointer = %p  n=%f \n", thisHist, thisHist->GetSumOfWeights());
        vbprintf("  drawing diff plot for bin %i, \n", thisPair.first);
        TCanvas* thisCanv =  new TCanvas(thisHist->GetName(), thisHist->GetName(), 450, 450);  
        //canvi.push_back(thisCanv);
        thisCanv->cd(0);
        gPad->SetLogy();
        gStyle->SetOptStat("");
        gStyle->SetOptFit(2);
        thisHist->Draw();
        //singleCutValCanv->Draw();
        //vbprintf("hist pointer = %p  n=%f \n", thisHist, thisHist->GetSumOfWeights());
        char outFilename[1024]; sprintf(outFilename, "%s/%s.png", outputDir, thisCanv->GetName());
        thisCanv->SaveAs(outFilename);
        sprintf(outFilename, "%s/%s.root", outputDir, thisCanv->GetName());
        thisCanv->SaveAs(outFilename);
        thisCanv = 0;
      }
    }
  }

  //////////////////////////////////////////////////////////////////////
  //
  // now calculate and histogram S and S_up, and their ratio S/S_up
  //
  // Note:
  //   - previously this was done using GSL solver
  //   - that is what the below large comment block is.
  //   - cant do that now that we are taking into account the background uncertanties.
  //   
  // New method:
  //   - for each possible value of cutVal
  //   -- fill a histogram using randomized values of fitSlope and fitIntercept
  //   --- random values thrown assuming gaus
  //   ---- the parameters are corrolated.  so we pick them from
  //        a 2d guas using the dartboard method.
  //   -- sum a large number of these hists into a signal likelihood hist
  //   --- this single hist will be smeared by the fit uncertanties due to the summing.
  //   --- L_smear(s) = Sum[L(s)] over randBackgrounds
  //   -- find the s_up that leaves only 10% of the skewed likelihood hist above s_up
  //   - pick the best s_up
  //   -- the best s_up will be the one that.... uhhh... ????
  //
  // -Jacob!
  //
  /////////////////////////////////////////////////////////////////////////

  // loop over different hpBins
  for (map<int, dataRow_t>::iterator thisPair=cutValHist.begin(); thisPair!=cutValHist.end(); ++thisPair) 
  {
    int thisHpBin = thisPair->first;
    //lprintf(" performing optimization calcs for bin %i \n", thisHpBin);
    dataRow_t& thisRow = thisPair->second;
    if (thisRow.status == 0) 
    {
      // our hist we will calc the background est (and S?) from. its simulated 90% data?
      TH1F* thisHist = thisRow.simCutHist;
      float binWidth = thisRow.cutHist->GetBinWidth(1);     //the same no matter what our sim spacing is.
      float bgNormFac = binWidth;

      lprintf("bin %i  fitSlope=%f+/-%f, fitIntercept=%f+/-%f,  \n", thisHpBin, thisRow.fitSlope, thisRow.fitSlopeError, thisRow.fitIntercept, thisRow.fitInterceptError);

      vbprintf("   [       sUp optimization percent completion       ] \n");
    
      // itterate through different possible values of cut vals
      for (int k = 1; k<=thisHist->GetNbinsX(); ++k) {

        if (k == 1) { vbprintf("   ["); fflush(stdout); }
        else if (k == thisHist->GetNbinsX()) { vbprintf("]\n"); }
        else if (k/(cutValSimBins/50) == floor(k/(cutValSimBins/50))) { vbprintf("="); fflush(stdout); }

        // sSignal is the amount of signal above the cut for a given cut val.
        float cutVal = thisHist->GetXaxis()->GetBinLowEdge(k);
        float sSignal = thisHist->Integral(k, thisHist->GetNbinsX());
        //float thisArg = thisRow.cutHist->GetBinLowEdge(k);
        float thisArg = cutVal;

        // define smear vector (had it as a hist, but we arent plotting it, so a vector is fine.)
        float sMin = 0;
        float sMax = 100;  // what would resonable values be?  100 should be more than enough...
        int sNumBins = 250;
        //TH1F* smearedL = new TH1F("name","title",sNumBins,sMin,sMax);
        vector<double> smearedL;

        // we need to find a reasonable value for sMax based on the current situation.
        // to do this, we will find the s location where the normalized L(s) is 
        // ~ equal to 0.01 using the actual background.
        //
        int numTests = 0;
        float testBG = -1.0/thisRow.fitSlope/bgNormFac * exp(thisRow.fitSlope * thisArg + thisRow.fitIntercept);
        sMax = testBG;
        double maxL = exp( testBG*log(testBG)-testBG-lgamma(testBG+1)-log( ROOT::Math::inc_gamma_c(testBG+1.0,testBG) )  );
        double testL = exp( testBG*log(sMax+testBG)-(sMax+testBG)-lgamma(testBG+1)-log( ROOT::Math::inc_gamma_c(testBG+1.0,testBG) )  );
        double percentL = 0.001;
        while ( ( testL/maxL < percentL*0.5 || percentL*1.5 < testL/maxL) && numTests < 1000 )
        {
          if (testL > percentL*maxL)
          {
            sMax *= 2;
          } else if (testL < percentL*maxL) {
            sMax *= 0.75;
          }
          testL = exp( testBG*log(sMax+testBG)-(sMax+testBG)-lgamma(testBG+1)-log(ROOT::Math::inc_gamma_c(testBG+1.0,testBG)) );
          numTests ++;
          //cout << sMax << " " << testL/maxL << endl;
        }
        //cout << "After " << numTests << " tests sMax reached at " << testL/maxL << " percent of maxL and sMax of " << sMax<< endl; 
        //sleep(1);

        double smearedLTotal = 0; 
        // iterate over 1000ish different backgrounds, add all of them to the smeared hist
        //int num_bg_exp = 10000;   // this value is stored at the start of the file
        for (int bg_i = 0; bg_i < num_bg_exp; bg_i++)
        {
          float thisBackground = 0;

          //check flag to find how we are setting the uncertainty from our fit
          if (true)  // we do this either way, we use the lower error no matter what, but the higher error only if fitErrorFlag == 0
          {

            // make our new random fitSlope and fitIntercept calculate the background
            // - we use a 2d gaus and use the dartboard method. 
            // - also accept only negative slopes, (ycalc)
            float sAvg = thisRow.fitErrorSlopeAvg;
            float iAvg = thisRow.fitErrorInterceptAvg;
            float sSigma = thisRow.fitErrorSlopeSigma;
            float iSigma = thisRow.fitErrorInterceptSigma;
            float row = thisRow.fitErrorRow;

            pair<float,float> fitParamsR;
            fitParamsR = getRandParamBackground ( sAvg, iAvg, sSigma, iSigma, row );
  
            float slopeR = fitParamsR.first; 
            float interceptR = fitParamsR.second;

            thisBackground = -1.0/slopeR/bgNormFac * exp(slopeR * thisArg + interceptR);
            thisBackground *= (0.9 / sampleFrac);
          }
          float choiceError = 0;
          if (fitErrorFlag == 1)
          {
            // pull fit Choice Error from a lognormal distribution.
            float expBG = -1.0/thisRow.fitSlope/bgNormFac * exp(thisRow.fitSlope * thisArg + thisRow.fitIntercept);
            expBG *= (0.9 / sampleFrac);
            float powBG = -thisRow.fitPowScale/(1+thisRow.fitPowOrder) * pow(thisArg,(1+thisRow.fitPowOrder));
            powBG *= (0.9 / sampleFrac);

            float avgBG = (expBG+powBG)/2.0;
            float bgVar = abs(avgBG-expBG);

            choiceError = expBG - getRandLogNormal(avgBG, bgVar);
            thisBackground += choiceError;
          }

          float inflow = 0;
          float outflow = 0;

          // We need to add in the background from neighboring bins flowing into this bin, 
          // and from this bin flowing out.  To do this find spillover into this bin and itterate through it.
          //  -each bin will have a % contribution that is something like %*thatBackground.
          //  -to randomize it multiply by a fractional poisson based on the number of events that
          //   contributed to that spillover.
          
          pair<float,float> outInflowR;

          outInflowR = getSpilloverOutIn( spillover, cutValHist, thisHpBin, thisArg);
          inflow = outInflowR.second;
          outflow = outInflowR.first;
          
          thisBackground += inflow;
          thisBackground -= outflow;
          if (thisBackground < 0) { thisBackground = 0; }

          //thisBackground =  -1.0/thisRow.fitSlope/bgNormFac * exp(thisRow.fitSlope * thisArg + thisRow.fitIntercept);
          //thisBackground *= (0.9 / sampleFrac);
          
          double kay = (double)gRandom->Poisson(thisBackground);
          double lnIncGamma = lgamma(kay+1.0)+log( ROOT::Math::inc_gamma_c(kay+1.0,(double)thisBackground) ); 

          // iterate over s, to fill the hist so we get L_smear(s). 
          for ( int s_i = 0; s_i <sNumBins ; s_i++)
          {
            float thisS = sMin+((sMax-sMin)/sNumBins)*(float)s_i;
            // we might need extra persision for the likelihood calc, so cast these guys as doubles.
            double lambda = (double)thisS+(double)thisBackground;
            //kay defined before loop to save time

            //calc log of L first, so stuff doesnt blow up.
            double logThisL_s;
            if ( kay == 0 && lambda == 0 )
            {
              logThisL_s = -lambda;
            } else {
              logThisL_s = kay*log(lambda)-lambda-lnIncGamma;
            }
            double thisL_s = exp(logThisL_s);
            //double thisL_s = pow(lambda,kay)*exp(lambda)/gamma(kay+1);

            //normalize for bin width later.  Because we integrate by summing down below.
            //normalize for number of bg trials later
 
            if(bg_i == 0)
            {
              smearedL.push_back( thisL_s );
            } else {
              smearedL[s_i] += thisL_s;
            }
            smearedLTotal += thisL_s;
          }
        }
      
        // we now have out L_smear(s) and we just need to find s_up.
        // s_up is the position of s where 90% of the data is below s, and 10% is above.
        double smearedLSum = 0;
        float sUp = 0;
        double normFactor = (sMax-sMin)/sNumBins / num_bg_exp;
        for (int s_i = 0; s_i < sNumBins; s_i++)
        {
          smearedLSum += smearedL[s_i] * normFactor;
          if (smearedLSum >= 0.9)
          {
            sUp = sMin+((sMax-sMin)/sNumBins)*(float)s_i;
            break;
          }
          if (s_i == sNumBins - 1)
          {
            //Note, this means sUp is actually
            //at a higher value and we failed to find it!            
            sUp = sMin+((sMax-sMin)/sNumBins)*(float)s_i;
          }

        }

        /*
        //print L_smear(s) to a file, for a sanaty check. 
        char outLFilename[1024];
        sprintf(outLFilename, "%s/L_%i.txt", outputDir, thisPair->first);
        FILE* LsmearOutFile = fopen(outLFilename, "w");
        ofstream LsmearOut(outLFilename, ios_base::app);
        fprintf(LsmearOutFile, "{s,L(s)}\n");
        for (int s_i = 0; s_i < sNumBins; s_i++)
        {
          fprintf(LsmearOutFile, "{%f,%f,%f},\n",sMin+((sMax-sMin)/sNumBins)*(float)s_i,smearedL[s_i],smearedLTotal);
        }
        fclose(LsmearOutFile);
        */

        //for a spicific cutval (10?) save L_smear Hist.
        if ( thisArg == 10 ) 
        {
          char smeartitle[1024]; char smearname[1024];
          sprintf(smeartitle, "Sample L_smear and L for bin %i, %cPol;s;L", thisPair->first, pol);
          sprintf(smearname, "L_smear%c%05i", pol, thisPair->first);
          TH1F* smearedLHist = new TH1F(smearname,smeartitle,sNumBins,sMin,sMax);
         
          // fill hist with data
          for (int iL = 0; iL < sNumBins; iL++)
          {
            smearedLHist->Fill(sMin+((sMax-sMin)/sNumBins)*((float)iL+0.5),smearedL[iL]/(float)num_bg_exp);
          }
          // add line for actual L.
          float expBG = -1.0/thisRow.fitSlope/bgNormFac * exp(thisRow.fitSlope * thisArg + thisRow.fitIntercept);
          expBG *= (0.9 / sampleFrac);
          double gammaNorm = lgamma(expBG+1.0)+log( ROOT::Math::inc_gamma_c(expBG+1.0,expBG) ); 
          TF1 * likeF = new TF1("likeF","exp( [0]*log([0]+x) - [0] - x - [1]  )",sMin,sMax);
          likeF->SetParameter(0,expBG);
          likeF->SetParameter(1,gammaNorm);

          // make actual plot things happen
   
          TCanvas* thisCanv =  new TCanvas(smearedLHist->GetName(), smearedLHist->GetName(), 400, 400);
          //thisCanv->Divide(2, 1);
          thisCanv->cd(1);
          //gPad->SetLogy();
          smearedLHist->Draw("HIST L");

          likeF->SetLineColorAlpha(kRed,1);
          likeF->Draw("SAME L");

          //make vertical line at s_up
          TLine * lineSup = new TLine(sUp,gPad->GetUymin(),sUp,gPad->GetUymax());
          lineSup->SetLineColor(kRed);
          lineSup->Draw("SAME L");

          char smearOutFile[1024];
          sprintf(smearOutFile, "%s/%s.png", outputDir, thisCanv->GetName());
          lprintf("saving canvas %s \n", smearOutFile);
          thisCanv->SaveAs(smearOutFile);
          sprintf(smearOutFile, "%s/%s.root", outputDir, thisCanv->GetName());
          thisCanv->Update();
          thisCanv->SaveAs(smearOutFile);
          thisCanv = 0;

          delete likeF;
          delete smearedLHist;
        }

        //cout << "cutval = " << cutVal << " thisArg = " << thisArg << " bg = " << testBG << " sSignal = " << sSignal << " sup = " << sUp << endl;
        //sleep(1);

        // cutvalOffset just makes sure the data falls into the right hist bin.
        // its not smart to put data in right on a bin edge!
        float cutValOffset = 0.5*thisHist->GetBinWidth(1);

        //now save the stuff
        //lprintf("  S_up  is %f \n", gslRoot);
        float s_sUp = (sUp > 0) ? sSignal / sUp : 0;
        thisRow.sHist->Fill(cutVal+cutValOffset ,sSignal);
        thisRow.sUpHist->Fill(cutVal+cutValOffset ,sUp);
        thisRow.s_sUpHist->Fill(cutVal+cutValOffset ,s_sUp);
        //lprintf("sHist pointer = %p \n", thisRow.sHist);

      }
      for (int k=0; k<=thisRow.sHist->GetNbinsX(); ++k) {
        thisRow.sHist->SetBinError(k, 0.001);
        thisRow.sUpHist->SetBinError(k, 0.001);
        thisRow.s_sUpHist->SetBinError(k, 0.001);
      }


  /*

  // set up the gsl solver
  gsl_root_fsolver * gslSolver = gsl_root_fsolver_alloc(gsl_root_fsolver_bisection);
  gsl_function gslFunc;
  struct f_parms parms;
  gslFunc.function = &likelihood0W;
  gslFunc.params = &parms;
  //lprintf("gsl solver instantiated \n");  
    
  parms.alpha = 0.9; // alpha  (actua1ly 1-alpha, careful now airstream driver)
  //for (pair<int, dataRow_t> thisPair : cutValHist) {
  for (map<int, dataRow_t>::iterator thisPair=cutValHist.begin(); thisPair!=cutValHist.end(); ++thisPair) {
  //   this gives us a S for each iteration
    int thisHpBin = thisPair->first;
    //lprintf(" performing optimization calcs for bin %i \n", thisHpBin);
    dataRow_t& thisRow = thisPair->second;
    if (thisRow.status == 0) {
      TH1F* thisHist = thisRow.simCutHist;
      // iterate through the hist we just made by cut value, accumulating S and calculating S/S_up as we go          
      //for (int k = thisRow.startBin; k<thisHist->GetNbinsX(); ++k) {sSignal += thisHist->GetBinContent(k);}
      //lprintf(" start bin = %i: %f \n", thisRow.startBin, thisHist->GetBinLowEdge(thisRow.startBin));
      float binWidth = thisRow.cutHist->GetBinWidth(1);
      float bgNormFac = binWidth;              
      lprintf("bin %i  fitSlope=%f, fitIntercept=%f,  \n", thisHpBin, thisRow.fitSlope, thisRow.fitIntercept);
      //for (int k = thisRow.startBin; k<=thisHist->GetNbinsX(); ++k) {
      for (int k = 1; k<=thisHist->GetNbinsX(); ++k) {
        float cutVal = thisHist->GetXaxis()->GetBinLowEdge(k);
        float sSignal = thisHist->Integral(k, thisHist->GetNbinsX());
        // we now have our S for this iteration
        //   now we need Sup(b): calculate it by integrating the fitted function to get background
        // which thing do we do?  I do not know!
        float thisArg = thisRow.cutHist->GetBinLowEdge(k);
        //float thisOptCutVal = thisRow.s_sUpHist->GetBinLowEdge(thisMaxBin); 
        float thisBackground = -1.0/thisRow.fitSlope/bgNormFac * exp(thisRow.fitSlope * thisArg + thisRow.fitIntercept);   
        //float thisBackground = -1.0/thisRow.fitSlope/bgNormFac * exp(thisRow.fitSlope * thisOptCutVal + thisRow.fitIntercept);
        thisBackground *= (0.9 / sampleFrac);
        //lprintf("startBin = %i, sSignal = %f  thisArg = %f \n", k, sSignal, thisArg);
        //lprintf(" b=%f \n", thisBackground);
        parms.b = thisBackground;
        gsl_root_fsolver_set(gslSolver, &gslFunc, 0, 100000);
        int solverStatus = GSL_CONTINUE;
        double gslRoot = 0;
        for (int k=0; k<100 && solverStatus==GSL_CONTINUE ; ++k) {
          solverStatus = gsl_root_fsolver_iterate(gslSolver);
          gslRoot = gsl_root_fsolver_root (gslSolver);
          double xLo = gsl_root_fsolver_x_lower(gslSolver);
          double xHi = gsl_root_fsolver_x_upper(gslSolver);
          solverStatus = gsl_root_test_interval(xLo, xHi, 0.001, 0.001);
          //lprintf("  %i iterations performed \n", k+1);
        }  
        float sUp = gslRoot;
        //lprintf("  S_up (gsl root) is %f \n", gslRoot);
        float s_sUp = (sUp > 0) ? sSignal / sUp : 0;
        thisRow.sHist->Fill(cutVal ,sSignal);
        thisRow.sUpHist->Fill(cutVal ,sUp);
        thisRow.s_sUpHist->Fill(cutVal ,s_sUp);
        //lprintf("sHist pointer = %p \n", thisRow.sHist);
      }
      for (int k=0; k<=thisRow.sHist->GetNbinsX(); ++k) {
        thisRow.sHist->SetBinError(k, 0.001);
        thisRow.sUpHist->SetBinError(k, 0.001);
        thisRow.s_sUpHist->SetBinError(k, 0.001);
      }
      
      */

      // get the optimized cut intercept from the maximum of the S/S_up plot
      int thisMaxBin = thisRow.s_sUpHist->GetMaximumBin();

      float thisOptCutVal = thisRow.optCutVal;
      //float thisOptCutVal = thisRow.s_sUpHist->GetBinLowEdge(thisMaxBin);

      lprintf("   optimized cut value %f \n", thisOptCutVal);
      // get the expected background by integrating the exponential fit function from the optimized cut value to infinity
      float thisBackground = -1.0/thisRow.fitSlope/binWidth * exp(thisRow.fitSlope * thisOptCutVal + thisRow.fitIntercept); 
      thisBackground *= (0.9 / sampleFrac);       
      //lprintf("   optimized background estimate %f \n", thisBackground);      //Moved these to after we get background errors! -Jacob
      //lprintf("   fit p-val %f \n", thisRow.logLikePVal);
      //float thisSamplePassing = thisRow.cutHist->GetBinContent(thisMaxBin);
      float thisSamplePassing = thisRow.cutHist->Integral(thisMaxBin, thisHist->GetNbinsX());
      float thisSimPassing = thisRow.simCutHist->Integral(thisMaxBin, thisHist->GetNbinsX());
      float thisSampleInFitRange = thisRow.cutHist->Integral(thisRow.startBin, thisRow.endBin);
      //float thisSampleInFitRangeTotal = thisRow.cutHistTotal->Integral(thisRow.startBin, thisRow.endBin);
      //float thisSampleInFitRangeWeighted = thisRow.cutHistWeighted->Integral(thisRow.startBin, thisRow.endBin);
      thisRow.expBackground = thisBackground;
      //thisRow.optCutVal = thisOptCutVal;
      thisRow.samplePassing = thisSamplePassing;
      thisRow.simPassing = thisSimPassing;
      thisRow.sampleInFitRange = thisSampleInFitRange;

      // Now that we have the correct background calculated we can find the 
      // expected errors in the background and save them.

      float thisBGLowError = 0;
      float thisBGHighError = 0;     

      if (true)
      {
        int numBGErrorTests = 1000;
        vector<float> thisBG;
        for (int bg_i = 0; bg_i < numBGErrorTests; bg_i++)
        {

          float sAvg = thisRow.fitErrorSlopeAvg;
          float iAvg = thisRow.fitErrorInterceptAvg;
          float sSigma = thisRow.fitErrorSlopeSigma;
          float iSigma = thisRow.fitErrorInterceptSigma;
          float row = thisRow.fitErrorRow;

          pair<float,float> fitParamsR;
          fitParamsR = getRandParamBackground ( sAvg, iAvg, sSigma, iSigma, row );
          float slopeR = fitParamsR.first;
          float interceptR = fitParamsR.second;
          float bgR = -1.0/slopeR/bgNormFac * exp(slopeR * thisOptCutVal + interceptR);
          bgR *= (0.9 / sampleFrac);

          //get the other uncertanties in here too!

          thisBG.push_back(bgR);
        }
        thisBGLowError = 0;
        thisBGHighError = 0;
        float CDFHere = 0;
        sort( thisBG.begin(), thisBG.end() );
        for ( int error_i = 0; error_i < thisBG.size(); error_i++)
        {
          CDFHere += 1;
          if (CDFHere >= 0.158*numBGErrorTests && thisBGLowError == 0)
          { 
            thisBGLowError = thisBG[error_i];
          } 
          if (CDFHere >= 0.84*numBGErrorTests && thisBGHighError == 0)
          {
            thisBGHighError = thisBG[error_i];
          }
        }
      }
      if ( fitErrorFlag == 1) // if we use this method, we only use it for the positive error bar!
      {
        float expBG = thisBackground;
        float powBG = -thisRow.fitPowScale/(1+thisRow.fitPowOrder) * pow(thisOptCutVal,(1+thisRow.fitPowOrder));
        powBG *= (0.9 / sampleFrac);
        float bgDiff = abs(expBG-powBG);
        thisBGHighError = expBG + bgDiff;
        //if (thisBGLowError < 0) { thisBGLowError = 0; }
      }


      // also add in the spillover error
      float outflow = 0;
      float inflow = 0;
      for (map<int, map<int,pair<double,int> > >::iterator thisSpill=spillover.begin(); thisSpill!=spillover.end(); ++thisSpill)
      {
        int fromHpBin = thisSpill->first;
        map<int,pair<double,int> > thisMap = thisSpill->second;
        for ( map<int,pair<double,int> >::iterator thatCont=thisMap.begin(); thatCont!=thisMap.end(); ++thatCont)
        {
          int toHpBin = thatCont->first;
          if (toHpBin == thisHpBin)
          {
            if (fromHpBin == thisHpBin)
            {
              float spill = 1.0-thatCont->second.first;
              outflow = spill*thisBackground;
            }
            else
            {
              map<int, dataRow_t>::iterator thatPair = cutValHist.find(fromHpBin);
              // make sure we actually have data and a fit for that bin
              if ( thatPair->second.status == 0 && thatPair != cutValHist.end())
              {
                float spill = thatCont->second.first;
                float thatSlope = thatPair->second.fitSlope;
                float thatIntercept = thatPair->second.fitIntercept;
                float thatBackground = -1.0/thatSlope/binWidth * exp(thatSlope * thisOptCutVal + thatIntercept);
                thatBackground *= (0.9 / sampleFrac);
                inflow += spill*thatBackground;
              }
            }
          }
        }
      }

      thisRow.fitErrorHigh = thisBGHighError - thisBackground;
      thisRow.fitErrorLow = thisBackground - thisBGLowError; 
      thisRow.outflow = outflow;
      thisRow.inflow = inflow;

      thisBGLowError -= outflow;
      thisBGHighError += inflow;

      //if (thisBgLowError < 0) { thisBgLowError = 0; }

      thisRow.expBGLowError = thisBGLowError;
      thisRow.expBGHighError = thisBGHighError;

      lprintf("   optimized background estimate %f -%f/+%f \n", thisBackground,thisBGLowError,thisBGHighError);
      lprintf("   fit p-val %f \n", thisRow.logLikePVal);

      // we also probably want to make a hist of this... for comparing. deal with that later?
      //cout << "background estimate = " << thisBackground << " -" << thisBGLowError << "/+" << thisBGHighError << endl;
      //sleep(10);

    }
    float thisSimEventsBeforeLDCut = thisRow.simCutHist->Integral(0, thisRow.simCutHist->GetNbinsX());
    thisRow.simEventsBeforeLDCut = thisSimEventsBeforeLDCut;
  }  

  // create output value for cutval and background estimates. -Jacob!!
  // # Bin, cutval, background, background error
  char cutvalFilename[1024];
  sprintf(cutvalFilename, "%s/cutvalTable.txt", outputDir);
  FILE*  cutvalFileOut = fopen(cutvalFilename, "w");
  ofstream cutvalFile(cutvalFilename, ios_base::app);
  fprintf(cutvalFileOut, "Bin #, Cutval, Signal, Background, - fit error, + fit error, - spill error, +spill error, bg lowbound, bg highbound \n");  
  for (map<int, dataRow_t>::iterator thisPair=cutValHist.begin(); thisPair!=cutValHist.end(); ++thisPair) {
  //for (pair<int, dataRow_t> thisPair : cutValHist) {
    int thisHpBin = thisPair->first;
    dataRow_t& thisRow = thisPair->second;
    if (thisRow.status == 0) {
      float thisCutval = thisRow.optCutVal;
      float thisSignal = thisRow.simPassing;
      float thisBackground = thisRow.expBackground;
      float thisBGLow = thisRow.expBGLowError;
      float thisBGHigh = thisRow.expBGHighError;
      fprintf(cutvalFileOut, "%5i, %6.4f, %6.4f, %6.4f, %6.4f, %6.4f, %6.4f, %6.4f, %6.4f, %6.4f  \n", thisHpBin, thisCutval, thisSignal, thisBackground, 
                       thisRow.fitErrorLow, thisRow.fitErrorHigh, thisRow.outflow, thisRow.inflow, thisBGLow, thisBGHigh);
    }
  }
  fclose(cutvalFileOut);

  
  //gsl_root_fsolver_free (gslSolver);
  
  lprintf("drawing \n");
  // display the S/SuP and other optimizing plots.   no processing here please
  for (pair<int, dataRow_t> thisPair : cutValHist) {
    dataRow_t& thisRow = thisPair.second;
    if (thisRow.status == 0) {
      if ( true ) {
      //if (thisPair.first == ourBin || !drawIt) {
        TLegend* sUpLegend = new TLegend(0.5, 0.8, 0.9, 0.9);
        sUpLegend->AddEntry(thisRow.sHist, "S");      
        sUpLegend->AddEntry(thisRow.sUpHist, "S_up");      
        //lprintf("sHist pointer = %p \n", thisRow.sHist);
        TCanvas* thisCanv =  new TCanvas(thisRow.sHist->GetName(), thisRow.sHist->GetName(), 800, 800);  
        thisCanv->Divide(2,2);
        int padNum = 0;
        thisCanv->cd(++padNum);
        gPad->SetLogy();
        thisRow.cutHist->Draw();
        //add uncertainty bands onto cutHist.
        //do this by making rand exponents and putting them on the plot
        /*
        TF1 * expoR = new TF1("expoR","exp([0]+[1]*x)",cutValMin,cutValMax);
        for (int i=0;i<1000;i++)
        {
          float slopeR = gRandom->Gaus(thisRow.fitSlope,thisRow.fitSlopeError);
          float interceptR = gRandom->Gaus(thisRow.fitIntercept,thisRow.fitInterceptError);
          int countsR = 0;
          while (slopeR < thisRow.fitSlope-thisRow.fitSlopeError || thisRow.fitSlope+thisRow.fitSlopeError < slopeR)
          {
            //cout << "slope while " << countsR << endl;
            slopeR = gRandom->Gaus(thisRow.fitSlope,thisRow.fitSlopeError);
            countsR++;
          }
          countsR = 0;
          while (interceptR < thisRow.fitIntercept-thisRow.fitInterceptError || thisRow.fitIntercept+thisRow.fitInterceptError < interceptR)
          {
            //cout << "intercept while " << countsR << endl;
            interceptR = gRandom->Gaus(thisRow.fitIntercept,thisRow.fitInterceptError);
            countsR++;
          }
          expoR->SetParameter(0,interceptR);
          expoR->SetParameter(1,slopeR);
          expoR->SetLineColorAlpha(kRed,1);
          expoR->Draw("SAME L");
        }
        //delete expoR;
        */
        char expoChar[1024];
        sprintf(expoChar,"exp(%f+%f + (%f-%f)*x)",thisRow.fitErrorInterceptAvg,thisRow.fitContourMax.first,thisRow.fitErrorSlopeAvg,thisRow.fitContourMax.second);
        cout << expoChar << endl;
        TF1 * expoMaxPlus = new TF1("expoMaxPlus",expoChar,cutValMin,cutValMax);
        sprintf(expoChar,"exp(%f-%f + (%f+%f)*x)",thisRow.fitErrorInterceptAvg,thisRow.fitContourMax.first,thisRow.fitErrorSlopeAvg,thisRow.fitContourMax.second);
        cout << expoChar << endl;
        TF1 * expoMaxMinus = new TF1("expoMaxMinus",expoChar,cutValMin,cutValMax);
        sprintf(expoChar,"exp(%f+%f + (%f+%f)*x)",thisRow.fitErrorInterceptAvg,thisRow.fitContourMin.first,thisRow.fitErrorSlopeAvg,thisRow.fitContourMin.second);
        cout << expoChar << endl;
        TF1 * expoMinPlus = new TF1("expoMinPlus",expoChar,cutValMin,cutValMax);
        sprintf(expoChar,"exp(%f-%f + (%f-%f)*x)",thisRow.fitErrorInterceptAvg,thisRow.fitContourMin.first,thisRow.fitErrorSlopeAvg,thisRow.fitContourMin.second);
        cout << expoChar << endl;
        TF1 * expoMinMinus = new TF1("expoMinMinus",expoChar,cutValMin,cutValMax);
        expoMaxPlus->SetLineColorAlpha(kRed,0.3);
        expoMinPlus->SetLineColorAlpha(kRed,0.3);
        expoMaxMinus->SetLineColorAlpha(kRed,0.3);
        expoMinMinus->SetLineColorAlpha(kRed,0.3);
        expoMaxPlus->Draw("SAME L");
        expoMinPlus->Draw("SAME L");
        expoMaxMinus->Draw("SAME L");
        expoMinMinus->Draw("SAME L");
        //and cutval vert line
        thisCanv->Update();
        TLine * lineOptCut = new TLine(thisRow.optCutVal,pow(10,gPad->GetUymin()),thisRow.optCutVal,pow(10,gPad->GetUymax()));
        lineOptCut->SetLineColor(kRed);
        lineOptCut->Draw("SAME L");
        gStyle->SetOptStat("");
        thisCanv->cd(++padNum);
        thisRow.s_sUpHist->Draw();
        gStyle->SetOptStat("");
        thisCanv->cd(++padNum);
        gPad->SetLogy();
        //thisRow.sHist->Draw();
        thisRow.sHist->GetYaxis()->SetRangeUser(1, 10000);
        thisRow.sHist->Draw();
        thisRow.sUpHist->SetLineColor(kRed);
        thisRow.sUpHist->Draw("same");
        sUpLegend->Draw();
        gStyle->SetOptStat("");
        thisCanv->Draw();
        thisCanv->cd(++padNum);
        //thisRow.logLikeHist->Draw();
        thisRow.corrPeakSnrHist->Draw("COLZ");
        thisRow.corrPeakSnrSimHist->Draw("same");
        thisCanv->Draw();
        //TLine* cutLine = new TLine(thisRow.logL, gPad->GetUymin(), thisRow.logL, gPad->GetUymax());
        float cx0 = gPad->GetUxmin();
        float cy0 = cutSlope * cx0 + thisRow.optCutVal;
        float cy1 = gPad->GetUymin();
        float cx1 = (cy1 - thisRow.optCutVal)/cutSlope;        
        lprintf("bin %i      (%f,%f) (%f,%f) \n", thisPair.first, cx0, cy0, cx1, cy1);
        //TLine* cutLine = new TLine(thisRow.logL, gPad->GetUymin(), thisRow.logL, gPad->GetUymax());
        TLine* cutLine = new TLine(cx0, cy0, cx1, cy1);        
        cutLine->SetLineColor(kRed);
        cutLine->Draw();
        char pValText[32]; sprintf(pValText, "p-val=%5.4f", thisRow.logLikePVal);
        TText* thisText = new TText(0.6, 0.8, pValText);
        thisText->SetNDC(true);
        thisText->SetTextSize(0.04);
        thisText->SetTextFont(42);
        //thisText->Draw();
        //thisCanv->Draw();

        char outFilename[1024]; sprintf(outFilename, "%s/%s.png", outputDir, thisCanv->GetName());
        lprintf("saving canvas %s \n", outFilename);
        thisCanv->SaveAs(outFilename);
        sprintf(outFilename, "%s/%s.root", outputDir, thisCanv->GetName());
        thisCanv->Update();
        thisCanv->SaveAs(outFilename);
        //sprintf(outFilename, "%s/%s.pdf", outputDir, thisCanv->GetName());
        //thisCanv->Update();
        //thisCanv->SaveAs(outFilename);

        thisCanv = 0;

        thisCanv =  new TCanvas(thisRow.simCutHist->GetName(), thisRow.simCutHist->GetName(), 400, 400);  
        //thisCanv->Divide(2, 1);
        thisCanv->cd(1);
        gPad->SetLogy();
        thisRow.simCutHist->Draw();

        sprintf(outFilename, "%s/%s.png", outputDir, thisCanv->GetName());
        lprintf("saving canvas %s \n", outFilename);
        thisCanv->SaveAs(outFilename);
        sprintf(outFilename, "%s/%s.root", outputDir, thisCanv->GetName());
        thisCanv->Update();
        thisCanv->SaveAs(outFilename);
        thisCanv = 0;        
      }
    }
  }

  vector<pair<int, dataRow_t> > vecMap(0);
  for (pair<int, dataRow_t> thisPair : cutValHist) {vecMap.push_back(thisPair);}
  sort(vecMap.begin(), vecMap.end(), sortFunc);

  // Mark the bins we dont want to use

  //We want to mark the last 1% of bins to not be used.
  double totalSimPreLD = 0;
  double totalSimPassing = 0;
  //for (pair<int, dataRow_t> thisPair : vecMap) {
  for (int i = 0; i< vecMap.size(); i++) {
    pair<int, dataRow_t> thisPair = vecMap[i];
    //sleep(1);

    dataRow_t& thisRow = thisPair.second;
    //cout << "bin " << thisPair.first;
    if (thisRow.status != 0){
      //cout << " cut for lack of data, binStatus = " << thisRow.binStatus << endl;
      continue;
    }
    if (thisRow.logLikePVal >= 1.0 ) {
      thisRow.binStatus = 4;
      vecMap[i].second.binStatus = 4;
      //cout << " cut for pval==1, binStatus = " << thisRow.binStatus << endl;
      continue;
    }
    if (thisRow.logLikePVal < 0.05) { 
      thisRow.binStatus = 3; 
      vecMap[i].second.binStatus = 3;
      //cout << " cut for pval<0.05, binStatus = " << thisRow.binStatus << endl;
      continue; 
    }
    //if (thisRow.expBackground > 1.0 ) { 
    //  thisRow.binStatus = 5;
    //  vecMap[i].second.binStatus = 5;
    //  //cout << " cut for bg>1, binStatus = " << thisRow.binStatus << endl;
    //  continue; 
    //}
    //cout << " all good for now, binStatus = " << thisRow.binStatus  << endl;
    //totalSimPreLD += thisRow.simEventsBeforeLDCut;
    //totalSimPassing += thisRow.simPassing;
  }
  /*  
  double sumSimPreLD = 0;
  double sumSimPassing = 0;
  //for (pair<int, dataRow_t> thisPair : vecMap) {
  for (int i = 0; i< vecMap.size(); i++) {
    pair<int, dataRow_t> thisPair = vecMap[i];
    dataRow_t& thisRow = thisPair.second;    
    if (thisRow.status != 0 || thisRow.binStatus != 0){
      continue;
    }
    sumSimPreLD += thisRow.simEventsBeforeLDCut/totalSimPreLD;
    sumSimPassing += thisRow.simPassing/totalSimPassing;
    //if (sumSimPreLD > 0.99)
    if (sumSimPassing > 0.99)
    {
      printf("Bin %i is being cut for low sensitivity, CDF sim = %f \n",thisPair.first,sumSimPassing);
      vecMap[i].second.binStatus = 6;
      thisRow.binStatus = 6;
    }
  }
  */

  //check that binStatus got saved right?
  //for (pair<int, dataRow_t> thisPair : vecMap) {
  //  printf( "(bin,binStatus) = (%i,%i)\n", thisPair.first, thisPair.second.binStatus );
  //}
  //sleep(10);
  
  // now that we have decided the bins to cut... lets move the cutval up and down by a shift to optimize further

  float C_guess = 0;
  float cutValShift = 0;

  float lastBestC = 1;
  float lastBestTP = 1;
  int goingUpCheck = 0;  

  //for (int shift_i = 0; shift_i < 25; shift_i++)
  for (int shift_i = 0; shift_i < 1; shift_i++)
  {
    // shift is applied to all cutvals together.
    // it varies from -1.0 to 4.0? (for now),,, or -0.5 to 2.0
    //float thisShift = (float)(shift_i-5)*0.1;
    float thisShift = 0;
    
    float bestScaleFactor = 0;
    float bestTotalProb = 0;


    for (float scaleFac = startScaleFac; scaleFac > endScaleFac; scaleFac /= scaleFacInc)
    {
      float totalProb = 1.0;

      //loop over bins and calculate new signal, background, and prob.

      //for (pair<int, dataRow_t> thisPair : vecMap) {
      for (int i = 0; i< vecMap.size(); i++) {

        //cout << i << endl;

        pair<int, dataRow_t> thisPair = vecMap[i];

        int thisHpBin = thisPair.first;
        dataRow_t& thisRow = thisPair.second;
        if (thisRow.binStatus != 0) { continue; }

        //cout << thisHpBin << " " << thisRow.binStatus << endl;

        float thisCutVal = thisRow.optCutVal+thisShift;

        float bgNormFac =  thisRow.cutHist->GetBinWidth(1);
        //float newBg = -1.0/thisRow.fitSlope/bgNormFac * exp(thisRow.fitSlope * thisCutVal + thisRow.fitIntercept);
        //newBg *= (0.9 / sampleFrac);
        
        int thisCutValBin = thisRow.simCutHist->GetXaxis()->FindBin(thisCutVal);
        float newSignal = thisRow.simCutHist->Integral(thisCutValBin, thisRow.simCutHist->GetNbinsX());
        //float newSignal2 = thisRow.simCutHist->Integral(thisCutVal*10, thisRow.simCutHist->GetNbinsX());
        //cout << "newSignal: " << newSignal << endl;

        newSignal *= scaleFac;


        // loop over random background pseudo experiments
        // to smear our probablity
        double poissonCdf = 0;
        for (int bg_i = 0; bg_i < num_bg_exp; bg_i++)
        {
          float thisBg = 0;

          // get starting rand background from paramter uncertanties

          float sAvg = thisRow.fitErrorSlopeAvg;
          float iAvg = thisRow.fitErrorInterceptAvg;
          float sSigma = thisRow.fitErrorSlopeSigma;
          float iSigma = thisRow.fitErrorInterceptSigma;
          float row = thisRow.fitErrorRow;
          pair<float,float> fitParamsR;

          fitParamsR = getRandParamBackground ( sAvg, iAvg, sSigma, iSigma, row );

          float slopeR = fitParamsR.first;
          float interceptR = fitParamsR.second;

          thisBg = -1.0/slopeR/bgNormFac * exp(slopeR * thisCutVal + interceptR);
          thisBg *= (0.9 / sampleFrac);

          if (true)
          {
            //cout << "        bin: " << thisHpBin << " paramBG:              "<< thisBg << " slopeR: " << slopeR << " intR: " << interceptR << endl;
          }

          // get fit choice uncertanties:

          float expBG = -1.0/thisRow.fitSlope/bgNormFac * exp(thisRow.fitSlope * thisCutVal + thisRow.fitIntercept);
          expBG *= (0.9 / sampleFrac);
          float powBG = -thisRow.fitPowScale/(1+thisRow.fitPowOrder) * pow(thisCutVal,(1+thisRow.fitPowOrder));
          powBG *= (0.9 / sampleFrac);

          float avgBG = (expBG+powBG)/2.0;
          float bgVar = abs(avgBG-expBG);

          float choiceError = expBG - getRandLogNormal(avgBG, bgVar);
          thisBg += choiceError;

          if (true)
          {
            //cout << "        bin: " << thisHpBin << " paramBG+choice:       " << thisBg << endl;
          }

          // get spillover uncertanties:

          pair<float,float> outInflowR;
          outInflowR = getSpilloverOutIn( spillover, cutValHist, thisHpBin, thisCutVal);

          float inflow = outInflowR.second;
          float outflow = outInflowR.first;
          
          thisBg += inflow;
          thisBg -= outflow;

          if (true)
          {
            //cout << "        bin: " << thisHpBin << " paramBG+choice+spill: " << thisBg << endl;
          }

          if (thisBg <= 0) { thisBg = 0; }

          // now to take into account the statistical uncertainties on the background estaimte
          // we pass 'kay' (thisBg) as poissonD(thisBg)

          double kay = gRandom->PoissonD(thisBg);          

          double thisPoissonCdf = 0.0;
          if (thisBg > pow(10,15)) {
            thisPoissonCdf = 1.0;
          } else {
            thisPoissonCdf = ( ROOT::Math::inc_gamma_c(kay+1.0, newSignal+thisBg)
                                   /ROOT::Math::inc_gamma_c(kay+1.0, thisBg)  );
          }

          //cout << "    s= " << newSignal << " b= " << thisBg << " kay= " << kay << " cdf= " << thisPoissonCdf << endl; 
          //if ( isnan(thisPoissonCdf) ) {
          //  sleep(1);
          //}        

          poissonCdf += thisPoissonCdf/num_bg_exp;
        }
        totalProb *= poissonCdf;
        //cout << " here " << endl;
      }
      //save the totalProb and shift if they are the new best
      cout << "  totalProb=" << totalProb << " scaleFac=" << scaleFac << "; bestProb=" << bestTotalProb << " bestC=" << bestScaleFactor << endl;
      if (totalProb < 0.1 && totalProb > bestTotalProb)
      {
        bestScaleFactor = scaleFac;
        bestTotalProb = totalProb;
      }
    }
    // we have the best scaleFactor for that shift now
    // so now save that to a vector? I guess?
    cout << "For cutValShift of " << thisShift << " best c=" << bestScaleFactor << " with totalProb=" << bestTotalProb << endl;
    //sleep(1);
    if (bestScaleFactor < lastBestC)
    {
      //cout << "  New Best C" << endl;
      cutValShift = thisShift;
      C_guess = bestScaleFactor;
      lastBestC = bestScaleFactor;
      lastBestTP = bestTotalProb;
      //if ( goingUpCheck > 0 ) { goingUpCheck = 0; }
    }
    else if (bestScaleFactor == lastBestC && bestTotalProb < lastBestTP)
    {
      //cout << "  New Best Prob" << endl;
      cutValShift = thisShift;
      C_guess = bestScaleFactor;
      lastBestC = bestScaleFactor;
      lastBestTP = bestTotalProb;
    }
    else if (bestScaleFactor > lastBestC)
    {
      goingUpCheck += 1;
      if (goingUpCheck > 1) 
      {
        //cout << "  C increasing, best C found, ending search." << endl;
        break;  //this saves much time.
      }
    }
  }
  
  cout << "\nOptimized overall shift to optCutVal of " << cutValShift << "\n" << endl;


  // now we have the amount each optimized y-intercept is being shifted by, 
  // so we will more finely sample in c to get an accurate value.
  //
  //  -While we do this we will also update/save our final optCutVal and 
  //   background estimate + uncertainty.
  //  -we can also make a hist of each bin's background estimates.
  //  -also, save the values we need to output for this new cutVal to a new struct

  int num_bg_exp_2 = 20000;  //normally this is 4*num_bg_exp
  map<int,dataRow_final> finalMap;
  map<int,vector<float> > randBackgroundDist;
  map<int,vector<float> > randSingletBgDist;
  for (pair<int, dataRow_t> thisPair : vecMap) {
    int thisHpBin = thisPair.first;
    dataRow_t& thisRow = thisPair.second;
    //if (thisRow.binStatus != 0) { continue; }   
    //if (thisRow.status != 0) { continue; }
    // we will actually calcuate these numbers for bad bins, with fits
    // but when we find the overall shift we wont use those bins.

    dataRow_final tempData;

    //float thisCutVal = thisRow.optCutVal + cutValShift;
    float thisCutVal = 0;
    float thisIntercept = 0;
    float thisSlope = 0;

    // these are bins with naturally calculated cutvals.
    if ( thisHpBin == 3014 ) { thisCutVal =  11.0; thisIntercept = 7.710; thisSlope = -1.119; }   //
    if ( thisHpBin == 3013 ) { thisCutVal =  10.0; thisIntercept = 7.668; thisSlope = -1.185; }   //
    if ( thisHpBin == 3037 ) { thisCutVal =  12.3; thisIntercept = 4.894; thisSlope = -0.7610; }   //
    if ( thisHpBin == 3015 ) { thisCutVal =  10.9; thisIntercept = 9.854; thisSlope = -1.389; }   //
    if ( thisHpBin == 3016 ) { thisCutVal =   9.4; thisIntercept = 9.685; thisSlope = -1.546; }   //
    if ( thisHpBin == 2971 ) { thisCutVal =  11.1; thisIntercept = 8.001; thisSlope = -1.2100; }   //
    if ( thisHpBin == 2970 ) { thisCutVal =  12.0; thisIntercept = 7.759; thisSlope = -1.113; }   //
    if ( thisHpBin == 2936 ) { thisCutVal =   9.3; thisIntercept = 9.715; thisSlope = -1.461; }   //
    if ( thisHpBin == 2989 ) { thisCutVal =   9.9; thisIntercept = 9.955; thisSlope = -1.515; }   //
    if ( thisHpBin == 2937 ) { thisCutVal =  10.2; thisIntercept = 7.238; thisSlope = -1.119; }   //
    if ( thisHpBin == 2990 ) { thisCutVal =  10.0; thisIntercept = 8.117; thisSlope = -1.212; }   //
    if ( thisHpBin == 2988 ) { thisCutVal =  10.1; thisIntercept = 8.826; thisSlope = -1.439; }   //
    if ( thisHpBin == 2991 ) { thisCutVal =   8.7; thisIntercept = 9.920; thisSlope = -1.725; }   //
    if ( thisHpBin == 3029 ) { thisCutVal =  10.9; thisIntercept = 8.694; thisSlope = -1.289; }   //
    if ( thisHpBin == 3003 ) { thisCutVal =  14.4; thisIntercept = 6.901; thisSlope = -0.8453; }   //
    if ( thisHpBin == 3030 ) { thisCutVal =   9.9; thisIntercept = 8.427; thisSlope = -1.301; }   //
    if ( thisHpBin == 2938 ) { thisCutVal =  10.3; thisIntercept = 5.842; thisSlope = -1.099; }   //
    if ( thisHpBin == 3011 ) { thisCutVal =   9.9; thisIntercept = 9.734; thisSlope = -1.561; }   //
    if ( thisHpBin == 3010 ) { thisCutVal =   9.4; thisIntercept = 10.44; thisSlope = -1.727; }   //
    if ( thisHpBin == 3008 ) { thisCutVal =  10.4; thisIntercept = 8.767; thisSlope = -1.331; }   //
    if ( thisHpBin == 2939 ) { thisCutVal =  12.1; thisIntercept = 6.908; thisSlope = -1.107; }   //
    if ( thisHpBin == 2901 ) { thisCutVal =  10.8; thisIntercept = 5.470; thisSlope = -0.9245; }   //
    if ( thisHpBin == 2977 ) { thisCutVal =  11.2; thisIntercept = 7.163; thisSlope = -1.041; }   //
    //if ( thisHpBin == 3001 ) { thisCutVal =  11.4; }   //  sideband
    //if ( thisHpBin == 2978 ) { thisCutVal =  10.0; }   //  sideband
    //if ( thisHpBin == 2900 ) { thisCutVal =  13.7; }   //  sideband


    
    // these are bins with adjusted cutVals
    if ( thisHpBin == 3012 ) { thisCutVal =  12.5; thisIntercept = 7.439; thisSlope = -1.011; }   //12.5
    if ( thisHpBin == 3018 ) { thisCutVal =  10.1; thisIntercept = 7.292; thisSlope = -1.226; }   //10.1
    if ( thisHpBin == 2968 ) { thisCutVal =  16.2; thisIntercept = 3.418; thisSlope = -0.5673; }   //16.2
    if ( thisHpBin == 2997 ) { thisCutVal =   9.3; thisIntercept = 7.455; thisSlope = -1.339; }   //9.3
    if ( thisHpBin == 3019 ) { thisCutVal =  13.5; thisIntercept = 4.971; thisSlope = -0.7723; }   //13.5
    if ( thisHpBin == 3017 ) { thisCutVal =  14.6; thisIntercept = 4.248; thisSlope = -0.6751; }   //14.6
    if ( thisHpBin == 2935 ) { thisCutVal =  14.3; thisIntercept = 4.000; thisSlope = -0.6749; }   //14.3
    if ( thisHpBin == 3004 ) { thisCutVal =  11.7; thisIntercept = 8.691; thisSlope = -1.174; }   //11.7
    if ( thisHpBin == 2998 ) { thisCutVal =  12.2; thisIntercept = 4.977; thisSlope = -0.8474; }   //12.2
    if ( thisHpBin == 2996 ) { thisCutVal =  57.4; thisIntercept = -0.7038; thisSlope = -0.116; }   //57.4
    if ( thisHpBin == 3007 ) { thisCutVal =  11.2; thisIntercept = 8.921; thisSlope = -1.248; }   //11.2
    if ( thisHpBin == 3021 ) { thisCutVal =  26.3; thisIntercept = 5.424; thisSlope = -0.4362; }   //26.3
    if ( thisHpBin == 2979 ) { thisCutVal =  13.1; thisIntercept = 4.557; thisSlope = -0.7685; }   //13.1
    if ( thisHpBin == 2999 ) { thisCutVal =  30.8; thisIntercept = 4.871; thisSlope = -0.3605; }   //30.8
    if ( thisHpBin == 2951 ) { thisCutVal = 168.6; thisIntercept = -0.6156; thisSlope = -0.05241; }   //168.6
    if ( thisHpBin == 2973 ) { thisCutVal =  11.5; thisIntercept = 7.846; thisSlope = -1.131; }   //11.5
    if ( thisHpBin == 2952 ) { thisCutVal = 114.0; thisIntercept = -0.7287; thisSlope = -0.06341; }   //114.0
    
    ///if (thisCutVal == 0) {continue;}

    tempData.optCutVal = thisCutVal;

    if (thisCutVal == 0) {continue;}

    int thisCutValBin = thisRow.simCutHist->GetXaxis()->FindBin(thisCutVal);
    float newSignal = thisRow.simCutHist->Integral(thisCutValBin, thisRow.simCutHist->GetNbinsX());
    tempData.simPassing = newSignal;

    float bgNormFac =  thisRow.cutHist->GetBinWidth(1);

    // loop over random background pseudo experiments
    // to generate a background distribution to use when calcualting 
    // the best C and making bg plots
    float maxBg = 0;
    vector<float> tempBgVec = vector<float>(num_bg_exp_2,0);
    for (int bg_i = 0; bg_i < num_bg_exp_2; bg_i++)
    {
      float thisBg = 0;

      // get starting rand background from paramter uncertanties

      float sAvg = thisRow.fitErrorSlopeAvg;
      float iAvg = thisRow.fitErrorInterceptAvg;
      float sSigma = thisRow.fitErrorSlopeSigma;
      float iSigma = thisRow.fitErrorInterceptSigma;
      float row = thisRow.fitErrorRow;
      pair<float,float> fitParamsR;
      fitParamsR = getRandParamBackground ( sAvg, iAvg, sSigma, iSigma, row );
      float slopeR = fitParamsR.first;
      float interceptR = fitParamsR.second;

      thisBg = -1.0/slopeR/bgNormFac * exp(slopeR * thisCutVal + interceptR);
      thisBg *= (0.9 / sampleFrac);

      // get fit choice uncertanties:

      float expBG = -1.0/thisRow.fitSlope/bgNormFac * exp(thisRow.fitSlope * thisCutVal + thisRow.fitIntercept);
      expBG *= (0.9 / sampleFrac);
      float powBG = -thisRow.fitPowScale/(1+thisRow.fitPowOrder) * pow(thisCutVal,(1+thisRow.fitPowOrder));
      powBG *= (0.9 / sampleFrac);

      float avgBG = (expBG+powBG)/2.0;
      float bgVar = abs(avgBG-expBG);

      float choiceError = expBG - getRandLogNormal(avgBG, bgVar);
      thisBg += choiceError;

      // get spillover uncertanties:

      pair<float,float> outInflowR;
      outInflowR = getSpilloverOutIn( spillover, cutValHist, thisHpBin, thisCutVal);

      float outflow = outInflowR.first;
      float inflow = outInflowR.second;

      thisBg += inflow;
      thisBg -= outflow;

      if (thisBg <= 0) { thisBg = 0; }

      tempBgVec[bg_i] = thisBg;
      if (thisBg > maxBg) { maxBg = thisBg; }
    }
    thisRow.fitSlope = thisSlope;
    thisRow.fitIntercept = thisIntercept;
      
    randBackgroundDist.insert(pair<int,vector<float> >(thisHpBin,tempBgVec));
    randSingletBgDist.insert(pair<int,vector<float> >(thisHpBin,vector<float>(tempBgVec.size(),0)));

    // save peak of hist as central bg est

    // find and save high and low error bounds too
    float thisBGLowError = 0;
    float thisBGHighError = 0;
    float thisBgEst = 0;
    float CDFHere = 0;
    sort( tempBgVec.begin(), tempBgVec.end() );
    for ( int error_i = 0; error_i < tempBgVec.size(); error_i++)
    {
      CDFHere += 1;
      if (CDFHere >= 0.158*num_bg_exp_2 && thisBGLowError == 0)
      {
        thisBGLowError = tempBgVec[error_i];
      }
      if (CDFHere >= 0.84*num_bg_exp_2 && thisBGHighError == 0)
      {
        thisBGHighError = tempBgVec[error_i];
      }
      if (CDFHere >= 0.50*num_bg_exp_2 && thisBgEst == 0)
      { 
        thisBgEst = tempBgVec[error_i];
      }
    }
    tempData.expBGHighError = thisBGHighError;
    tempData.expBGLowError = thisBGLowError;

    // make a histogram of this background distribution here

    char title[128]; sprintf(title, "Background Distribution, bin %i, run %i-%i, %cPol;background;count", thisHpBin, startRun, endRun, pol);
    char name[128]; sprintf(name, "backgroundHist%c%05i", pol, thisHpBin);
    TH1F* bgHist = new TH1F(name, title, 100, 0, thisBGHighError*4);

    for (int bg_i=0; bg_i<num_bg_exp_2; bg_i++)
    {
      bgHist->Fill(tempBgVec[bg_i]);
    }
    //int bgEstBin = bgHist->GetMaximumBin();

    tempData.expBackground = thisBgEst;
    tempData.bgHist = bgHist;


    finalMap.insert(pair<int,dataRow_final>(thisHpBin,tempData));
  }

  // calcualte the singlet bg (distributed between all the bins)
  for (int bg_i = 0; bg_i < num_bg_exp_2; bg_i++)
  {
    double numSingletBg =  gRandom->PoissonD(2.22733);
    
    //find which bin these signlets fall in (H of V)
    while ( numSingletBg > 0.5 ) 
    {
      // in H and V we have 69 bins (40 V-pol, 28 H-pol)
      int binIndex = floor(gRandom->Uniform(0,69));
      if (binIndex < vecMap.size() && finalMap[vecMap[binIndex].first].optCutVal != 0)
      {
        //so... technically I think we are allowing 1 or 2 bins that dont pass to be counted... its probably fine.
        // to fix we want to itterate through randSingletBgDist.size() or something?
        
        randSingletBgDist[vecMap[binIndex].first][bg_i] += 1;

        //cout << binIndex << " " << vecMap[binIndex].first << " " << randSingletBgDist[vecMap[binIndex].first][bg_i] << " " << numSingletBg << endl;
      }
      numSingletBg -= 1.0;
    }
  }

  


  float bestTotalP = 0.0;
  float optC = C_guess;

  for (float scaleFac = C_guess; scaleFac > C_guess/scaleFacInc; scaleFac /= scaleFacIncFine)
  {
    double totalProb = 1.0;
    for (pair<int, dataRow_t> thisPair : vecMap) {
      int thisHpBin = thisPair.first;
      dataRow_t& thisRow = thisPair.second;
      float thisCutVal = finalMap[thisHpBin].optCutVal;
      if (thisCutVal == 0) {continue;}
      //if (thisRow.binStatus != 0 || thisRow.status != 0) { continue; }    //if binStatus != 0, then status shouldnt either, but just to be safe include both.

      //float thisCutVal = thisRow.optCutVal + cutValShift;

      float newSignal = finalMap[thisHpBin].simPassing;

      newSignal *= scaleFac;

      double poissonCdf = 0.0;

      for (int bg_i = 0; bg_i < num_bg_exp_2; bg_i++)
      {
        float thisBg = randBackgroundDist[thisHpBin][bg_i];

        double kay = gRandom->PoissonD(thisBg);

        double thisPoissonCdf = 0.0;
        if (thisBg > pow(10,12)) {
          thisPoissonCdf = 1.0;
        } else {
          thisPoissonCdf = ( ROOT::Math::inc_gamma_c(kay+1.0, newSignal+thisBg)
                            /ROOT::Math::inc_gamma_c(kay+1.0, thisBg)  );
        }

        poissonCdf += thisPoissonCdf/num_bg_exp_2;
      }
      totalProb *= poissonCdf;
    }
    cout << "C = " << scaleFac << " totalProb = " << totalProb << endl;
    if ( totalProb < 0.1 && totalProb > bestTotalP )
    {
      bestTotalP = totalProb;
      optC = scaleFac;
    }
  }


  ////////////////////////////////////////////////////////////////////////////
  //
  //    Do a dumb thing to figure out the limit we can set with 90% confidence.
  //     -This really should be a seperate code
  //     -But I am just putting it here for now for simplicity....
  //
  ////////////////////////////////////////////////////////////////////////////

  cout << endl;

  double realC = C_guess;
  for (float scaleFac = C_guess; scaleFac > C_guess/(scaleFacInc*2); scaleFac /= scaleFacIncFine)
  {
    double totalProb = 1.0;

    float bSum = 0;
    float sSum = 0;
    float NSum = 0;

    for (pair<int, dataRow_t> thisPair : vecMap) {
      int thisHpBin = thisPair.first;
      dataRow_t& thisRow = thisPair.second;

      float thisCutVal = finalMap[thisHpBin].optCutVal;
      if (thisCutVal == 0) {continue;}

      //if (thisRow.binStatus != 0 || thisRow.status != 0) { continue; }    //if binStatus != 0, then status shouldnt either, but just to be safe include both.

      //float thisCutVal = thisRow.optCutVal + cutValShift;

      float newSignal = finalMap[thisHpBin].simPassing;

      newSignal *= scaleFac;

      double N = 0;  // N is the number of events passing in a bin

      // add in passing signlets 
      if (thisHpBin == 3037 || thisHpBin == 2998) { N += 1; }

      /*
      // add in clusters
      if (thisHpBin == 3018 || thisHpBin == 3004) { N += 2; }
      if (thisHpBin == 3016 || thisHpBin == 3029) { N += 4; }
      if (thisHpBin == 3003) { N +=  6; }
      if (thisHpBin == 2979) { N += 10; }
      if (thisHpBin == 3017) { N += 37; }
      */

      /*
      // add in payload blasts and other crap?
      if (thisHpBin == 3019 || thisHpBin == 2998) { N += 1; }
      if (thisHpBin == 3012 || thisHpBin == 3014) { N += 1; }
      if (thisHpBin == 2936 || thisHpBin == 2937) { N += 1; }
      if (thisHpBin == 2990 || thisHpBin == 3011) { N += 1; }
      if (thisHpBin == 3010 || thisHpBin == 2939) { N += 1; }
      if (thisHpBin == 3037 || thisHpBin == 2998) { N += 1; }
      if (thisHpBin == 3015 || thisHpBin == 2988) { N += 2; }
      if (thisHpBin == 3013 || thisHpBin == 3037 || thisHpBin == 3029) { N += 3; }
      if (thisHpBin == 3016) { N += 5; }
      */

      double poissonCdf = 0.0;

      //for (int bg_i = 0; bg_i < 1; bg_i++)
      for (int bg_i = 0; bg_i < num_bg_exp_2; bg_i++)
      {
        
        float bgNormFac =  thisRow.cutHist->GetBinWidth(1);
        float thisthisBg = -0.9/thisRow.fitSlope/bgNormFac/sampleFrac * exp(thisRow.fitSlope * finalMap[thisHpBin].optCutVal + thisRow.fitIntercept);
        //float thisBg = gRandom->PoissonD(thisthisBg) + randSingletBgDist[thisHpBin][bg_i];
        float thisBg = randBackgroundDist[thisHpBin][bg_i] + randSingletBgDist[thisHpBin][bg_i];

        double thisPoissonCdf = 0.0;
        if (thisBg > pow(10,12)) {
          thisPoissonCdf = 1.0;
        } else {
          thisPoissonCdf = ( ROOT::Math::inc_gamma_c(N+1.0, newSignal+thisBg)
                            /ROOT::Math::inc_gamma_c(N+1.0, thisBg)  );
        }
        if (isnan(thisPoissonCdf)) { thisPoissonCdf = 1.0; }
        //cout << "  " << thisHpBin << " " << N << " " << thisBg << " " << newSignal << " " << thisPoissonCdf << endl;

        poissonCdf += thisPoissonCdf/num_bg_exp_2;

        bSum += thisthisBg;
        NSum += N;
        sSum += newSignal;

      }
      //cout << " " << thisHpBin << " " << N << " " << newSignal << " " << poissonCdf << " " << totalProb << endl;
      totalProb *= poissonCdf;
    }

    //cout << " N = " << NSum << endl;
    //cout << " s = " << sSum << endl;
    //cout << " b = " << bSum << endl;

    //totalProb = ( ROOT::Math::inc_gamma_c(NSum+1.0, sSum+bSum)
    //             /ROOT::Math::inc_gamma_c(NSum+1.0, bSum)  );

    cout << "C = " << scaleFac << " totalProb = " << totalProb << endl;
    if ( totalProb < 0.1 && totalProb > bestTotalP )
    {
      bestTotalP = totalProb;
      realC = scaleFac;
    }
  }
  cout << endl << endl << "THE REAL OPTIMIZED C IS: " << realC << " with a prob of " << bestTotalP << endl << endl;


  // Do a check to see what the probablity of seeing only 2 passing events is? (and print the result)

  vector<double> totalPassingVector = vector<double>(num_bg_exp_2,0);
  int lessPassing = 0;

  for (int bg_i = 0; bg_i < num_bg_exp_2; bg_i++)
  {

    double totalPassing = 0;

    for (pair<int, dataRow_t> thisPair : vecMap) 
    {
      int thisHpBin = thisPair.first;
      dataRow_t& thisRow = thisPair.second;

      float thisCutVal = finalMap[thisHpBin].optCutVal;
      if (thisCutVal == 0) {continue;}

      //if (thisRow.binStatus != 0 || thisRow.status != 0 || finalMap[thisHpBin].optCutVal==0 ) { continue; }      
     
      float thisBg = randBackgroundDist[thisHpBin][bg_i] + randSingletBgDist[thisHpBin][bg_i];
      double randPassing = gRandom->PoissonD(thisBg);
      totalPassing += randPassing;
    }
    //cout << bg_i << " " << totalPassing << endl;
    totalPassingVector[bg_i] = totalPassing;
    if (totalPassing < 2.5) { lessPassing += 1; }
  }

  float fracLessPassing = (float)lessPassing / (float)num_bg_exp_2;

  cout << endl << endl << "Fraction of pseudo exp with less than 2 events passing: " << fracLessPassing << endl;

  //totalPassingVector = vector<double>(num_bg_exp_2,0);
  lessPassing = 0;
    
  for (int bg_i = 0; bg_i < num_bg_exp_2; bg_i++)
  {

    double totalPassing = 0;

    for (pair<int, dataRow_t> thisPair : vecMap)
    {
      int thisHpBin = thisPair.first;
      dataRow_t& thisRow = thisPair.second;

      float thisCutVal = finalMap[thisHpBin].optCutVal;
      if (thisCutVal == 0) {continue;}

      //if (thisRow.binStatus != 0 || thisRow.status != 0 || finalMap[thisHpBin].optCutVal == 0) { continue; }

      float bgNormFac =  thisRow.cutHist->GetBinWidth(1);
      float thisthisBg =  -0.9/thisRow.fitSlope/bgNormFac/sampleFrac * exp(thisRow.fitSlope * finalMap[thisHpBin].optCutVal + thisRow.fitIntercept);
      float thisBg = thisthisBg + randSingletBgDist[thisHpBin][bg_i];

      double randPassing = gRandom->PoissonD(thisBg);
      totalPassing += randPassing;
      //cout << " " << thisHpBin << " " << thisBg << " " << randPassing << endl;
      
      //if (bg_i == 0) { cout << "   " << thisHpBin << " " << thisthisBg << " " << finalMap[thisHpBin].simPassing << endl; }

    }

    
    //cout << bg_i << " " << totalPassing << " " << totalPassingVector[bg_i] << endl;

    totalPassingVector[bg_i] = totalPassing;
    if (totalPassing < 2.5) { lessPassing += 1; }
  }

  fracLessPassing = (float)lessPassing / (float)num_bg_exp_2;

  cout << "Fraction of pseudo exp with less than 2 events passing (without systematics): " << fracLessPassing << endl << endl;





  // Do the cut on bins with background > 1 and low sensitivity
  for (int i = 0; i< vecMap.size(); i++) {
    pair<int, dataRow_t> thisPair = vecMap[i];
    dataRow_t& thisRow = thisPair.second;
    //cout << "bin " << thisPair.first;
    //oindree changed below line
    //if (thisRow.status != 0 || thisRow.binStatus != 0){
    if (thisRow.binStatus == 1 || thisRow.binStatus == 2 || thisRow.binStatus == 3 || thisRow.binStatus == 4){
      //cout << " cut for lack of data, binStatus = " << thisRow.binStatus << endl;
      continue;
    }
    //OINDREE INTERVENING
    //if (thisRow.expBackground > 1.0 ) {

    if (finalMap[thisPair.first].expBackground > 1.0 ) {
      thisRow.binStatus = 5;
      vecMap[i].second.binStatus = 5;
      //oindree trying to save these bins, by tuning cut for the bin to be higher to make background 0.5 
      //finalMap[thisPair.first].expBackground = 0.1; 
      //finalMap[thisPair.first].optCutVal = ( ( log(-(finalMap[thisPair.first].expBackground) * sampleFrac * thisRow.fitSlope * (thisRow.cutHist->GetBinWidth(1)) / 0.9) ) - thisRow.fitIntercept ) / thisRow.fitSlope; 
      //cout << " cut for bg>1, binStatus = " << thisRow.binStatus << endl;
      continue; // oindree commented
    }
    else {continue;}
    //cout << " all good for now, binStatus = " << thisRow.binStatus  << endl;
    totalSimPreLD += thisRow.simEventsBeforeLDCut;
    totalSimPassing += thisRow.simPassing;
  } 
  double sumSimPreLD = 0;
  double sumSimPassing = 0;
  //for (pair<int, dataRow_t> thisPair : vecMap) {
  for (int i = 0; i< vecMap.size(); i++) {
    pair<int, dataRow_t> thisPair = vecMap[i];
    dataRow_t& thisRow = thisPair.second;
    if (thisRow.status != 0 || thisRow.binStatus != 0){
      continue;
    }
    sumSimPreLD += thisRow.simEventsBeforeLDCut/totalSimPreLD;
    sumSimPassing += thisRow.simPassing/totalSimPassing;
    //if (sumSimPreLD > 0.99)
    if (sumSimPassing > 0.99)
    {
      printf("Bin %i is being cut for low sensitivity, CDF sim = %f \n",thisPair.first,sumSimPassing);
      vecMap[i].second.binStatus = 6;
      thisRow.binStatus = 6;
    }
  }


  TFile* simCountFile = new TFile("/users/PAS0174/osu8620/anita/SamsAnalysis/analysisSoftware/results/plots/simCountFile_801_808.root");
  TTree* simCountTree = (TTree*) simCountFile->Get("simCountTree");
  int iBinNum;
  float iEventCount;
  simCountTree->SetBranchAddress("binNum", &iBinNum);
  simCountTree->SetBranchAddress("eventCount", &iEventCount);
  simCountTree->BuildIndex("binNum");




  // Make some last minute plots
  // -First super cool bin Status Map
  // -Then hist of rand backgrouds

  TH2F* binStatusHist = 0;
  char binStatusTitle[1024];
  char binStatusName[1024];
  sprintf(binStatusTitle, "HpBin Status Map: run %i-%i, %cPol;ea(km);no(km) ",
          startRun, endRun, whichPol);
  sprintf(binStatusName, "BinStatusHist%1i", whichPol);
  binStatusHist = new TH2F(binStatusName, binStatusTitle, 600, -3000, 3000, 600, -3000, 3000);


  BedmapReader* bedmap = BedmapReader::Instance(false);


  // make up the plot of antarctica to use
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




  for (int gBin = 0; gBin < binStatusHist->GetNcells()-1; ++gBin) {
    int xBin, yBin, zBin;
    binStatusHist->GetBinXYZ(gBin, xBin, yBin, zBin);
    float eag = binStatusHist->GetXaxis()->GetBinCenter(xBin) * 1000;
    float nog = binStatusHist->GetYaxis()->GetBinCenter(yBin) * 1000;
    double thisLat, thisLon;
    bedmap->EaNoToLonLat(eag, nog, thisLon, thisLat);
    double thisTheta = (-thisLat+90.0+hpThOffset) * M_PI/180.0;
    double thisPhi = (thisLon+hpPhiOffset) * M_PI/180.0;
    //while (thisPhi < -M_PI) thisPhi += 2*M_PI;
    pointing point(thisTheta, thisPhi);
    int binNo = healpix->ang2pix(point);

    int gotCount = simCountTree->GetEntryWithIndex(binNo);
    float simCount = 0;
    if (gotCount > 0) {
      simCount = iEventCount;
    }

    for (int i=0; i<vecMap.size(); i++)
    {
      if( vecMap[i].first !=  binNo ) { continue; }

      if (simCount < 0.01 && vecMap[i].second.binStatus == 1) { continue; }
      if ( vecMap[i].second.binStatus == 0 || vecMap[i].second.binStatus == 5 || vecMap[i].second.binStatus == 6)
      //{ binStatusHist->Fill(eag/1000, nog/1000, vecMap[i].second.binStatus+7); }
      { binStatusHist->Fill(eag/1000, nog/1000, 10); }
      else
      { binStatusHist->Fill(eag/1000, nog/1000, vecMap[i].second.binStatus); }
    }
  }

  TCanvas* hpCanv = 0;
  hpCanv = new TCanvas("hpBinStatusCanv", "binStatus", 600, 600);  
  //gStyle->SetPalette(55);  //rainbow
  gStyle->SetPalette(91);    //pastel
  //TColor::InvertPalette();

  // make a legend
  TLegend* binStatusL = new TLegend(0.6, 0.75, 0.9, 0.9);

  TH1F * dummyBox2 = new TH1F("","",1,0,1);
  dummyBox2->SetFillColor(gStyle->GetColorPalette(36));
  binStatusL->AddEntry( dummyBox2,"< 5 bins with data || < 5 events in fit", "F");              // binStatus = 1

  TH1F * dummyBox3 = new TH1F("","",1,0,1);
  dummyBox3->SetFillColor(gStyle->GetColorPalette(73));
  binStatusL->AddEntry( dummyBox3,"Expo Fit Failed", "F");                    // binStatus = 2

  TH1F * dummyBox4 = new TH1F("","",1,0,1);
  dummyBox4->SetFillColor(gStyle->GetColorPalette(109));
  binStatusL->AddEntry( dummyBox4,"Pval < 0.05", "F");                        // binStatus = 3

  TH1F * dummyBox5 = new TH1F("","",1,0,1);
  dummyBox5->SetFillColor(gStyle->GetColorPalette(146));
  binStatusL->AddEntry( dummyBox5,"Pval > 0.999", "F");                       // binStatus = 4

  //TH1F * dummyBox6 = new TH1F("","",1,0,1);
  //dummyBox6->SetFillColor(gStyle->GetColorPalette(182));
  //binStatusL->AddEntry( dummyBox6,"Background Est > 1.0", "F");               // binStatus = 5
  //oindree doing this
  //dummyBox6->SetFillColor(gStyle->GetColorPalette(254));
  //binStatusL->AddEntry( dummyBox6,"Bin kept after cut stricter", "F");               // binStatus = 5

  //TH1F * dummyBox7 = new TH1F("","",1,0,1);
  //dummyBox7->SetFillColor(gStyle->GetColorPalette(218));
  //binStatusL->AddEntry( dummyBox7,"Low Sensitivity", "F");                    // binStatus = 6

  TH1F * dummyBox1 = new TH1F("","",1,0,1);
  dummyBox1->SetFillColor(gStyle->GetColorPalette(254));
  binStatusL->AddEntry( dummyBox1,"Bin Kept", "F");                           // binStatus = 0

  hpCanv->cd(1);
  //gStyle->SetOptStat(0);
  coastLineGr->SetMarkerColor(14);
  coastLineGr->SetFillColor(18);
  coastLineGr->SetLineWidth(1);
  coastLineGr->SetFillStyle(1001);

  binStatusHist->SetStats(kFALSE);
  binStatusHist->Draw("COL");

  coastLineGr->Draw("C");
  double eag, nog;  // for bedmap conversion
  for (pair<int, dataRow_t> thisEntry : cutValHist) {

    int gotCount = simCountTree->GetEntryWithIndex(thisEntry.first);
    float simCount = 0;
    if (gotCount > 0) {
      simCount = iEventCount;
    }
    if (simCount < 0.01 && thisEntry.second.binStatus!=0) { continue; }

    pointing thisPoint = healpix->pix2ang(thisEntry.first);
    double thisLat = 90.0 - (thisPoint.theta * 180.0/M_PI) + hpThOffset ;
    double thisLon = thisPoint.phi * 180.0/M_PI - hpPhiOffset;
    bedmap->LonLattoEaNo(thisLon, thisLat,  eag, nog);
    eag /= 1000;
    nog /= 1000;
    vbprintf("hp bin %i lat,lon is %f,%f  ea,no is %f,%f \n", thisEntry.first, thisLat, thisLon, eag, nog);
    char hpBinText[8]; sprintf(hpBinText, "%i", thisEntry.first);
    vbprintf("drawing hp bin number at %f, %f, %f, %f \n", eag-125, nog, eag+125, nog+100);
    bool posQuad = (eag*nog > 0);
    double textEa = posQuad ? eag : eag-75;
    double textNo = posQuad ? nog-125 : nog+100;
    TText* thisText = new TText(textEa, textNo, hpBinText);
    thisText->SetTextSize(0.025);
    thisText->SetTextAngle(posQuad ? 55 : -55);
    thisText->SetTextFont(42);
    thisText->Draw();
  }
  binStatusL->Draw("S");

  sprintf(filename, "%s/hpBinStatus_%03i_%03i_%c.png", outputDir, startRun, endRun, whichPol);
  hpCanv->SaveAs(filename);
  sprintf(filename, "%s/hpBinStatus_%03i_%03i_%c.root", outputDir, startRun, endRun, whichPol);
  hpCanv->SaveAs(filename);

  // now the bg hists.
  for (pair<int, dataRow_final> thisEntry : finalMap) {
    int hpBin = thisEntry.first;
    TH1F * thisHist = thisEntry.second.bgHist;
    if (thisEntry.second.optCutVal != 0) {
    
      TCanvas* thisCanv =  new TCanvas(thisHist->GetName(), thisHist->GetName(), 800, 800);
      thisHist->Draw();
            
      sprintf(filename, "%s/backgroundHist%c0%i.png", outputDir, whichPol, hpBin);
      thisCanv->SaveAs(filename);
      sprintf(filename, "%s/backgroundHist%c0%i.root", outputDir, whichPol, hpBin);
      thisCanv->SaveAs(filename);
    }
  }

  //sleep(10);
  char summaryString[65536] = "results summary: ";
  sprintf(summaryString+strlen(summaryString), "circPeakSepThreshold =%2.1f, circPeakStrengthThreshold =%4.3f cut slope =%f \n", circPeakSepThreshold, circPeakStrengthThreshold, cutSlope);
  sprintf(summaryString+strlen(summaryString),"                             sim-events     sim-events   sim-events     \n");
  sprintf(summaryString+strlen(summaryString),"                           pre-rotated-cut    passing    before-cuts    \n");
  //for (float scaleFac = 1.0; scaleFac > 0.001; scaleFac /= pow(10, 0.2)) {
  //bool firstShortCycle = true;
  sprintf(filename, "%s/optShortSummary.txt", inputDir); 
  FILE* shortFile = fopen(filename, "a");
  sprintf(filename, "%s/optLongSummary.txt", inputDir); 
  FILE* longFile = fopen(filename, "a");
  bool firstShortCycle = true;
  bool firstLongCycle = true;
  float prevProb = 0;
  for (int i=0; i<2; i++)
  //for (float scaleFac = startScaleFac; scaleFac > endScaleFac; scaleFac /= scaleFacInc) 
  {
    float scaleFac = optC;

    lprintf("\nresults for flux scale factor %f \n", scaleFac);
    //TFile* simCountFile = new TFile("results/plots/simCountFile_801_808.root");  //just moved to earlier
    //TTree* simCountTree = (TTree*) simCountFile->Get("simCountTree");
    //int iBinNum;
    //float iEventCount;
    //simCountTree->SetBranchAddress("binNum", &iBinNum);
    //simCountTree->SetBranchAddress("eventCount", &iEventCount);
    //simCountTree->BuildIndex("binNum");
    int binsFitted = 0;
    int binsAccepted = 0;
    float totalProb = 1.0;
    float simEventsBeforeLDCutSum = 0;
    float simPassingTotal = 0;
    float simCountTotal = 0;
    float simEventsBeforeLDCutSumAcc = 0;
    float simPassingTotalAcc = 0;
    float simCountTotalAcc = 0;
    char spreadSheetStr[65536];
    if ( i == 0 ) {
      sprintf (spreadSheetStr, "Optimization results: slope=%f  simulation scale factor=%f \n", cutSlope, scaleFac);
      sprintf(spreadSheetStr+strlen(spreadSheetStr),"  bin    status    bin   total-events   events-in    sim-events       at-least    optimized       bin-fit  p-value  sim-events     expected      sim-events   Poisson   \n");
      sprintf(spreadSheetStr+strlen(spreadSheetStr)," number   code    status     sample      fit-range  pre-rotated-cut   5-fit-bins  cut-intercept    p-value   >0.05    passing      background     before-cuts    CDF     \n");
    } else {
      sprintf (spreadSheetStr, "Optimization results: slope=%f  simulation scale factor=%f \n", cutSlope, scaleFac);
      sprintf(spreadSheetStr+strlen(spreadSheetStr),"  bin    status    bin   total-events   events-in    sim-events       at-least    optimized       bin-fit  p-value  sim-events     expected      expBG High      expBG Low      sim-events   Poisson   \n");
      sprintf(spreadSheetStr+strlen(spreadSheetStr)," number   code    status     sample      fit-range  pre-rotated-cut   5-fit-bins  cut-intercept    p-value   >0.05    passing      background         Error        Error      before-cuts    CDF     \n");
    }
    //float simCount = 0;
    for (pair<int, dataRow_t> thisPair : vecMap) {
      dataRow_t& thisRow = thisPair.second;
      int gotCount = simCountTree->GetEntryWithIndex(thisPair.first);
      float simCount = 0;
      if (gotCount > 0) {
        simCount = iEventCount*scaleFac;
      }
      float simPassing = 0;
      if (i == 0) {
        simPassing = thisRow.simPassing;
      } else {
        simPassing = finalMap[thisPair.first].simPassing;
      }
      simPassing *= scaleFac;
      float simEventsBeforeLDCut = thisRow.simEventsBeforeLDCut; 
      simEventsBeforeLDCut *= scaleFac;
      int enoughFitBins = (thisRow.status == 0) ? 1 : 0;
      //printf("-------------pValue = %f \n", thisRow.logLikePVal);
      int pValueOkay = (thisRow.logLikePVal > minPVal && thisRow.logLikePVal < 1.0) ? 1 : 0;
      //printf("-------------bin %i   pValue = %f  minPVal = %f    pValueOkay =%i \n", thisPair.first, thisRow.logLikePVal, minPVal, pValueOkay);
      //char inChar;
      //std::cin >> inChar;
      if (enoughFitBins && pValueOkay && thisRow.optCutVal>0) {
        simEventsBeforeLDCutSumAcc += simEventsBeforeLDCut;
        simPassingTotalAcc += simPassing;
        simCountTotalAcc += simCount;
      }
      simEventsBeforeLDCutSum += simEventsBeforeLDCut;
      simPassingTotal += simPassing;
      simCountTotal += simCount;

      // poissonCdf = 1 - integral ( poisson(s+b,b) ) from 0 to s.

      //double poissonCdf = 1.0 - TMath::Gamma(thisRow.expBackground+1.0, thisRow.simPassing*scaleFac);
      // this is wrong, 
      //   -for one, why is simPassing being muliplied by scaleFac again.
      //   -for two, integral of poisson for us =/= the Gamma function, 
      //    because the s+b changes the limits of integration.

      double poissonCdf = 0;
      if (i == 0) {
        poissonCdf = (  ROOT::Math::inc_gamma_c(thisRow.expBackground+1.0, simPassing+thisRow.expBackground)
                       /ROOT::Math::inc_gamma_c(thisRow.expBackground+1.0, thisRow.expBackground)  );
      } else {
        poissonCdf = (  ROOT::Math::inc_gamma_c(finalMap[thisPair.first].expBackground+1.0, simPassing+finalMap[thisPair.first].expBackground)
                       /ROOT::Math::inc_gamma_c(finalMap[thisPair.first].expBackground+1.0, finalMap[thisPair.first].expBackground)  );
      }        
      
      // this has the fixed normalization and corrected simPassing (+bg)

      if (i == 0)
      {
        sprintf(spreadSheetStr+strlen(spreadSheetStr)," %04i     %2i      %2i  %14f  %12f     %12f       %i      %12f  %12f     %i %12f %15f %12f %12f \n", 
                thisPair.first, thisRow.status, thisRow.binStatus, thisRow.numEvents, thisRow.sampleInFitRange, simEventsBeforeLDCut, enoughFitBins, 
                thisRow.optCutVal, thisRow.logLikePVal,  pValueOkay, simPassing, thisRow.expBackground, simCount, poissonCdf);
      } else {
        sprintf(spreadSheetStr+strlen(spreadSheetStr)," %04i     %2i      %2i  %14f  %12f     %12f       %i      %12f  %12f     %i %12f %15f %15f %15f %12f %12f \n",
                thisPair.first, thisRow.status, thisRow.binStatus, thisRow.numEvents, thisRow.sampleInFitRange, simEventsBeforeLDCut, enoughFitBins,
                finalMap[thisPair.first].optCutVal, thisRow.logLikePVal,  pValueOkay, simPassing, finalMap[thisPair.first].expBackground, finalMap[thisPair.first].expBGHighError, finalMap[thisPair.first].expBGLowError, simCount, poissonCdf);
      }

      if (thisRow.status > -1) ++binsFitted;
      //if (thisRow.logLikePVal < 0.05) { thisRow.binStatus = 3; }
      //else if (thisRow.logLikePVal >= 1.0) { thisRow.binStatus = 4; }
      //else if (thisRow.expBackground > 1.0) { thisRow.binStatus = 5; }
      //if ( thisRow.binStatus == 0 ) {     // if we make it to here, we're keeping the bin, at least for now
      if ( thisRow.binStatus == 0 || thisRow.binStatus == 5) {     // oindree doing this to accept bins with high background after tuning cut to make background 0.5 
        totalProb *= poissonCdf;
        ++binsAccepted;
      }
    }
    //char totalString[16384] = "";
    sprintf(spreadSheetStr+strlen(spreadSheetStr),"\n");
    sprintf(spreadSheetStr+strlen(spreadSheetStr),"  TOTAL                        %12f                                                    %12f                 %12f    \n", 
            simEventsBeforeLDCutSum, simPassingTotal, simCountTotal);
    sprintf(spreadSheetStr+strlen(spreadSheetStr),"  ACCEPTED                     %12f                                                    %12f                 %12f    \n", 
            simEventsBeforeLDCutSumAcc, simPassingTotalAcc, simCountTotalAcc);
    sprintf(spreadSheetStr+strlen(spreadSheetStr),"  FRAC                         %12f                                                    %12f                 %12f    \n", 
            simEventsBeforeLDCutSumAcc/simEventsBeforeLDCutSum, simPassingTotalAcc/simPassingTotal, simCountTotalAcc/simCountTotal);

    sprintf(summaryString+strlen(summaryString), "scale factor = %f \n", scaleFac);
    sprintf(summaryString+strlen(summaryString),"  TOTAL                   %12f   %12f   %12f \n", 
            simEventsBeforeLDCutSum, simPassingTotal, simCountTotal);
    sprintf(summaryString+strlen(summaryString),"  ACCEPTED                %12f   %12f   %12f \n", 
            simEventsBeforeLDCutSumAcc, simPassingTotalAcc, simCountTotalAcc);
    sprintf(summaryString+strlen(summaryString),"  FRAC                    %12f   %12f   %12f \n", 
            simEventsBeforeLDCutSumAcc/simEventsBeforeLDCutSum, simPassingTotalAcc/simPassingTotal, simCountTotalAcc/simCountTotal);
    
    sprintf(spreadSheetStr+strlen(spreadSheetStr),"%i bins fitted \n", binsFitted);
    sprintf(spreadSheetStr+strlen(spreadSheetStr),"%i bins accepted \n", binsAccepted);
    sprintf(spreadSheetStr+strlen(spreadSheetStr)," total probability = %f \n\n", totalProb);

    sprintf(summaryString+strlen(summaryString),"%i bins fitted \n", binsFitted);
    sprintf(summaryString+strlen(summaryString),"%i bins accepted \n", binsAccepted);
    sprintf(summaryString+strlen(summaryString),"total prob = %f \n", totalProb);
    sprintf(summaryString+strlen(summaryString),"--------------------\n\n");
    //sprintf(spreadSheetStr+strlen(spreadSheetStr), totalString);
    //sprintf(summaryString+strlen(summaryString), totalString);
    //lprintf(spreadSheetStr);
    FILE* spreadSheetFile;
    if ( i == 0 )
    {
      sprintf(filename, "%s/oindree_optimization_sl_%02.0f_sf_%01.6f.txt", outputDir, -cutSlope, scaleFac);
    } else {      
      sprintf(filename, "%s/oindree_optimization_final_sl_%02.0f_sf_%01.6f.txt", outputDir, -cutSlope, scaleFac);
    }
    lprintf("spreadsheet filename is %s \n\n", filename);
    spreadSheetFile = fopen(filename, "w");
    if (spreadSheetFile == NULL) {printf("error opening file\n");}
    //fprintf(spreadSheetFile, "test \n");
    fprintf(spreadSheetFile, spreadSheetStr);
    fflush(spreadSheetFile);
    fclose(spreadSheetFile);
    char formatStr[256] = "%12.1f   %12.3f   %7.2f    %7.6f    %12.6f    %7.6f    %12.6f     %7.6f    %14.6f     %7.6f     %8.6f ";
    // write the summary file entries   -------------------------
    if (firstLongCycle) {fprintf(longFile, "\n");}
    fprintf(longFile, formatStr, circPeakSepThreshold, circPeakStrengthThreshold, cutSlope, scaleFac, 
            simEventsBeforeLDCutSumAcc, simEventsBeforeLDCutSumAcc/simEventsBeforeLDCutSum, simPassingTotalAcc, simPassingTotalAcc/simPassingTotal, 
            simCountTotalAcc, simCountTotalAcc/simCountTotal, totalProb);
    if ((totalProb <= 0.1 && prevProb >= 0.1) || (totalProb >= 0.1 && prevProb <= 0.1)) {
      if (!firstLongCycle) {fprintf(longFile, "x");}
    }
    fprintf(longFile, "\n");

    if (totalProb > 0.00001 && totalProb < 0.99999) {
      if (firstShortCycle) {fprintf(shortFile, "\n");}
      fprintf(shortFile, formatStr, circPeakSepThreshold, circPeakStrengthThreshold, cutSlope, scaleFac, 
            simEventsBeforeLDCutSumAcc, simEventsBeforeLDCutSumAcc/simEventsBeforeLDCutSum, simPassingTotalAcc, simPassingTotalAcc/simPassingTotal, 
            simCountTotalAcc, simCountTotalAcc/simCountTotal, totalProb);
              if ((totalProb <= 0.1 && prevProb >= 0.1) || (totalProb >= 0.1 && prevProb <= 0.1)) {
                if (!firstShortCycle) {fprintf(shortFile, "x");}
              }
      fprintf(shortFile, "\n");
      firstShortCycle = false;
    }
    // ----------------------------------------------------------

    firstLongCycle = false;
    // write to summary file
    prevProb = totalProb;
  }
  fclose(shortFile);
  lprintf(summaryString);
  //char filename[1024];
  sprintf(filename, "%s/oindree_optSummary_%2.1f_%4.3f_sl_%02.0f.txt", outputDir, circPeakSepThreshold, circPeakStrengthThreshold, -cutSlope);
  FILE* summaryFile;  
  summaryFile = fopen(filename, "w");
  fprintf(summaryFile, summaryString);
  fclose(summaryFile);
  lprintf("Goodbye world \n");  
  fclose(logFile);
  if (drawIt)  app->Run();
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
  if (keyword.compare("--CUT_SLOPE")==0) {cutSlope=stof(parm);}
  else if (keyword.compare("--CIRC_PEAK_SEP_THRESHOLD")==0) {circPeakSepThreshold=stof(parm);}
  else if (keyword.compare("--CIRC_PEAK_STRENGTH_THRESHOLD")==0) {circPeakStrengthThreshold=stof(parm);}
  else if (keyword.compare("--PHI_HP_OFFSET")==0) {hpPhiOffset = stof(parm);}
  else if (keyword.compare("--THETA_HP_OFFSET")==0) {hpThOffset = stof(parm);}
  // TODO INPUT_FILE or INPUT_DIR   use {strcpy(inputFilePath, parm.c_str());}
  //else if (keyword.compare("--USE_ABBYS_CUT_VALUES")==0) {}
  else result = -1;
  printf("keyword parameter %s    result code %i \n", kwp, result);
  return result;
}

double likelihood0W(double arg, void* parms) {
  struct f_parms* p = (struct f_parms*)parms;
  double b = (p->b);
  double alpha = (p->alpha);
  return likelihood0(b, arg) - alpha;
}

bool sortFunc(pair<int, dataRow_t> p1, pair<int, dataRow_t> p2) {
  dataRow_t r1 = p1.second;
  dataRow_t r2 = p2.second;
  if (r1.status != r2.status) {
    return r1.status > r2.status;
  } else {
    bool pValOk1 = (r1.logLikePVal > minPVal && r1.logLikePVal != 1.0);  //add in pval==1.0 being bad?
    bool pValOk2 = (r2.logLikePVal > minPVal && r2.logLikePVal != 1.0);
    if (pValOk1 != pValOk2) {
      return (pValOk1 > pValOk2);
    } else {
      if (r1.simEventsBeforeLDCut != r2.simEventsBeforeLDCut) {
        return r1.simEventsBeforeLDCut > r2.simEventsBeforeLDCut;
      } else {
        bool bgOk1 = (r1.expBackground < 0.1);
        bool bgOk2 = (r2.expBackground < 0.1);
        return bgOk1 > bgOk2;
      }
    }
  }
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

/////////////////////////////////////////////////////////////////////
//
// Returns a random number from a lognormal distribution obtained using the dart board method.
//
// logNormal(x,mu,sigma) = 1/(x*sigma*sqrt(2pi)) * exp(-(log(x)-mu)^2/(2*sigma^2) )
//
/////////////////////////////////////////////////////////////////////
float getRandLogNormal( float mean, float var ) {

  float sigma = log( mean*mean/sqrt(var+mean*mean) );
  float mu = sqrt( log(var/(mean*mean)+1) );


  float xR = 0;
  float yR = 0;
  float yCalc = 0;
  int trys = 0; 
  while( yR < yCalc && trys < 50 ) {
    while(xR == 0){ xR = gRandom->Rndm(); };
    yR = gRandom->Rndm();

    xR *= 5*mean;       //rescale the random numbers
    yR *= 2/mean;  

    yCalc = 1/(xR*sigma*sqrt(2*M_PI)) * exp( -(log(xR)-mu)*(log(xR)-mu)/(2*sigma*sigma) );
  
    trys++;
  }
  return xR;
}

// function that returns a random background based on given exponetial fit
// parameter central values and uncertanties.
//  -requires fitSlopeAvg, fitSlopeSigma, fitInterceptAvg, fitInterceptSigma, 
//   and row (the correlation between the fit parameters).
//  -Returns random fitSlope and fitIntercept (slopeR,interceptR)
pair<float,float> getRandParamBackground ( float sAvg, float iAvg, float sSigma, float iSigma, float row )
{
  // make our new random fitSlope and fitIntercept calculate the background
  // - we use a 2d gaus and use the dartboard method.
  // - also accept only negative slopes, (ycalc)
  float xR = 0;
  float yR = 0;
  float zR = 0;
  float xCalc = 1;
  float yCalc = 1;
  float zCalc = 1;
  int tryParams = 0;
  while (zR > zCalc || yCalc >= 0)
  {
    xR = gRandom->Rndm();
    yR = gRandom->Rndm();
    zR = gRandom->Rndm();
    // find x and y such that they go out to the 5 sigma level above and below the mean
    xCalc = 5*2*iSigma*(xR-0.5)+iAvg;
    yCalc = 5*2*sSigma*(yR-0.5)+sAvg;
    zCalc = exp(-1.0/(2.0-2.0*row*row)*
               (       (xCalc-iAvg)*(xCalc-iAvg)/(iSigma*iSigma)
                +2*row*(xCalc-iAvg)*(yCalc-sAvg)/(iSigma*sSigma)
                +      (yCalc-sAvg)*(yCalc-sAvg)/(sSigma*sSigma) ));
    tryParams++;
  }

  pair<float,float> paramsR = pair<float,float>(yCalc,xCalc);
  //float slopeR = yCalc; //gRandom->Gaus(thisRow.fitSlope,thisRow.fitSlopeError);
  //float interceptR = xCalc; //gRandom->Gaus(thisRow.fitIntercept,thisRow.fitInterceptError);

  //debug stuff
  //if ( isnan(yCalc) )
  //{
  //  cout << "getRandParamBackground Error:" << endl;
  //  cout << "  sAvg,sSigma,iAvg,iSigma,row : (" << sAvg << "," << sSigma << "," << iAvg << "," << iSigma << "," << row << ")" << endl; 
  //  sleep(1);
  //}

  return paramsR;
}

// function to get spillover uncertanties
// -needs to be given spillover map, cutValHist map, hpBin, and cutval
// -returns (outflow,inflow), already randomized.

pair<float,float> getSpilloverOutIn( map<int, map<int,pair<double,int> > > spillover,  map<int, dataRow_t> cutValHist, int thisHpBin, float cutVal)
{

  map<int, dataRow_t>::iterator thisPair = cutValHist.find(thisHpBin);
  float bgNormFac = thisPair->second.cutHist->GetBinWidth(1);     //the same no matter what our sim spacing is.

  float thisSlope = thisPair->second.fitSlope;
  float thisIntercept = thisPair->second.fitIntercept;
  float thisBackground = -1.0/thisSlope/bgNormFac * exp(thisSlope * cutVal + thisIntercept);
  thisBackground *= (0.9 / sampleFrac);

  float inflow = 0;
  float outflow = 0;

  for (map<int, map<int,pair<double,int> > >::iterator thisSpill=spillover.begin(); thisSpill!=spillover.end(); ++thisSpill)
  {
    int fromHpBin = thisSpill->first;
    map<int,pair<double,int> > thisMap = thisSpill->second;
    for ( map<int,pair<double,int> >::iterator thatCont=thisMap.begin(); thatCont!=thisMap.end(); ++thatCont)
    {
      int toHpBin = thatCont->first;
      if (toHpBin == thisHpBin)
      {
        if (fromHpBin == thisHpBin)
        {
          float spill = 1.0-thatCont->second.first;
          //int spillN = thatCont->second.second;
          float spillR = abs(gRandom->Gaus(0,spill));
          outflow = spillR*thisBackground;
          if ( outflow >= thisBackground ) { outflow = thisBackground; }
        }
          else
        {
          map<int, dataRow_t>::iterator thatPair = cutValHist.find(fromHpBin);
          // make sure we actually have data and a fit for that bin
          if ( thatPair->second.status == 0 && thatPair != cutValHist.end())
          {
            float spill = thatCont->second.first;
            //float spillN = thatCont->second.second;
            float spillR = abs(gRandom->Gaus(0,spill));
            float thatSlope = thatPair->second.fitSlope;
            float thatIntercept = thatPair->second.fitIntercept;
            float thatBackground = -1.0/thatSlope/bgNormFac * exp(thatSlope * cutVal + thatIntercept);
            thatBackground *= (0.9 / sampleFrac);
            inflow += spillR*thatBackground;
          }
          // not sure what to do if we dont have a fit? or the bin has no data?
          // -If the bin has no data, just set that term to 0 (dont add anything) because
          //  that means we cant get any spill over from it
          // -If the bin was rejected for not enough bins with data then also do nothing
          //  because that means it had very few events and should hav low spillover
          // -If the event had a bad pval on the fit, just use it anyway
          // -If the fit failed and there wa plenty of data in the bin....?
          //  Im not sure, thats the case I am confused by. right not we do nothing and skip the bin.
        }
        //cout << "thisHpBin = " << thisHpBin << " toHpBin = " << toHpBin << " fromHpBin = " << fromHpBin << endl;
        //cout << " outflow = " << outflow << " inflow = " << inflow << endl;
        //cin.ignore();
      }
    }
  }
  pair<float,float> outInR = pair<float,float>(outflow,inflow);
  return outInR;
}


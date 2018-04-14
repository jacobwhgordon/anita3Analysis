///////////////////////////////////////////////////////////////////////
//
// Little file to calculate spillover.
//
// -Takes healpix orintation and index as input
//  - index is like a spacing thing. it sets skips events alot of the time.
// -Reads in analysisOutput_2.root file
// -Does rebinning on events it cares about
// -Calc spillover
// -Output to file.
//
////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"

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

#include <mpi.h>

vector<float> surfaceEllipseParms1(float phi, float theta, float dist, float errPh, float errTh,
        UsefulAdu5Pat* gps, BedmapReader* bedmap);


int processKeywordParm(char* kwp);

float cutSlope = -6.0;
float circPeakSepThreshold = 46.0;
float circPeakStrengthThreshold = 0.015;
//float sampleFrac = 0.1;
//float minPVal = 0.05;

float hpThOffset = 0.0;
float hpPhiOffset = 0.0;


int main(int argc, char* argv[]) {
  printf("Hello World:  %s \n", argv[0]);
  char filename[1024];
  int startRun = 175;
  int endRun = 439;
  char whichPol = 'H';

  TH2F* hpCountsHist = 0;
  TH2F* spilloverHist = 0;

  char inputDir[1024] = "/users/PAS0174/osu8620/anita/SamsAnalysis/results";

  bool noisyMode = false;
  bool makePlots = true;

  int pNum = 0;
  int index = 0;

  //MPI_Init(&argc, &argv);

  for (int i=1; i<argc; ++i) 
  {
    if (argv[i][0] == '-') 
    {
      if (strlen(argv[i]) > 1) 
      {
        if (argv[i][1] == 'v') {verboseMode = true;}
        else if (argv[i][1] == 'V') {verboseMode = true; noisyMode = true;}
        else if (argv[i][1] == '-') {
          processKeywordParm(argv[i]); }
      }
    } else {
      if (++pNum==1) {index = atoi(argv[i]);}
    }
  }

  //MPI_Comm_rank(MPI_COMM_WORLD, &index);

  if (index > 9)
  {
    index -= 10;
    hpThOffset += 0.56;
  }

  
  char outputDir[1024];
  sprintf(outputDir, "%s/spillover/phiOff%03.4f/thOff%03.4f", inputDir, hpPhiOffset, hpThOffset);
  //char sumOutputDir[1024];
  //sprintf(sumOutputDir, "%s/%2.1f/%1.3f", inputDir, circPeakSepThreshold, circPeakStrengthThreshold);
  printf("theta healpix offset = %f \n", hpThOffset);
  printf("phi healix offset = %f \n", hpPhiOffset);
  printf("index = %i \n", index);
  char outputDirCmd[1024];
  sprintf(outputDirCmd, "mkdir -p %s", outputDir);
  printf("output directory command is %s \n", outputDirCmd);
  system(outputDirCmd);

  BedmapReader* bedmap = BedmapReader::Instance(false); 

  char inFilename[1024]; sprintf(inFilename, "%s/analysisOutput_2_%03i_%03i.root", inputDir, startRun, endRun);
  printf("input filename is %s \n", inFilename);
  TFile* inFile = new TFile(inFilename);
  TTree* tree0 = (TTree*)inFile->Get("analysisOutput0");
  printf("input tree 0 has %lli entries \n", tree0->GetEntries());
  TTree* tree1 = new TTree("analysisOutput1", "analysis1");

  int eventNumber = 0;
  int runNumber = 0;
  char pol = '0';
  float hPeak = 0;
  float cPeak = 0;
  float cohSnr = 0;
  float cohSnr1 = 0;
  float cohSnr2 = 0;
  float linPolFrac = 0;
  float minCircPeakVal = 0;
  float cPolDist = 0;
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
  tree0->SetBranchAddress("linPolFrac", &linPolFrac);
  tree0->SetBranchAddress("minCircPeakVal", &minCircPeakVal);
  tree0->SetBranchAddress("cPolDist", &cPolDist);


  tree1->Branch("runNumber", &runNumber, "runNumber/I");
  tree1->Branch("eventNumber", &eventNumber, "eventNumber/I");
  tree1->Branch("pol", &pol, "pol/B");
  tree1->Branch("hpBin", &hpBin, "hpBin/I");
  tree1->Branch("hpWeight", &hpWeight, "hpWeight/F");


  tree0->SetBranchAddress("gpsRaw", &gpsRaw);                                         // raw payload gps info
  //tree0->SetBranchAddress("gps", &gps);                                               // payload gps info
  tree0->SetBranchAddress("theta", &theta);                                           // event theta
  tree0->SetBranchAddress("phi", &phi);                                               // event phi
  tree0->SetBranchAddress("ea", &ea);                                                 // event easting
  tree0->SetBranchAddress("no", &no);                                                 // event northing
  tree0->SetBranchAddress("lat", &lat);                                               // event latitude
  tree0->SetBranchAddress("lon", &lon);                                               // event longitude
  tree0->SetBranchAddress("alt", &alt);                                               // event altitude


  tree0->BuildIndex("eventNumber", "pol");

  TTreeIndex* tree0Index = (TTreeIndex*)tree0->GetTreeIndex();

  /////////////////////////////////////////////////////////
  //
  // Define some stuff for hists we want to make
  //
  ////////////////////////////////////////////////////////////

  if (makePlots)
  {
    char title[1024];
    char name[1024];
    sprintf(title, "Binned Continent Events: run %i-%i, %cPol;ea(km);no(km) ",
            startRun, endRun, whichPol);
    sprintf(name, "hpCountHist%1i", whichPol);
    hpCountsHist = new TH2F(name, title, 600, -3000, 3000, 600, -3000, 3000);

    sprintf(title, "Binned Spillover Events: run %i-%i, %cPol;ea(km);no(km) ",
            startRun, endRun, whichPol);
    sprintf(name, "spilloverHist%1i", whichPol);
    //spilloverHist = new TH2F(name, title, 400, -2000, 2000, 400, -2000, 2000);

    spilloverHist = new TH2F(name, title, 5000, 000, 500, 5000, -1400, -900);

  }



  ///////////////////////////////////////////
  //
  // Do the rebinning!
  //
  ///////////////////////////////////////////////////

  // new healpix for healpix stuff

  Healpix_Ordering_Scheme scheme = Healpix_Ordering_Scheme::RING;
  T_Healpix_Base<int>* healpix = new T_Healpix_Base<int>(healPixN, scheme);
  T_Healpix_Base<int>* healpix_dense = new T_Healpix_Base<int>(healPixN+3, scheme);

  //lprintf("healpix map has %i pixels \n", healpix->Npix());


  bool doRebinning = true;
  int prevEvent = 0;
  char prevPol = 0;


  float thError = 0.25;     //this might depend on SNR later?
  float phiError = 0.50;

  if ( doRebinning )
  {
    //printf("Rebinnning index %i \n",index);
    
    vector<float> errEllipseParms1 = {0, 0, 0, 0};
    vector<vector<double> > errHexVertices = vector<vector<double> >(0);
    vector<vector<double> > errHexPixels = vector<vector<double> >(0);

    for (int e = index; e<tree0Index->GetN(); e+=10)
    {
      int entNum = tree0Index->GetIndex()[e];
      tree0->GetEntry(entNum);
      if (pol != whichPol) {continue;}
      int getRes = tree0->GetEntryWithIndex(eventNumber, pol);
      bool useThisEntry = true;
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


  //////////////////////////////////////////////////////////////////
  //
  //  Spillover calc
  //
  ///////////////////////////////////////////////////////////////////

  

  // Pre calculate a vector of vectors of probablilities to save time later



  int nSteps = 100;
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


  // spillover explained!
  // (BinFrom, (BinTo, (Contribution, Count)))
  map<int,map<int,pair<double,int> > > spillover;


  // special storage for some extra info about a few events
  //  stores event number
  //   as int
  //  stores event spillover PDF  ( ea, no, p )
  //   as vector< pair< pair<float,float>,float > >
  //  we want the error ellipse too, but we will make that later
  //  we also want the event ea and no, we will get that later too

  map<int,vector<pair<pair<float,float>,float> > > spilloverEvents;

  map<int,map<int,int> > lastevent;

  // bedmap for bedmap stuff
  //BedmapReader* bedmap = BedmapReader::Instance(false);

  // for geom stuff? not sure we need this
  //AnitaGeomTool* geom = AnitaGeomTool::Instance();

  //printf("Starting Spillover Loop for index %i \n" , index);
  
  //loop over events
  for (int e=0; e<tree1Index->GetN(); e+=800) {
  //for (int e=0; e<10; ++e) {
    int entNum = tree1Index->GetIndex()[e];
    tree1->GetEntry(entNum);
    if (pol != whichPol) {continue;}
    vbprintf("event %i in hpBin %i \n", eventNumber, hpBin);
    // do circPol cuts and get info if new event
    if (eventNumber != prevEvent || pol != prevPol) {
      int getRes = tree0->GetEntryWithIndex(eventNumber, pol);
      bool useThisEntry = true;
      if (getRes <= 0) {lprintf("database error event %i pol %c bin %i - skipping \n", eventNumber, pol, hpBin); useThisEntry = false;}
      if (circPeakSepThreshold>0 && cPolDist>circPeakSepThreshold) {useThisEntry = false;}
      if (circPeakStrengthThreshold>0 && minCircPeakVal<circPeakStrengthThreshold) {useThisEntry = false;}
      if (!useThisEntry) {continue;}
      //lprintf(" event number %i %cPol \n", eventNumber, pol) ;
      //lprintf(" indexed tree get result is %i \n", getRes);
      prevEvent = eventNumber;
      prevPol = pol;
    }
    //printf("event# %i for index %i\n", e, index);
    if (spillover.find(hpBin) == spillover.end())
    {
      map<int,pair<double,int> > thisMap;
      spillover.insert(pair<int,map<int,pair<double,int> > >(hpBin, thisMap) );
      map<int,int> thisMap2;
      lastevent.insert(pair<int,map<int,int> >(hpBin, thisMap2) );
    }

    if(makePlots) {
      //if(e = ??)
      vector<pair<pair<float,float>,float> > tempVec;
      spilloverEvents.insert(pair<int,vector<pair<pair<float,float>,float> > >(eventNumber, tempVec) );
    }

    UsefulAdu5Pat* gps = new UsefulAdu5Pat(gpsRaw);

    double eventTh = theta;
    double eventPhi = phi;

    int nSteps = 100;


    //printf("event# %i for index %i start of grid loop\n", e, index);
    for (int xi = 0; xi < nSteps; xi++)
    {
      double thisTh = 0.1*(50-(float)xi)*thError + eventTh;
      for (int yi = 0; yi < nSteps; yi++)
      {
        //printf("xi=%i, yi=%i, event=%i, index=%i \n",xi,yi,e,index);
        double thisPhi = 0.1*(50-(float)yi)*phiError + eventPhi;

        // Find P for this grid point

        //double thisP = exp( -1*(thisPhi-eventPhi)*(thisPhi-eventPhi)/(2*phiError*phiError) - (thisTh-eventTh)*(thisTh-eventTh)/(2*thError*thError) );
        //thisP *= 1/(2*M_PI*phiError*thError);       // normailize to 1
        //thisP *= hpWeight;                          // normailize to hpWeight
        double thisP = probArray[xi][yi];
        if (thisP < minProb) {continue;}
        thisP *= hpWeight;

        // Find which bin this grid point falls in.

        int thisHpBin = 0;

        //  Trace back to continent to find lon and lat
        double thisLat = 0; double thisLon= 0; double thisAlt = 0; double thisAdj = 0;
        int traceStatus = gps->traceBackToContinent(thisPhi*M_PI/180.0, thisTh*M_PI/180.0, &thisLon, &thisLat, &thisAlt, &thisAdj);
        //gps->getSourceLonAndLatAtDesiredAlt(thisPhi*M_PI/180.0, thisTh*M_PI/180.0, thisLon, thisLat, alt);
        //gps->getSourceLonAndLatAtAlt(thisPhi*M_PI/180, thisTh*M_PI/180, thisLon, thisLat, 0.0);

        if (xi == 50 && yi == 50){
          //printf("event# %i for index %i middle of grid loop\n", e, index);
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

        if(makePlots)
        {
          //if(e==??)
          //get easting and northing
          double thisEa; double thisNo;
          bedmap->LonLattoEaNo(thisLon,thisLat,thisEa,thisNo);
          pair<float,float> thisEaNo = pair<float,float>(thisEa,thisNo);
          pair<pair<float,float>,float> thisG = pair<pair<float,float>,float>(thisEaNo,thisP);

          spilloverEvents[eventNumber].push_back(thisG);
        }

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
    delete gps;
    gps = 0;

    //printf("End of loop event# %i index %i \n",e,index);

    //map<int,pair<float,int> > = spillover[hpBin];    
    //for ( map<int,pair<double,int> >::iterator thisCont=thisMap.begin(); thisCont!=thisMap.end(); ++thisCont) {
    //  cout << " Bin: << hpBin << " to Bin: " 

  }

  //add a sleep statment to stager the different nodes? It might help?
  //printf("sleeping %i", index+1);
  //sleep(index);


  // We still need to normalize the spillover per bin.
  //  -Each bin should have a total of 1 in all the bins in it's map.
  //  -While we are here, lets also make an output file.

  /*

  char spillFilename[1024];
  sprintf(spillFilename, "%s/spillOverTable%c_%i.txt", outputDir, whichPol, index);
  //printf("filename: %s \n",spillFilename);
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
      //cout << "Bin " << hpBin << " contributes " << thisCont->second.first << " percent of its background to bin ";
      //cout << inHpBin << " from " << thisCont->second.second << " events." << endl;
      fprintf(spillFileOut, "%5i, %5i, %6.4f, %5i \n", hpBin, inHpBin, thisCont->second.first, thisCont->second.second);
    }
  }
  fclose(spillFileOut);
  //cin.ignore();

  */

  ////////////////////////////////////////////////////////////////
  //
  //  Lets make those two hists now
  //
  //   - One just comes from a loop over all events
  //   - The other gets created from data in the spillover loop
  //   -- Its pretty slow... so.. yeah, deal with that while testing, lol.
  //
  /////////////////////////////////////////////////////////////////

  if( makePlots )
  {
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

    // We want bin outlines too.
    TH2F* hpBinOutline = new TH2F("name", "title", 600, -3000, 3000, 600, -3000, 3000);

    //basically, just make a big huge list of how many events are in each bin.  Use a map.

    map<int,int> denseBinCount;
    map<int,int> binCount;


    for (int e = 0; e<tree0Index->GetN(); e+=1)
    {
      int entNum = tree0Index->GetIndex()[e];
      tree0->GetEntry(entNum);
      if (pol != whichPol) {continue;}  //do we want it for all pols, or just one? idk...
      int getRes = tree0->GetEntryWithIndex(eventNumber, pol);
      bool useThisEntry = true;
      if (getRes <= 0) {lprintf("database error event %i pol %c - skipping \n", eventNumber, pol); useThisEntry = false;}
      if (circPeakSepThreshold>0 && cPolDist>circPeakSepThreshold) {useThisEntry = false;}
      if (circPeakStrengthThreshold>0 && minCircPeakVal<circPeakStrengthThreshold) {useThisEntry = false;}
      if (!useThisEntry) { continue; }
      vbprintf("Mapping event %i \n", eventNumber);
      
      // we have event lat and lon saved, so we can use it here.
      double thisETheta = (-lat+90.0+hpThOffset) * M_PI/180.0;  //dont worry about shifts for now
      double thisEPhi = (lon+hpPhiOffset) * M_PI/180.0;
      pointing point(thisETheta, thisEPhi);
      int binNo_dense = healpix_dense->ang2pix(point);
      int binNo = healpix->ang2pix(point);

      map<int, int>::iterator thisEntry = binCount.find(binNo);
      if (thisEntry == binCount.end())
      {
        binCount.insert(pair<int,int>(binNo,1));
      } else {
        binCount[binNo] += 1;
      }
      map<int, int>::iterator thisEntry_dense = denseBinCount.find(binNo_dense);
      if (thisEntry_dense == denseBinCount.end())
      {
        denseBinCount.insert(pair<int,int>(binNo_dense,1));
      } else {
        denseBinCount[binNo_dense] += 1;
      }
    } 

    // now plot the events
    
    //BedmapReader* bedmap = BedmapReader::Instance(false);

    int lastHpBin = 0;
    for (int gBin = 0; gBin <hpCountsHist->GetNcells()-1; ++gBin) {
      int xBin, yBin, zBin;
      hpCountsHist->GetBinXYZ(gBin, xBin, yBin, zBin);
      float ea = hpCountsHist->GetXaxis()->GetBinCenter(xBin) * 1000;
      float no = hpCountsHist->GetYaxis()->GetBinCenter(yBin) * 1000;
      double thisLat, thisLon;
      bedmap->EaNoToLonLat(ea, no, thisLon, thisLat);
      double thisTheta = (-thisLat+90.0+hpThOffset) * M_PI/180.0;
      double thisPhi = (thisLon+hpPhiOffset) * M_PI/180.0;
      //while (thisPhi < -M_PI) thisPhi += 2*M_PI;
      pointing point(thisTheta, thisPhi);
      int binNo = healpix_dense->ang2pix(point);
      int binNoBig = healpix->ang2pix(point);

      map<int, int>::iterator thisEntry = denseBinCount.find(binNo);
      float eventCount = 0;
      
      if (thisEntry != denseBinCount.end()) {
        eventCount = denseBinCount[binNo];
      }
      hpCountsHist->Fill(ea/1000, no/1000, eventCount);
      //printf(" ea=%8.0f, no=%8.0f, lat=%6.3f, lon=%6.3f,  theta=%4.3f, phi=%4.3f,  pixelNo=%i, count=%8.1f \n",
      //        ea, no, thisLat, thisLon, thisTheta, thisPhi, pixelNo, eventCount);
      //if (count%100 == 0) {printf("  %8i entries processed", count);}

      if (binNoBig != lastHpBin)
      {
        lastHpBin = binNoBig;
        hpBinOutline->Fill(ea/1000, no/1000, 1);
      }

    }
  

    // now we have to draw the hists

    TCanvas* hpCountsCanv = 0;
    hpCountsCanv = new TCanvas("hpCountsCanv", "hpCounts", 600, 600);

    hpCountsCanv->cd(1);
    //gStyle->SetOptStat(0);
    coastLineGr->SetMarkerColor(14);
    coastLineGr->SetFillColor(18);
    coastLineGr->SetLineWidth(1);
    coastLineGr->SetFillStyle(1001);

    hpCountsHist->SetStats(kFALSE);
    hpCountsHist->Draw("COLZ");

    gStyle->SetPalette(73);
    //hpBinOutline->SetMarkerStyle(20);
    //hpBinOutline->SetMarkerColor(14);
    //hpBinOutline->SetMarkerSize(5);
    hpBinOutline->SetStats(kFALSE);
    hpBinOutline->Draw("SAME COL");

    gPad->SetLogz();
    coastLineGr->Draw("C");
    double ea, no;  // for bedmap conversion
    for (pair<int, int> thisEntry : binCount) {
      pointing thisPoint = healpix->pix2ang(thisEntry.first);
      double thisLat = 90.0 - (thisPoint.theta * 180.0/M_PI) - hpThOffset ;
      double thisLon = thisPoint.phi * 180.0/M_PI - hpPhiOffset;
      bedmap->LonLattoEaNo(thisLon, thisLat,  ea, no);
      ea /= 1000;
      no /= 1000;
      vbprintf("hp bin %i lat,lon is %f,%f  ea,no is %f,%f \n", thisEntry.first, thisLat, thisLon, ea, no);
      char hpBinText[8]; sprintf(hpBinText, "%i", thisEntry.first);
      vbprintf("drawing hp bin number at %f, %f, %f, %f \n", ea-125, no, ea+125, no+100);
      bool posQuad = (ea*no > 0);
      double textEa = posQuad ? ea : ea-75;
      double textNo = posQuad ? no-125 : no+100;
      TText* thisText = new TText(textEa, textNo, hpBinText);
      thisText->SetTextSize(0.025);
      thisText->SetTextAngle(posQuad ? 55 : -55);
      thisText->SetTextFont(42);
      thisText->Draw();
    }

    sprintf(filename, "%s/hpCounts_%03i_%03i_%c.png", outputDir, startRun, endRun, whichPol);
    hpCountsCanv->SaveAs(filename);
    sprintf(filename, "%s/hpCounts_%03i_%03i_%c.root", outputDir, startRun, endRun, whichPol);
    hpCountsCanv->SaveAs(filename);


    //Now make the plot for the spillover!
    //
    // we have the event number we want to look at, so loop over them and add the points to the histogram
    for (pair<int, vector<pair<pair<float,float>,float> > > thisPair : spilloverEvents)
    {
      int thisEventNumber = thisPair.first;
      cout << "plotting spillover PDF for event " << thisEventNumber << endl;
      vector<pair<pair<float,float>,float> > thisVec = thisPair.second;
      for (int iV=0; iV < thisVec.size(); iV++)
      {
        spilloverHist->Fill(thisVec[iV].first.first/1000,thisVec[iV].first.second/1000,thisVec[iV].second);
      }
    }

    //Now plot the thingss
    // -Event PDFs
    // -Bin Outlines
    // -Antarctica
    // - We want error ellipse too, maybe later?

    TCanvas* spillCanv = 0;
    spillCanv = new TCanvas("spillCanv", "event PDF", 800, 800);

    spillCanv->cd(1);
    //gStyle->SetOptStat(0);
    coastLineGr->SetMarkerColor(14);
    coastLineGr->SetFillColor(18);
    coastLineGr->SetLineWidth(1);
    coastLineGr->SetFillStyle(1001);

    spilloverHist->SetStats(kFALSE);
    spilloverHist->Draw("COLZ");

    gStyle->SetPalette(73);
    //hpBinOutline->SetMarkerStyle(20);
    //hpBinOutline->SetMarkerColor(14);
    //hpBinOutline->SetMarkerSize(5);
    hpBinOutline->SetStats(kFALSE);
    hpBinOutline->Draw("SAME COL");

    //gPad->SetLogz();
    coastLineGr->Draw("C");
    //double ea, no;  // for bedmap conversion
    /*
    for (pair<int, int> thisEntry : binCount) {
      pointing thisPoint = healpix->pix2ang(thisEntry.first);
      double thisLat = 90.0 - (thisPoint.theta * 180.0/M_PI) - hpThOffset ;
      double thisLon = thisPoint.phi * 180.0/M_PI - hpPhiOffset;
      bedmap->LonLattoEaNo(thisLon, thisLat,  ea, no);
      ea /= 1000;
      no /= 1000;
      vbprintf("hp bin %i lat,lon is %f,%f  ea,no is %f,%f \n", thisEntry.first, thisLat, thisLon, ea, no);
      char hpBinText[8]; sprintf(hpBinText, "%i", thisEntry.first);
      vbprintf("drawing hp bin number at %f, %f, %f, %f \n", ea-125, no, ea+125, no+100);
      bool posQuad = (ea*no > 0);
      double textEa = posQuad ? ea : ea-75;
      double textNo = posQuad ? no-125 : no+100;
      TText* thisText = new TText(textEa, textNo, hpBinText);
      thisText->SetTextSize(0.025);
      thisText->SetTextAngle(posQuad ? 55 : -55);
      thisText->SetTextFont(42);
      thisText->Draw();
    }
    */
    sprintf(filename, "%s/spilloverEventHist_%c.png", outputDir, whichPol);
    spillCanv->SaveAs(filename);
    sprintf(filename, "%s/spilloverEventHist_%c.root", outputDir, whichPol);
    spillCanv->SaveAs(filename);




  }

  //printf("Closing for index %i \n",index);

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
  if (keyword.compare("--PHI")==0) {hpPhiOffset = stof(parm);}
  else if (keyword.compare("--THETA")==0) {hpThOffset = stof(parm);}
  // TODO INPUT_FILE or INPUT_DIR   use {strcpy(inputFilePath, parm.c_str());}
  //else if (keyword.compare("--USE_ABBYS_CUT_VALUES")==0) {}
  else result = -1;
  printf("keyword parameter %s    result code %i \n", kwp, result);
  return result;
}


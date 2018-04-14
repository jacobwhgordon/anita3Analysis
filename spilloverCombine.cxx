// This just combines the partial spillover files into one single file, also re averges.
// need to run it after running spillover and before running the optimization. 
//  (spillover is just split up because its really really slowwww)

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

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>



int main(int argc, char** argv) 
{
  char inputDir[1024] = "../results";
  char pol = 'H';

  // loop over phi and theta.
  for (int it=0; it<10; it++)
  {
    float theta = -0.56*it;

    for (int ip = 0; ip<10; ip++)
    {
      float phi = 0.56*ip;

      if (ip == 0 && it == 0)
      {
        theta = -0.0;
      }
      else if (it == 0) { theta = 0.0; }

      char fileDir[1024];
      sprintf(fileDir, "%s/spillover/phiOff%03.4f/thOff%03.4f", inputDir, phi, theta);
      
      // spillover explained!
      // (BinFrom, (BinTo, (Contribution, Count)))
      map<int,map<int,pair<double,int> > > spillover;

      // get the ten files and add up the things
      for (int i=0; i<10; i++) {
        int index = i;
        char spillFilename[1024];
        sprintf(spillFilename, "%s/spillOverTable%c_%i.txt", fileDir, pol, index);
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

            //spillStr >> hpBin >> inHpBin >> cont >> numEvents;
            //printf("bin #%i  inBin=%i  cont=%f  numEvents=%i \n", hpBin, inHpBin, cont, numEvents);
 
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
      //sleep(1);
      char spillFilename[1024];
      sprintf(spillFilename, "%s/spillover%c.txt", fileDir, pol);
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
    }
  }
}

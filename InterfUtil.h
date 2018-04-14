
// convenience definitions for conditional message output
#define vcout if(verboseMode) cout
#define vbprintf if(verboseMode) printf
#define lprintf(format, args...)  {printf(format, ## args);  fprintf(logFile, format, ## args); }

//#define vbprintlf(format, args...)      // you must add continuation backslashes if you uncomment this macro
//          if (verboseMode) { 
//            printf(format, ## args);  
//            fprintf(logFile, format, ## args);  
//          }

#define HEADER_FILENAME_ID "timedHeadFile"
//#define HEADER_FILENAME_ID "decBlindHeadFileV1_"
//#define HEADER_FILENAME_ID "decimatedHeadFile"

//int[32] qualityCutMask = {0};
//int[32] analysisCutMask = {0};
#include <vector>
#include <math.h>
#include <iostream>

#include "AnitaConventions.h"
#include "UsefulAdu5Pat.h"
#include "TMath.h"

using namespace std;

const float EARTH_POLAR_RADIUS = 6378000;

bool verboseMode = false;
FILE* logFile;

float sunCutDistPhi = 5.0;
float sunCutDistTheta = 5.0;
float sunReflCutDistPhi = 5.0;
float sunReflCutDistTheta = 5.0;

int minTrigPhisectors = 2;
int maxTrigPhisectors = 6;
float surfSatThreshold = 1000;
float dcOffsetThreshold = 50;
int minWfLength = 240;
float nadirNoiseThreshold = 0.5;
float minPeakTheta = -35.0;        // Abby's = -35      Brians = -35     my original=-35
float maxPeakTheta = -6.0;        // Abby's = none      Brians =-6     
float maxHwAngle = 33.75;        
float minLinPolFrac = 0.3;      // Abby's = 0.3      Brian's = 0.3         my original = 0.3
float minAnyPolFrac = 0.5;       
float maxCircPolFrac = 0.5;
float corrPeakThresh = 0.040;     // Abby's = 0.075   Brian's = 0.075      my original = 0.04
float hilbPeakThresh = 25;       // Abby's = 15;    Brian's = 15           my original = 25
float corrPeakDiagThresh = 0.08;  // TODO change to slope-intercept form
float hilbPeakDiagThresh = 80;
float peakRatioThreshold = 0.9;  // Abby's = 0.9    Brian's = 0.85         my original = 0.9
float snrThreshold = 3.0;        // 3.0 is arbitrary
int trigPhisecProx = 1;
int nsDiffMin = -1000;
int nsDiffMax = 1500;
float cPolPeakSeparationThreshold = 20.0;
float cPolPeakRatioThreshold = 2.0;
float cPolLesserPeakThreshold = 0.0;
float cPeakStrengthThreshold = 0.5;
float cPolNumPeaks = 2;

float ldSlope = -38.;
float ldInterceptThreshold = 0;   // unsupported
int healPixN = 4;

float southCutDist = 45;
float wedgeCutMid = 0;
float wedgeCutEdge = 60;
float minSnr = 3.0;
float misRecDistPhi = 2.5;
float misRecDistTheta = 1.0;
// TODO enable multiple rectangle cuts
// TODO enable elliptical cut
const float skyRectCut[4] = {184, -8, 192, -6};
float simENuCut[2] = {19.5, 20.5};
double sineSubThreshold = 0.075;
int sineSubFreqBands = 1;
int timeDomainPadFactor = 1;
float deadTimeThreshold = 0.9;
float minCircPeakThreshold = 0.5;

vector<double> MCMURDO = {-77.8639, 167.0368};
vector<double> VOSTOK = {-78.46, 106.84};
vector<double> WAIS_DIVIDE = {AnitaLocations::LATITUDE_WAIS, AnitaLocations::LONGITUDE_WAIS};
vector<double> SOUTH_POLE = {-90, 0};
vector<double> HALLEY = {-75.58, -26.66};
vector<double> DAVIS = {-68.57, 77.97};
vector<double> HOT_SPOT_1 = {-76, -100};


int pulserNsLdb2 = 1500;       //LDB seavey vPol 6kV
int pulserNsLdb1 = 50001575;   //LDB seavey hPol 10kV
int pulserNsLdb3 = 50001575;   //LDB seavey hPol 10kV
int pulserNsWais = -40;        //WAIS hpol

const int FILTER_OPTION_NONE = 0;
const int FILTER_OPTION_ABBY = 1;
const int FILTER_OPTION_SINE_SUBTRACT = 2;
const int FILTER_OPTION_ADAPTIVE_NOTCH = 3;
const int FILTER_OPTION_GEOMETRIC = 4;

//float ellipCutCenterLatLon[2] = {(float)WAIS_DIVIDE[0], (float)WAIS_DIVIDE[2]};
//float ellipCutRadii[2] = {50, 100};
//float ellipCutAngle = 2;

// for normal error ellipse distribution
//vector<float> pointingErr = {0.25, 0.5}; // the other numbers were backwards?  fixed it!
vector<float> pointingErr = {0.5, 0.25}; // assumed pointing error in phi/theta, for drawing error ellipses
int errorPixelN = 4;

// for "no" error ellipse distribution (and btw does not work)
//vector<float> pointingErr = {0.005, 0.0025}; // assumed pointing error in theta/phi, for drawing error ellipses
//int errorPixelN = 1;

vector<float> surfaceEllipseParms0(float theta, float dist, float alt, float errPh, float errTh) {
// project error cone on the tangent plane at the localization point
// tangent plane is taken to a sphere of polar earth radius, at payload-to-peak separation angle (theta) 
  float altF = abs(dist * tan(theta*M_PI/180.0));
  float distNear = abs(altF / tan((theta + errTh)*M_PI/180.0));
  //float distMid = abs(altF / tan(theta*M_PI/180.0));
  float distFar = abs(altF / tan((theta - errTh)*M_PI/180.0));
  float centerD = 0.5 * (distFar+distNear)-dist;
  float rAxis = 0.5 * (distFar-distNear);
  //printf("  distFar = %f, distNear = %f \n", distFar, distNear);
  float thAxis = dist * errPh*M_PI/180.0;   // assuming small angle
  // deliver source-to-center adjustment, and semi-axes
  return vector<float>({centerD, rAxis, thAxis});
}

vector<vector<double> > hexagonVertices(float x, float y, float rAxis, float thAxis, float phi) {
  const float eta = 1.1;
  const float sin60 = sqrt(3)*0.5;
  vector<vector<double> > p(0);
  // set up the points origin-centered and axis-aligned
  p.push_back({eta*rAxis, eta*rAxis*0.5, -eta*rAxis*0.5, -eta*rAxis, -eta*rAxis*0.5, eta*rAxis*0.5});
  p.push_back({0, eta*thAxis*sin60, eta*thAxis*sin60, 0, -eta*thAxis*sin60, -eta*thAxis*sin60});
  //for (int k=0; k<6; ++k) {printf("       %f, %f", p[k][0], p[k][1]);}

  // rotate through phi and translate to x,y
  for (int k=0; k<6; ++k) {
    float thisR = sqrt(p[0][k]*p[0][k] + p[1][k]*p[1][k]);
    float thisTh = atan2(p[1][k], p[0][k]);
    thisTh += phi;
    p[0][k] = thisR*cos(thisTh) + x;
    p[1][k] = thisR*sin(thisTh) + y;
  }
  return p;  
}

// this version shoots every point on to the continent
vector<vector<double> > errorPixelsRayTrace(float phi, float theta, float errPhi, float errTheta, UsefulAdu5Pat* gps, BedmapReader* bedmap) {
  const float eta = 1.2;
  int n=errorPixelN;
  vector<vector<double> > p(3);
  vector<vector<double> > q(3);
  // set up the points origin-centered, unit-scaled and axis-aligned 
  float pdfSum = 0;
  for (int m=-n; m<=n; m+=1) {
    for (int k=-(2*n-abs(m)); k<=(2*n-abs(m)); k+=2) {
      float tx = 0.5*k/(float)(n);
      p[0].push_back(tx);
      float ty = (float)(m)/(float)n *sqrt(3)*0.5;
      p[1].push_back(ty);
      float dist2 = tx*tx + ty*ty;
      float pdfRaw = exp(-dist2);
      p[2].push_back(pdfRaw);
      pdfSum += pdfRaw;
    }
  }  
  // rescale and translate to phi,theta 
  for (int k=0; k<p[0].size(); ++k) {p[0][k] *= errPhi*eta; p[1][k] *= errTheta*eta;}
  for (int k=0; k<p[0].size(); ++k) {p[0][k] += phi; p[1][k] += theta;}
  for (int k=0; k<p[0].size(); ++k) {p[2][k] /= pdfSum;}
  // shoot each point to the continent 
  for (int k=0; k<p[0].size(); ++k) {
    double thisLon, thisLat, thisAlt, thisAdj;
    gps->traceBackToContinent(p[0][k] *M_PI/180.0, p[1][k] *M_PI/180.0, &thisLon, &thisLat, &thisAlt, &thisAdj);
    // convert to Ea, No
    double thisEa, thisNo;
    bedmap->LonLattoEaNo(thisLon, thisLat, thisEa, thisNo);
    //q[0].push_back(thisLat);
    //q[1].push_back(thisLon);
    q[0].push_back(thisEa);
    q[1].push_back(thisNo);
    q[2].push_back(p[2][k]);
  }
  return q;
}

vector<vector<double> > ellipsePixels(float x, float y, float rAxis, float thAxis, float phi) {
  const float eta = 1.2;
  int n=errorPixelN;
  vector<vector<double> > p(3);
  // set up the points origin-centered, unit-scaled and axis-aligned 
  float pdfSum = 0;
  for (int m=-n; m<=n; m+=1) {
    for (int k=-(2*n-abs(m)); k<=(2*n-abs(m)); k+=2) {
      float tx = 0.5*k/(float)(n);
      p[0].push_back(tx);
      float ty = (float)(m)/(float)n *sqrt(3)*0.5;
      p[1].push_back(ty);
      float dist2 = tx*tx + ty*ty;
      float pdfRaw = exp(-dist2);
      p[2].push_back(pdfRaw);
      pdfSum += pdfRaw;
    }
  }
  //printf("untransformed pixel points: \n"); for (int k=0; k<p[0].size(); ++k) {printf("   %f, %f \n", p[0][k], p[1][k]);}}
  // scale to axis size
  for (int k=0; k<p[0].size(); ++k) {p[0][k] *= rAxis*eta; p[1][k] *= thAxis*eta;}
  // rotate through phi and translate to x,y
  for (int k=0; k<p[0].size(); ++k) {
    float thisR = sqrt(p[0][k]*p[0][k] + p[1][k]*p[1][k]);
    float thisTh = atan2(p[1][k], p[0][k]);
    thisTh += phi;
    p[0][k] = thisR*cos(thisTh) + x;
    p[1][k] = thisR*sin(thisTh) + y;
  }
  // normalize the PDF coefficients
  for (int k=0; k<p[0].size(); ++k) {p[2][k] /= pdfSum;}
  return p;
}

float reflectionAngle(float h, float el) {
  float eld = -9999;
  const float EARTH_POLAR_RADIUS = 6378000;
  if (el < -5.0) {
    float tr = M_PI / 2.0 - (el * M_PI / 180.0);
    float r = EARTH_POLAR_RADIUS;
    float rh = r + h;
    float z = cos(tr);
    float d = -rh*z - sqrt(rh*rh*z*z - 2*r*h - h*h);
    float a = d * sin(tr) / r;
    float td = M_PI + 2*a -tr;
    eld = (90.0 - 180.0*td/M_PI);
  }
  return eld;
}

float likelihood0(float b, float a) {
  float term1 = TMath::Gamma(1+b, a+b);
  float term2 = TMath::Gamma(1+b, b);
  return (term1 - term2)/(1.0 - term2);
}

float vectorAngle(float phi0, float phi1, float theta0, float theta1) {
  float ph0 = phi0 * M_PI / 180.;
  float ph1 = phi1 * M_PI / 180.;
  float th0 = (theta0 + 90.) * M_PI / 180.;
  float th1 = (theta1 + 90.) * M_PI / 180.;
  float ang = acos(cos(th0)*cos(th1) + sin(th0)*sin(th1)*cos(ph0-ph1));
  ang *= (180. / M_PI);
  return ang;
}


//Here we put functions that handle fitting and error evaluation for 
//the Anita 3 Binned analysis.
//
//By: Jacob Gordon (Unless there are problems... then... who 
//                  knows who wrote this)
//
///////////////////////////////////////////////////////////////////////

//includes

#include <cassert>
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnContours.h"
#include "Minuit2/MnPlot.h"

#include "fitFunc.h"

using namespace ROOT::Minuit2;

//functions

//actual def of the operator minuit calls
double FitFCN::operator()(const std::vector<double>& par) const {

  assert(par.size() == 2);
  FitFunction fit(par[0],par[1]);

  double ll_val = 0.;
  for (unsigned int n=0; n < theMeasurements.size(); n++) 
  {
    /* // dropping the exponential term for minimization reasons
    ll_val += (-1)*theMeasurements[n]*
                 (fit(thePositions[n])
                  +theMeasurements[n]*log(fit(thePositions[n]))
                  );
    */
    

    /* what I think it should be, all terms
    ll_val += (-1)*theMeasurements[n]*
                 (fit(thePositions[n])
                  +theMeasurements[n]*log(fit(thePositions[n]))
                  -lgamma(theMeasurements[n]+1) );
    */
    //made negative so we can minimize it.

    // what root is doing (check!)
    if (theMeasurements[n] == 0)
    {
      ll_val += (fit(thePositions[n])-theMeasurements[n]);
    } else {
      ll_val += (fit(thePositions[n])-theMeasurements[n]+theMeasurements[n]*(log(theMeasurements[n]/fit(thePositions[n]))));
    }
    

    /* // what root does *w
    if (theMeasurements[n] == 0)
    {
      ll_val += theMeasurements[n]*(fit(thePositions[n])-theMeasurements[n]);
    } else {
      ll_val += theMeasurements[n]*(fit(thePositions[n])-theMeasurements[n]+theMeasurements[n]*(log(theMeasurements[n]/fit(thePositions[n]))));
    }
    */

  }

  return ll_val;
}


//Make one to try fitting a power law
double PowerFCN::operator()(const std::vector<double>& par) const {

  assert(par.size() == 2);
  PowerFunction fit(par[0],par[1]);

  double ll_val = 0.;
  for (unsigned int n=0; n < theMeasurements.size(); n++)
  {
    if (theMeasurements[n] == 0)
    {
      ll_val += (fit(thePositions[n])-theMeasurements[n]);
    } else {
      ll_val += (fit(thePositions[n])-theMeasurements[n]+theMeasurements[n]*(log(theMeasurements[n]/fit(thePositions[n]))));
    }
  }
  return ll_val;
}




// lets also make a chi2FCN
double Chi2FCN::operator()(const std::vector<double>& par) const {

  assert(par.size() == 2);
  FitFunction fit(par[0],par[1]);

  double chi2 = 0.;
  for (unsigned int n=0; n < theMeasurements.size(); n++)
  {
    chi2 += (theMeasurements[n]-fit(thePositions[n]))*(theMeasurements[n]-fit(thePositions[n]))/(theMVariances[n]*theMVariances[n]);
  }
  return chi2;
}



//////////////////////////////////
// fuctions for local utility/////
//////////////////////////////////

// quick easy fuction to return the poisson probability for a given weight and prodiction
double FitFunc::poissonProb(double data, double model)
{
  return exp(data*log(model)-model-lgamma(data+1));
}


// maybe a function to do the stuff needed to get a minimum?
// -takes in the fitFCN we are working with, and starting values for the yint and slope.
// -finds and prints the minimum 

FunctionMinimum FitFunc::findFit( FitFCN theFCN, MnUserParameters upar ) 
{
  // create 
  MnMigrad migrad(theFCN, upar);

  // minimize
  FunctionMinimum min = migrad();  

  // output?
  //std::cout << "minimum: " << min << std::endl;
  
  // we also need to get back the fit parameter? they might be saved in theFCN?
  return min;
}

FunctionMinimum FitFunc::findFit( Chi2FCN theFCN, MnUserParameters upar )
{
  // create
  MnMigrad migrad(theFCN, upar);

  // minimize
  FunctionMinimum min = migrad();

  // output?
  //std::cout << "minimum: " << min << std::endl;

  // we also need to get back the fit parameter? they might be saved in the$
  return min;
}

FunctionMinimum FitFunc::findFit( PowerFCN theFCN, MnUserParameters upar )
{
  // create
  MnMigrad migrad(theFCN, upar);

  // minimize
  FunctionMinimum min = migrad();

  // output?
  //std::cout << "minimum: " << min << std::endl;

  // we also need to get back the fit parameter? they might be saved in the$
  return min;
}


// function to fill the covariance matrix
// currently deviding by pos.size to normalize (1/N factor) but does pos.size() == N??
void FitFunc::findCovM( double & M11, double & M12, double & M21, double & M22, std::vector<double> pos, std::vector<double> meas, double yint, double slope )
{
  FitFunction fit(yint, slope);
  
  for (int i=0; i<pos.size(); i++)
  {
    //use the numerically stable version of the poisson dist.
    double Prob = poissonProb(meas[i],fit(pos[i]));
    //std::cout << Prob << std::endl;
    M11 += meas[i]*meas[i]*(meas[i]+fit(pos[i]))*(meas[i]+fit(pos[i]))*Prob / pos.size();
    M12 += pos[i]*meas[i]*meas[i]*(meas[i]+fit(pos[i]))*(meas[i]+fit(pos[i]))*Prob / pos.size();
    M21 += pos[i]*meas[i]*meas[i]*(meas[i]+fit(pos[i]))*(meas[i]+fit(pos[i]))*Prob / pos.size();
    M22 += pos[i]*pos[i]*meas[i]*meas[i]*(meas[i]+fit(pos[i]))*(meas[i]+fit(pos[i]))*Prob / pos.size();
  }
}


// fuction takes in parameters and returns there errors
// p0 is the yint is a
// p1 is the slope is b
void FitFunc::findErrors( double & p0error, double & p1error, std::vector<double> pos, std::vector<double> meas, std::vector<double> var, double p0, double p1 )
{
  FitFunction fit(p0, p1);

  double top0 = 0;
  double bottom0 = 0;
  double top1 = 0;
  double bottom1 = 0;
  double err = 0;
  for (int i=0; i<pos.size(); i++)
  {

    /* //for James's extra factor of w_i
    top0    += meas[i]*meas[i]*meas[i]*fit(pos[i]);
    bottom0 += meas[i]*meas[i]*fit(pos[i]);
    top1    += meas[i]*meas[i]*meas[i]*fit(pos[i])*pos[i]*pos[i];
    bottom1 += meas[i]*meas[i]*fit(pos[i])*pos[i]*pos[i];
    */

    //for no extra factor of w_i, and... eq on pg 195
    top0    = 1;
    bottom0 += fit(pos[i]);
    top1    = 1;
    bottom1 += fit(pos[i])*pos[i]*pos[i];

    err += 1/(var[i]*var[i]);

  }
  /* //for James's extra w_i
  p0error = sqrt( pos.size()*top0/(bottom0*bottom0) );
  p1error = sqrt( pos.size()*top1/(bottom1*bottom1) );
  */

  p0error = sqrt(1/bottom0);
  p1error = sqrt(1/bottom1);

  //normalize ??
  p0error *= sqrt(pos.size());
  p1error *= sqrt(pos.size());

}


// fuction takes in parameters and returns there errors
// p0 is the yint is a
// p1 is the slope is b
void FitFunc::findChi2Errors( double & p0error, double & p1error, std::vector<double> pos, std::vector<double> meas, std::vector<double> var, double p0, double p1 )
{
  FitFunction fit(p0, p1);

  double top0 = 0;
  double bottom0 = 0;
  double top1 = 0;
  double bottom1 = 0;
  double err = 0;
  for (int i=0; i<pos.size(); i++)
  {
    top0 = 1;
    top1 = 1;

    /* // was this wrong?... idk...
    bottom0 += 0.5*fit(pos[i]);
    bottom1 += 0.5*fit(pos[i])*pos[i]*pos[i];
    */

    bottom0 += (2*fit(pos[i])*fit(pos[i])-fit(pos[i])*meas[i])/meas[i];
    bottom1 += pos[i]*pos[i]*(2*fit(pos[i])*fit(pos[i])-fit(pos[i])*meas[i])/meas[i];

    err += 1/(var[i]*var[i]);

  }
  p0error = sqrt( 1/(bottom0) ); 
  p1error = sqrt( 1/(bottom1) );

  //normalize by number in sum and bin size?
  p0error *= sqrt(pos.size());
  p1error *= sqrt(pos.size());


}


 


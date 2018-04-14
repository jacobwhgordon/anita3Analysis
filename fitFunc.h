//header file so minuit has that stuff it needs
#ifndef MN_FitFunction_H_
#define MN_FitFunction_H_

#include <math.h>

//model we are giving minuit, y(x)=e^(a+b*x)
class FitFunction {

public:

  FitFunction(double a, double b) :
    theA(a), theB(b) {}

  ~FitFunction() {}

  double a() const {return theA;}
  double b() const {return theB;}

  double operator()(double x) const { return exp(a()+b()*x); }

private:

  double theA;
  double theB;

};
#endif //  MN_FitFunction_H_

// Alternate model to pass minuit, a power law: y(x) = ax^(b)
class PowerFunction {

public:

  PowerFunction(double a, double b) :
    theA(a), theB(b) {}

  ~PowerFunction() {}

  double a() const {return theA;}
  double b() const {return theB;}

  double operator()(double x) const { return a()*pow(x,b()); }

private:

  double theA;
  double theB;

};



//now a thing for the FCN minuit wants

#ifndef MN_FitFcn_H_
#define MN_FitFcn_H_

#include "Minuit2/FCNBase.h"   //note, currently referencing the Minuit2 .h file, might need to be the Minuit one?

#include <vector>

class FitFCN : public ROOT::Minuit2::FCNBase {

public:

  FitFCN(
    const std::vector<double>& meas,
    const std::vector<double>& pos,
    const std::vector<double>& mvar) :
    theMeasurements(meas), thePositions(pos),
    theMVariances(mvar), theErrorDef(0.5) {}

  ~FitFCN() {}


  virtual double Up() const {return theErrorDef;}
  virtual double operator()(const std::vector<double>&) const;

  std::vector<double> measurements() const {return theMeasurements;}
  std::vector<double> positions() const {return thePositions;}
  std::vector<double> variances() const {return theMVariances;}

  void setErrorDef(double def) {theErrorDef = def;}

private:

  std::vector<double> theMeasurements;
  std::vector<double> thePositions;
  std::vector<double> theMVariances;
  double theErrorDef;
};

#endif //MN_FitFcn_H_ 


#ifndef MN_PowerFcn_H_
#define MN_PowerFcn_H_

#include "Minuit2/FCNBase.h"   //note, currently referencing the Minuit2 .h file, might need to be the Minuit one?
#include <vector>

class PowerFCN : public ROOT::Minuit2::FCNBase {

public:

  PowerFCN(
    const std::vector<double>& meas,
    const std::vector<double>& pos,
    const std::vector<double>& mvar) :
    theMeasurements(meas), thePositions(pos),
    theMVariances(mvar), theErrorDef(0.5) {}

  ~PowerFCN() {}


  virtual double Up() const {return theErrorDef;}
  virtual double operator()(const std::vector<double>&) const;

  std::vector<double> measurements() const {return theMeasurements;}
  std::vector<double> positions() const {return thePositions;}
  std::vector<double> variances() const {return theMVariances;}

  void setErrorDef(double def) {theErrorDef = def;}

private:

  std::vector<double> theMeasurements;
  std::vector<double> thePositions;
  std::vector<double> theMVariances;
  double theErrorDef;
};

#endif //MN_PowerFcn_H_


#ifndef MN_Chi2Fcn_H_
#define MN_Chi2Fcn_H_

#include "Minuit2/FCNBase.h"   //note, currently referencing the Minuit2 .h file, might need to be the Minuit one?

#include <vector>

class Chi2FCN : public ROOT::Minuit2::FCNBase {

public:

  Chi2FCN(
    const std::vector<double>& meas,
    const std::vector<double>& pos,
    const std::vector<double>& mvar) :
    theMeasurements(meas), thePositions(pos),
    theMVariances(mvar), theErrorDef(1.0) {}

  ~Chi2FCN() {}


  virtual double Up() const {return theErrorDef;}
  virtual double operator()(const std::vector<double>&) const;

  std::vector<double> measurements() const {return theMeasurements;}
  std::vector<double> positions() const {return thePositions;}
  std::vector<double> variances() const {return theMVariances;}

  void setErrorDef(double def) {theErrorDef = def;}

private:

  std::vector<double> theMeasurements;
  std::vector<double> thePositions;
  std::vector<double> theMVariances;
  double theErrorDef;
};

#endif //MN_FitFcn_H_


//now some more random local functions

#ifndef FitFunc_H_
#define FitFunc_H_

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
using namespace ROOT::Minuit2;

class FitFunc {


public:

  FunctionMinimum findFit(FitFCN, MnUserParameters);
  FunctionMinimum findFit(Chi2FCN, MnUserParameters);
  FunctionMinimum findFit(PowerFCN, MnUserParameters);
  void findCovM(double &, double &, double &, double &, std::vector<double>, std::vector<double>, double, double);
  void findErrors(double &, double &, std::vector<double>, std::vector<double>, std::vector<double>, double, double);
  void findChi2Errors(double &, double &, std::vector<double>, std::vector<double>, std::vector<double>, double, double);
  double poissonProb(double, double);

private:

};

#endif //FitFunc_H_

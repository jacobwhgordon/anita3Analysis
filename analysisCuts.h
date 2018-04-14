using namespace std;

#include <vector>

const int PRE_Q_CUT_TRIGGER_TYPE = 0;
const int PRE_Q_CUT_PAYLOAD_BLAST = 1;
const int PRE_Q_CUT_CAL_PULSER = 2;
const int PRE_Q_CUT_SAT_STRIPE = 3;
const int PRE_Q_CUT_CPOL_SEP = 4;
const int PRE_Q_CUT_NUM_PHISECTOR = 5;
const int PRE_Q_CUT_WF_SATURATION = 6;
const int PRE_Q_CUT_WF_DC_OFFSET = 7;
const int PRE_Q_CUT_SHORT_WF = 8;

const int Q_CUT_TRIGGER_TYPE = 0;
const int Q_CUT_NUM_PHISECTOR = 1;
const int Q_CUT_PAYLOAD_BLAST = 2;
const int Q_CUT_WF_SATURATION = 3;
const int Q_CUT_WF_DC_OFFSET = 4;
const int Q_CUT_SHORT_WF = 5;
const int Q_CUT_NADIR_NOISE = 6;
const int Q_CUT_NO_TRIGGER = 7;

const int A_CUT_SOLAR_DIRECT = 1;
const int A_CUT_SOLAR_REFL = 2;
const int A_CUT_GEOSTAT_DIRECT = 3;
const int A_CUT_GEOSTAT_REFL = 4;
//const int A_CUT_THERMAL = 5;
const int A_CUT_CONTINENT = 6;
const int A_CUT_THETA = 11;
const int A_CUT_HW_ANGLE = 12;
const int A_CUT_PEAK_RATIO = 13;
const int A_CUT_LINPOL_FRAC = 14;
const int A_CUT_CORR_PEAK = 15;
const int A_CUT_HILBERT_PEAK = 16;
const int A_CUT_DIAGONAL = 17;     // hilbert peak / correlation peak diagonal
const int A_CUT_WAIS = 18;
const int A_CUT_LDB = 19;
const int A_CUT_HICAL = 20;
const int A_CUT_SOUTH = 21;
const int A_CUT_ANYPOL_FRAC = 22;
const int A_CUT_CIRCPOL_FRAC = 23;
const int A_CUT_SNR = 24;  
const int A_CUT_WEDGE = 25;        // centered wedge on payload-north sky
const int A_CUT_RECT = 26;         // rectangle in the payload-north sky
const int A_CUT_CONT_CIRC = 27;    // continent dist from a point
const int A_CUT_SKY_CIRC = 28;     //  sky pl-coords dist from a point
const int A_CUT_L3_TRIG_DIR = 29;     //  sky pl-coords dist from a point
const int A_CUT_SIM_E_NU = 30;     //  simulated neutrino energy
const int A_CUT_CM3 = 31;
const int A_CUT_CPOL_PEAK_SEPARATION = 32;
const int A_CUT_CPOL_PEAK_RATIO = 33;
const int A_CUT_CPOL_LESSER_PEAK = 34;
const int A_CUT_LPOL_THETA = 35;
const int A_CUT_DEADTIME = 36;
const int A_CUT_CPOL_PEAK_STRENGTH = 37;
const int A_CUT_HP_BIN_REJECTED = 38;
//const int A_CUT_BIN_BACKGROUND = 39;
//const int A_CUT_BINS_IN_FIT = 40;
//const int A_CUT_EVENTS_IN_FIT_REGION = 41;
const int A_CUT_SATALLITE_AREA = 42;

const char NUM_PRE_Q_CUTS = 32;
const char NUM_Q_CUTS = 32;
const char NUM_A_CUTS = 64;

const char PRE_Q_CUT_DESC[NUM_PRE_Q_CUTS][32] =
                                {"trigger type",                //0
                                 "payload blast",               //1
                                 "cal pulser",                  //2
                                 "satellite stripe",            //3
                                 "cpol peak seperation",        //4
                                 "L3 triggering phisectors",    //5
                                 "waveform saturation",         //6
                                 "waveform dc offset",          //7
                                 "short waveform",              //8
                                 "","","","","",""};


const char Q_CUT_DESC[NUM_Q_CUTS][32] = {"trigger type",        // 0
                                 "L3 triggering phisectors",    // 1
                                 "payload blast",               // 2
                                 "waveform saturation",         // 3
                                 "waveform dc offset",          // 4
                                 "short waveforms",             // 5
                                 "nadir noise",                 // 6
                                 "no trigger",                  // 7
                                 "","","","","","","",""};
// there's a better way
const char A_CUT_DESC[NUM_A_CUTS][32] = {"undefined",                   // 0
                                         "solar direct",                // 1
                                         "solar reflection",            // 2
                                         "geostationary (inactive)",    // 3
                                         "geostat refl (inactive)",     // 4
                                         "thermal (deprecated)",        // 5
                                         "reconstruct to continent",    // 6
                                         "undefined",                   // 7
                                         "undefined",                   // 8
                                         "undefined",                   // 9
                                         "undefined",                   // 10
                                         "elevation angle",             // 11
                                         "hw angle",                    // 12
                                         "2nd/1st peak ratio",          // 13
                                         "linear pol fraction",         // 14
                                         "correlation peak",            // 15
                                         "hilbert peak",                // 16
                                         "linear discriminant cut",     // 17
                                         "wais pulse",                  // 18
                                         "LDB pulse (inactive)",        // 19
                                         "HICAL pulse (inactive)",      // 20
                                         "south pointing",              // 21
                                         "any pol fraction",            // 22
                                         "circ pol fraction",           // 23
                                         "snr",                         // 24
                                         "wedge",                       // 25
                                         "rectangle",                   // 26
                                         "continent circle",            // 27
                                         "sky circle",                  // 28
                                         "hardware trigger direction",  // 29
                                         "simulated neutrino energy",   // 30      
                                         "central moment 3",            // 31
                                         "cpol peak separation",        // 32
                                         "cpol peak ratio",             // 33
                                         "cpol lesser peak",            // 34
                                         "lpol theta",                  // 35
                                         "high deadtime",               // 36
                                         "cpol peak strength",          // 37
                                         "healpix bin rejected",        // 38
                                         "bin background",              // 39
                                         "background fit",              // 40
                                         "events in fit region",        // 41
                                         "satallite stripe cut"         // 42
                                        };

vector<int> aCutOrderStage1 = { 
                                    A_CUT_SOLAR_REFL,
                                    A_CUT_CONTINENT,
                                    A_CUT_THETA,
                                    A_CUT_L3_TRIG_DIR,
                                    A_CUT_WAIS,
                                    A_CUT_LDB
                              };                    

vector<int> aCutOrderStage2 =  {  // these cuts will be calculated and applied in this program
                                    A_CUT_PEAK_RATIO,
                                    A_CUT_CORR_PEAK,
                                    A_CUT_HILBERT_PEAK,
                                    A_CUT_LINPOL_FRAC,
                                    A_CUT_DEADTIME,
                                    A_CUT_CPOL_PEAK_SEPARATION,
                                    A_CUT_CPOL_PEAK_STRENGTH
                                    // A_CUT_DIAGONAL
                               };

vector<int> qCutOrder = {
                          Q_CUT_TRIGGER_TYPE, 
                          Q_CUT_NUM_PHISECTOR,
                          Q_CUT_PAYLOAD_BLAST,
                          Q_CUT_WF_SATURATION,
                          Q_CUT_WF_DC_OFFSET,
                          Q_CUT_SHORT_WF,
                          Q_CUT_NADIR_NOISE,
                          Q_CUT_NO_TRIGGER
                        };

vector<int> preQCutOrder = {
                             PRE_Q_CUT_TRIGGER_TYPE,
                             PRE_Q_CUT_PAYLOAD_BLAST,
                             PRE_Q_CUT_CAL_PULSER,
                             PRE_Q_CUT_WF_SATURATION,
                             PRE_Q_CUT_WF_DC_OFFSET,
                             PRE_Q_CUT_SHORT_WF,
                             PRE_Q_CUT_SAT_STRIPE,
                             PRE_Q_CUT_CPOL_SEP,
                             PRE_Q_CUT_NUM_PHISECTOR
                            };

# anita3Analysis

This repository stores code used in Jacob Gordon's thesis research.  

Code requires ROOT HEALPix, and anitaTools (OSU edit of anitaTools created by Sam Stafford, to be more precise) to run properly.  

runInterferometry.cxx analyses raw .root data files storing ANITA-III data. (data also not included)  
  The file also preforms interferometry on the data and calcualted quality cuts. 
  (Written by Sam Stafford, Lightly altered by Jacob Gordon)
  
runAnalysisStage01.cxx preforms the first round of analysis cuts, which focus on cuts that are directional or location based.
  The file takes in the .root file runInterferomerty outputs
  (Written by Sam Stafford, Edited heavily by Jacob Gordon to add cuts, new plots, new table printouts, and new functions)

runAnalysisStage02.cxx preforms the rest of the cuts.  Cuts that depend on the HEALPix binning to be preformed.
  The file can be run with or without final cuts enabled, but should be run WITHOUT final cuts enabled until optimizeLDCuts has been run.
  The file takes in the .root file runAnalysisStage01 outputs, and if finalCuts are enabled, a .txt file with bin by bun cut information
  output by optimizeLDCuts
  (Written by Sam Stafford, Edited heavily by Jacob Gordon to add cuts, new plots, new table printouts, and new features.)
  
optimizeLDCuts.cxx preforms a healpix bin by healpix bin optimiztion of the final cut paramaters, including the 
  binned analysis's primary cut, the Linear Discriminate cut.  
  The file takes in the .root file runAnalysisStage02 outputs.
  (Structure based on old code by Sam Stafford, but nearly every aspect of the file was re-written and improved by Jacob Gordon)
  ((Some large commented out sections may remain...))
  
spillover.cxx is designed to preform the complicated calcualtions required to assess the effects of using a limited data set to 
  'train' our final data cuts.  Spicifically to test the systematic uncertainty created by the assumption that the events in our
  training sample have a distribution representitive of the full data set.  Its outputs are used in optimizeLDCuts with default settings.
  A flag within optimizeLDCut can be changed so that the spillover calculations are preformed inside optimizeLDCut, but 
  they are very time consuming, so it is not recommended.
  (Written by Jacob Gordon)
  
.h files are used by multiple stages of the analysis code.
  (fitFunc.h written by Jacob Gordon, the other .h files are by Sam Stafford)

.job and .sh files are scripts written to easly and efficiently submit jobs.  They are intended for use on the Ohio Supercomputing Cluster Oakley,
  but should be easly modifiable.  They are mostly included here because I thought this would be an easy place to store them for future
  reference.  
  (Written by Jacob Gordon)

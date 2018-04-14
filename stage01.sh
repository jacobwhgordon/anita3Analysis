#we need to run stage01 for all the healpix bins...
#so, thats what this does.

#base command:
#
#./runAnalysisStage01 -CA -D/fs/scratch/PAS0174/anita/2015_05_19/sample_90/geomFilter -i0 175 439


##for i in {175..439}; 
##for i in {175..175};


##for run in {421..439}
## should be 0-7 for full segmented run numbers, or just 10 for all runs in one job.
for run in {0..7}
do
 
  runL=$(echo "176+33*$run" | bc)
  runH=$(echo "175+33+33*$run" | bc)
  if [ ${run} == 0 ]
  then
    runL=$(echo "175" | bc)
  fi
  if [ ${run} == 10 ]
  then
    runL=$(echo "175" | bc)
    runH=$(echo "439" | bc)
  fi

  for k in {0..9}
  do

    ##if [ ${k} != 10 ] && [ ${k} != 20 ] && [ ${k} != 30 ] && [ ${k} != '40' ] && [ ${k} != '50' ] && [ ${k} != '60' ] && [ ${k} != '70' ]&&[ ${k} != '80' ]&&[ ${k} != '90' ]
    ##then
      echo ${runL} ${runH}
      sed -i -e "s/\(-i\).*/\1${k} ${runL} ${runH}/" \
           -e "s/\(-N stage01\).*/\1_${run}_${k}/" stage01.job
      qsub stage01.job
      sleep .1
    ##fi

  done

done


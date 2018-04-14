#we need to run runInterferoetry for all the healpix bins...
#so, thats what this does.

##for i in {175..439}; 
##for i in {175..175};
##do

##  ./runInterferometry -PB -Ooutput ${i};

##done



##for run in {1..500}
for run in {1..500}
do

  sed -i -e "s/\(-S1\).*/\1 ${run}/" \
         -e "s/\(-N interferometry\).*/\1_sim_${run}/" runInterferometrySim.job
  eval 'qsub runInterferometrySim.job'

done


##  ./runInterferometry -PB -s1 -O/fs/scratch/PAS0174/anita/2015_05_19/sample_90/geomFilter 175

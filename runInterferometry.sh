#we need to run runInterferoetry for all the healpix bins...
#so, thats what this does.

##for i in {175..439}; 
##for i in {175..175};
##do

##  ./runInterferometry -PB -Ooutput ${i};

##done



##for run in {421..439}
for run in {175..439}
do

for k in {9..9}
do

  if [ ${run} != 255 ]
  then
    echo ${k}
    sed -i -e "s/\(-s\).*/\1${k} ${run}/" \
           -e "s/\(-N interf\).*/\1_${run}_${k}/" runInterferometry.job
    eval 'qsub runInterferometry.job'
    sleep 0.1
  fi

done

done



##  ./runInterferometry -PB -s1 -O/fs/scratch/PAS0174/anita/2015_05_19/sample_90/geomFilter 175

#we need to run runInterferoetry for all the healpix bins...
#so, thats what this does.

##for i in {175..439}; 
##for i in {175..175};
##do

##  ./runInterferometry -PB -Ooutput ${i};

##done

##mkdir /fs/scratch/PAS0174/anita/2015_05_19/sample_90/geomFilter/run${run}

##for run in {421..439}
for run in {419..419}
do

mkdir -p /fs/scratch/PAS0174/anita/2015_05_19/sample_90/geomFilter/run${run}

for k in {5..99..10}
do

  #then
    echo ${k}
    sed -i -e "s/\(-O\).*/\1\/fs\/scratch\/PAS0174\/anita\/2015_05_19\/sample_90\/geomFilter\/run${run} -S100 -s${k} ${run}/" \
           -e "s/\(-N interf\).*/\1_${run}_${k}/" runInterferometryFine.job
    eval 'qsub runInterferometryFine.job'
    sleep 0.1
  #fi

done

done



##  ./runInterferometry -PB -s1 -O/fs/scratch/PAS0174/anita/2015_05_19/sample_90/geomFilter 175
##  ./runInterferometry -Noverlap -PB --FILTER_OPTION=4 -bbaselineSampleSmooth_1_2.00.root -O/fs/scratch/PAS0174/anita/2015_05_19/sample_90/geomFilter/run255 -S50 -s9 255



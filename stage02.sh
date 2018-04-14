#we need to run stage02 for all the healpix bins...
#so, thats what this does.

##for i in {175..439}; 
##for i in {175..175};




for k in {0..9}
do

  echo ${k}
  sed -i -e "s/\(-F\).*/\1analysisOutput90\/anlysisOutput_175_439_${k}.root /" \
         -e "s/\(-N stage02\).*/\1_${k}/" stage02.job
  eval 'qsub stage02.job'
  sleep .1  

done


##file to make optimizeHP.job scrpits and submit them


for i in {0..9}
do
  theta=$(echo "-0.56*$i" | bc)
  echo $theta
  
  for j in {0..9}
  do
    phi=$(echo "0.56*$j" | bc)
    echo $phi

    sed -i -e "s/\(-PHI_HP_OFFSET=\).*/\1${phi} --THETA_HP_OFFSET=${theta}/" \
           -e "s/\(-N optimizeH\).*/\1${i}${j}/" optimizeHP.job
    
    eval 'qsub optimizeHP.job'
    sleep .1
  
  done
done

## We also wanted to check the default params, so do that here
##theta=0.0
##phi=0.0
##
##sed -i -e "s/\(-PHI_HP_OFFSET=\).*/\1${phi} --THETA_HP_OFFSET=${theta}/" \
##       -e "s/\(-N optimize\).*/\1$100/" optimizeHP.job
##
##eval 'qsub optimizeHP.job'


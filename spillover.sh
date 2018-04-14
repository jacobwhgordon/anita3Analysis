##file to make spillover.job scrpits and submit them


for i in {0..9}
do
  theta=$(echo "-0.56*$i" | bc)
  echo $theta
  
  for j in {0..9}
  do
    phi=$(echo "0.56*$j" | bc)
    echo $phi

    sed -i -e "s/\(-PHI=\).*/\1${phi} --THETA=${theta} /" \
           -e "s/\(-N spillover\).*/\1${i}${j}V/" spillover.job
    
    eval 'qsub spillover.job'
    sleep .05
  done
done

## We also wanted to check the default params, so do that here
##theta=0.56
##phi=0.0

##for k in {0..9}
##do

##  sed -i -e "s/\(-PHI=\).*/\1${phi} --THETA=${theta} ${k}/" \
##         -e "s/\(-N spillover\).*/\110V_${k}/" spillover.job
##  eval 'qsub -t 0-9 spillover.job'
##
##done

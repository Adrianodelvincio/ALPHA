
for ru in 68506 68452 68465 68475 68481 68489 68498 
do
    echo $ru
    
    cat r${ru}*timewindows.csv | awk 'BEGIN {FS=","}{print $3}' | sort | grep ^[0-9] > teststop
    cat r${ru}*timewindows.csv | awk 'BEGIN {FS=","}{print $2}' | sort | grep ^[0-9] > teststart


    nu=`wc teststart | awk '{print $1}'`
    numu=`expr $nu - 1`

    tail -$numu  teststart > statest
    head -$numu  teststop > stotest

    diff statest stotest 

done
    
#  497  grep 43784.317443 r68489_*
#  498  grep 43784.317443 r68489_*.timewindows.csv
#  499  grep 43786.339472 r68489_*.timewindows.csv


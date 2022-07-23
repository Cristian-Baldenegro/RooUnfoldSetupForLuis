FNAME=bayes_distr.txt

if [ ! -f $FNAME ]; then
    wget http://www.roma1.infn.it/~dagos/$FNAME
fi

LINENUMBERS=($(grep -n "&&&" $FNAME | cut -d ":" -f1))

sed -n "${LINENUMBERS[0]},${LINENUMBERS[1]}p" $FNAME | tail -n +2 | head -n -2 > src/bayes.for         
sed -n "${LINENUMBERS[1]},${LINENUMBERS[2]}p" $FNAME | tail -n +2 | head -n -2 > src/bayes_c.for       
sed -n "${LINENUMBERS[2]},${LINENUMBERS[3]}p" $FNAME | tail -n +2 | head -n -2 > src/bayes_test.for    
sed -n "${LINENUMBERS[3]},${LINENUMBERS[4]}p" $FNAME | tail -n +2 | head -n -2 > src/bayes_test_2d.for 





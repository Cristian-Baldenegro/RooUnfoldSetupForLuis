#!/bin/bash                                                                             

number2=0
nbins=163;
while [ $number2 -lt $nbins ]; do
/opt/exp_soft/cms/t3/t3submit -short -singleout submitSplitResponse.sh $number2 
number2=$((number2+1))
done

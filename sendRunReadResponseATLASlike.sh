#!/bin/bash                                                                                                                                                     
number2=0
nbins=20;
while [ $number2 -lt $nbins ]; do
/opt/exp_soft/cms/t3/t3submit -long -singleout submitReadResponse_ATLASlike.sh $number2 
number2=$((number2+1))
done

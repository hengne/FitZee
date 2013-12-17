#!/bin/sh

NTIMES=2

./sub_EEP_${NTIMES}.sh &> sub_EEP_${NTIMES}.log &
./sub_EEN_${NTIMES}.sh &> sub_EEN_${NTIMES}.log &

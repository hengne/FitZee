#!/bin/sh

NTIMES=2

./sub8_EEP_${NTIMES}.sh &> sub8_EEP_${NTIMES}.log &
./sub8_EEN_${NTIMES}.sh &> sub8_EEN_${NTIMES}.log &

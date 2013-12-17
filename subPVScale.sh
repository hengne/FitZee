#!/bin/sh

NTIMES="1"
REGVERSION="V8Elec"
MODE=76
METHOD=7 
RUN="A"
bsub -R "pool>20000" -q cmscaf1nd -J fitEtaScale_${REGVERSION}_${RUN}_${MODE}_${METHOD} dofitEtaScale.sh ${REGVERSION} ${RUN} ${MODE} ${METHOD} ${NTIMES}
RUN="B"
bsub -R "pool>20000" -q cmscaf1nd -J fitEtaScale_${REGVERSION}_${RUN}_${MODE}_${METHOD} dofitEtaScale.sh ${REGVERSION} ${RUN} ${MODE} ${METHOD} ${NTIMES}
RUN="C"
bsub -R "pool>20000" -q cmscaf1nd -J fitEtaScale_${REGVERSION}_${RUN}_${MODE}_${METHOD} dofitEtaScale.sh ${REGVERSION} ${RUN} ${MODE} ${METHOD} ${NTIMES}
RUN="D"
bsub -R "pool>20000" -q cmscaf1nd -J fitEtaScale_${REGVERSION}_${RUN}_${MODE}_${METHOD} dofitEtaScale.sh ${REGVERSION} ${RUN} ${MODE} ${METHOD} ${NTIMES}




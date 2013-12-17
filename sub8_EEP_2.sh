#!/bin/sh

NTIMES=2
IZ="1"
METHOD=8

IDXS=`seq 1 1 100`;

for IDX in $IDXS;
do
  date ;
  bsub -R "pool>20000" -q cmscaf1nd -J fit${NTIMES}_IZ${IZ}_IX${IDX} dofitEE.sh ${NTIMES} ${IDX} ${IZ} ${METHOD};
done

# cmscaf1nd
# 8nh

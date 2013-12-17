#!/bin/sh

NTIME="2"
tag="_method23_EvenOdd0"

outdir="out_${NTIME}${tag}"
out="calibTable_out_${NTIME}_ABCD.dat"
ref="calibTableRef_${NTIME}_ABCD.dat"

cuts="0.5 2.0 0.0 1.0"

if [ -e tmp.dat ]; then
 rm tmp.dat
fi

touch tmp.dat

# EE+
IZ="1"
IDS=`seq 1 1 100`
for IDX in $IDS;
do
  TABLEOUT="calibTable_out_IZ${IZ}_IX${IDX}_${NTIME}_ABCD.dat"
  TABLEIN="calibTable_in_IZ${IZ}_IX${IDX}_$((NTIME+1))_ABCD.dat"
  ./checkcalibtable.exe  ${outdir}/${TABLEOUT} ${cuts}   ${outdir}/${TABLEIN}
  cat  ${outdir}/${TABLEOUT} >> tmp.dat
done

#EE-
IZ="-1"
IDS=`seq 1 1 100`
for IDX in $IDS;
do
  TABLEOUT="calibTable_out_IZ${IZ}_IX${IDX}_${NTIME}_ABCD.dat"
  TABLEIN="calibTable_in_IZ${IZ}_IX${IDX}_$((NTIME+1))_ABCD.dat"
  ./checkcalibtable.exe  ${outdir}/${TABLEOUT} ${cuts}   ${outdir}/${TABLEIN}
  cat  ${outdir}/${TABLEOUT} >> tmp.dat
done



./alinecalibtable.exe tmp.dat ${outdir}/${out}

./checkcalibtable.exe ${outdir}/${out} \
               ${cuts}  \
               ${outdir}/${ref} \
            &> checkcalibtable_${NTIME}${tag}.log


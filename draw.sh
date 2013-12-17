#!/bin/sh

tag1="_EE"
list="A B C D"
ref="Ref_3_ABCD_EtaScaleV2_rmBad"
scale="etascale_V8Elec_Mode74_Method7_ABCD_2.dat"

eoslocation="root://eoscms//eos/cms/store/caf/user/heli/2013Nov06/data"

reffile="calibTable${ref}.dat"

evenodd="0"

for tag in ${list} ;
do
  ./drawMee.exe \
  ${eoslocation}/CombRegV8ZeeData16${tag}${tag1}.root \
  ${reffile} \
  drawMee_out_CombRegV8ZeeData16${tag}${tag1}_${ref}.root \
  -1 7 1.0 ${evenodd} ${scale} \
  &> log${tag}&

done 





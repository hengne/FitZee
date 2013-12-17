#!/bin/sh

tag1="_EE"
ref="Ref_3_ABCD_EtaScaleV2_rmBad"

rm drawMee_out_CombRegV8ZeeData16ABCD${tag1}_${ref}.root

hadd \
drawMee_out_CombRegV8ZeeData16ABCD${tag1}_${ref}.root \
drawMee_out_CombRegV8ZeeData16A${tag1}_${ref}.root \
drawMee_out_CombRegV8ZeeData16B${tag1}_${ref}.root \
drawMee_out_CombRegV8ZeeData16C${tag1}_${ref}.root \
drawMee_out_CombRegV8ZeeData16D${tag1}_${ref}.root 

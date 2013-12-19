#!/bin/sh

date

#n fits
NTIME=$1

#job index
IDX=$2

#iz
IZ=$3

#fit method
METHOD=$4

#mode

MODE=$5

#Gaussian Resolution
GAUSRESO="1.0"

#Hit energy fraction
HITENERGYFRACTION="0.01"

#RegVersion
REGVERSION="V7Elec"

#Eta scale file
#ETASCALEFILE="etascale_V8Elec_Mode76_Method7_A_1.dat;etascale_V8Elec_Mode76_Method7_B_1.dat;etascale_V8Elec_Mode76_Method7_C_1.dat;etascale_V8Elec_Mode76_Method7_D_1.dat"
ETASCALEFILE="etascale_V8Elec_Mode76_Method7_Odd_A_1.dat;etascale_V8Elec_Mode76_Method7_Odd_B_1.dat;etascale_V8Elec_Mode76_Method7_Odd_C_1.dat;etascale_V8Elec_Mode76_Method7_Odd_D_1.dat"

#debug
DEBUG="1"

#do odd if 1 , even if 2, both if 0
DOEVENODD="1"

#initilize working area
TOP=$PWD
SOURCE=/afs/cern.ch/user/h/heli/work/private/calib/2013Nov06/FitEBEE/fitEEEtaScaleV4
DATA=/eos/cms/store/user/heli/2013Nov06/ShervinDataEE

alias eos='/afs/cern.ch/project/eos/installation/0.2.31/bin/eos.select'
 


ROOTFILEIN=CombRegV8ZeeData_IZ${IZ}_IX${IDX}_ABCD.root
ROOTFILEOUT=fitzee_out_IZ${IZ}_IX${IDX}_${NTIME}_ABCD.root
TABLEIN=emptyCalibTable_IZ${IZ}_IX${IDX}.dat
TABLEOUT=calibTable_out_IZ${IZ}_IX${IDX}_${NTIME}_ABCD.dat
TABLEREF=emptyCalibTable_IZ${IZ}_IX${IDX}.dat
LOGFILE=fitzee_out_IZ${IZ}_IX${IDX}_${NTIME}_ABCD.log

if [ $NTIME -ge 2 ];
then
  TABLEIN="calibTable_in_IZ${IZ}_IX${IDX}_${NTIME}_ABCD.dat"
  TABLEREF="calibTableRef_$((NTIME-1))_ABCD.dat"
fi

echo "setup" 
cp ${SOURCE}/env_CMSSW_5_3_7_patch5.sh .
source env_CMSSW_5_3_7_patch5.sh
cd ${TOP}


date
#

echo "copy codes"
cd ${SOURCE}
cp calibRecord.h  fitzee.cc functions.h  variables.h voigt.h ${TOP}/
cp ${ETASCALEFILE//;/ } ${TOP}/
cd ${TOP}
date

echo "build"
g++ -o fitzee.exe fitzee.cc \
      -pthread -m64 \
      -I${ROOTSYS}/include \
      -I./  -L${ROOTSYS}/lib \
      -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d \
      -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics \
      -lMathCore -lThread -lMinuit -lMinuit2 -lpthread \
      -lTreePlayer \
      -Wl,-rpath,${ROOTSYS}/lib -lm -ldl

date

#echo "copy data"
#eos cp ${DATA}/${ROOTFILEIN} ./
#date

# use eos directly
ROOTFILEIN="root://eoscms/${DATA}/${ROOTFILEIN}" 
date

echo "copy tables"
cp ${SOURCE}/table/${TABLEIN} ./
cp ${SOURCE}/table/${TABLEREF} ./
date

echo "run"
./fitzee.exe ${MODE} \
               ${ROOTFILEIN} \
               ${ROOTFILEOUT} \
               ${TABLEIN} \
               ${TABLEOUT} \
               ${TABLEREF} \
               0.9999 \
               ${METHOD} ${GAUSRESO} \
               ${HITENERGYFRACTION} \
               ${DEBUG} \
               ${DOEVENODD} \
               ${REGVERSION} \
               ${ETASCALEFILE} \
           &> ${LOGFILE}


date

echo "copy output"
OUTDIR="${SOURCE}/out_${NTIME}_method${METHOD}_mode${MODE}_evenodd${DOEVENODD}_reg${REGVERSION}"
if [ ! -d $OUTDIR ];
then
  mkdir $OUTDIR
fi
cp ${ROOTFILEOUT} ${OUTDIR}/
cp ${TABLEOUT} ${OUTDIR}/
cp ${LOGFILE} ${OUTDIR}/


date




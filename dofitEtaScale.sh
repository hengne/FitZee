#!/bin/sh

date

# $1: REGVERSION
# $2: RUN
# $3: MODE
# $4: METHOD
# $5: NTIMES

#NTIMES
NTIMES=$5

#fit mode
MODE=$3

#fit method
METHOD=$4

#Run pieriod
RUN=$2

#Regression version
REGVERSION=$1 #"V5Elec"
#Gaussian Resolution
GAUSRESO="1.0"

#
ETASCALEREF=$6

#debug
DEBUG="1"

#initilize working area
TOP=$PWD
SOURCE=/afs/cern.ch/user/h/heli/work/private/calib/2013Nov06/FitEtaScale
DATA=/eos/cms/store/caf/user/heli/2013Nov06/data

alias eos='/afs/cern.ch/project/eos/installation/0.2.31/bin/eos.select'
 

ROOTFILEIN=CombRegV8ZeeData16${RUN}.root
ROOTFILEOUT=fitzeescale_out_${REGVERSION}_Mode${MODE}_Method${METHOD}_${RUN}_${NTIMES}.root
LOGFILE=fitzeescale_out_${REGVERSION}_Mode${MODE}_Method${METHOD}_${RUN}_${NTIMES}.log

echo "setup" 
cp ${SOURCE}/env_CMSSW_5_3_7_patch5.sh .
source env_CMSSW_5_3_7_patch5.sh
cd ${TOP}


date
#

echo "copy codes"
cd ${SOURCE}
cp calibRecord.h  fitzeescale.cc functions.h  variables.h voigt.h ${TOP}/
cp ${ETASCALEREF} ${TOP}/
cd ${TOP}
date

echo "build"
g++ -o fitzeescale.exe fitzeescale.cc \
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

date

echo "run"
./fitzeescale.exe ${MODE} \
               root://eoscms/${DATA}/${ROOTFILEIN} \
               ${ROOTFILEOUT} \
               0.9999 \
               ${METHOD} ${GAUSRESO} \
               ${REGVERSION} \
               ${DEBUG} ${ETASCALEREF}\
           &> ${LOGFILE}


date

echo "copy output"
OUTDIR="${SOURCE}/out_$REGVERSION_mode${MODE}_method${METHOD}_${NTIMES}"
if [ ! -d $OUTDIR ];
then
  mkdir $OUTDIR
fi
cp ${ROOTFILEOUT} ${OUTDIR}/
cp ${LOGFILE} ${OUTDIR}/


date




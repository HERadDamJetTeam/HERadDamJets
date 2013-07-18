#!/bin/bash

DIR=`echo $1`
YEAR=`echo $2`
ENERGY=`echo $3`
TLUMI=`echo $4`
ILUMI=`echo $5`
GBLTG=`echo $6`
PART=`echo $7`

echo ""
echo ">> `/bin/date` Submitting Condor job(s)..."
#==========================================================
#==========================================================
#==========================================================
echo ">> `/bin/date` SinglePigun: parameters in $1"

mkdir ${DIR} -p
cat ./MonoJet_DIGI-RECO_temp_split_cfg.py \
| sed -e s/YEAR/${YEAR}/ \
| sed -e s/ENERGYIN/${ENERGY}/ \
| sed -e s/LUMIDRK/${TLUMI}/ \
| sed -e s/INSTLUMI/${ILUMI}/ \
| sed -e s/MYGBLTG/${GBLTG}/ \
| sed -e s/PNUMBER/${PART}/ \
> ./${DIR}/MonoJet_DIGI-RECO_20${YEAR}_${ENERGY}_lumi${TLUMI}_${PART}_cfg.py

#add job submission commands here

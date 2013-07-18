#!/bin/bash

DIR=`echo $1`
YEAR=`echo $2`
ENERGY=`echo $3`
PART=`echo $4`

echo ""
echo ">> `/bin/date` Submitting Condor job(s)..."
#==========================================================
#==========================================================
#==========================================================
echo ">> `/bin/date` SinglePigun: parameters in $1"

mkdir ${DIR} -p
cat ./MonoJet_GEN-SIM_temp_split_cfg.py \
| sed -e s/ENERGYIN/${ENERGY}/ \
| sed -e s/YEAR/${YEAR}/ \
| sed -e s/PNUMBER/${PART}/ \
> ./${DIR}/MonoJet_GEN-SIM_20${YEAR}_${ENERGY}_${PART}_cfg.py

#add job submission commands here

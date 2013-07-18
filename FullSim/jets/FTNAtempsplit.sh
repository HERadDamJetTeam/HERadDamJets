#!/bin/bash

DIR=`echo $1`
YEAR=`echo $2`
ENERGY=`echo $3`
LUMI=`echo $4`
GBLTG=`echo $5`

mkdir ${DIR} -p
cat ./fullsimtowernoiseanalyzer_temp_split_cfg.py \
| sed -e s/YEAR/${YEAR}/ \
| sed -e s/ENERGYIN/${ENERGY}/ \
| sed -e s/LUMIDRK/${LUMI}/ \
| sed -e s/MYGBLTG/${GBLTG}/ \
> ./${DIR}/fullsimtowernoiseanalyzer_${YEAR}_${ENERGY}_lumi${LUMI}_cfg.py

cmsRun ${DIR}/fullsimtowernoiseanalyzer_${YEAR}_${ENERGY}_lumi${LUMI}_cfg.py

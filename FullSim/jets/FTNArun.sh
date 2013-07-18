#!/bin/bash

for YEAR in 19
  do
    for ENERGY in 30 #50 100
      do
		  ./FTNAtempsplit_cool.sh jet ${YEAR} ${ENERGY} 0 STAR${YEAR}_61_V1A
		  ./FTNAtempsplit_cool.sh jet ${YEAR} ${ENERGY} 100 STAR${YEAR}_61_V1A
		  ./FTNAtempsplit_cool.sh jet ${YEAR} ${ENERGY} 200 STAR${YEAR}_61_V1A
		  ./FTNAtempsplit_cool.sh jet ${YEAR} ${ENERGY} 300 STAR${YEAR}_61_V6A
		  ./FTNAtempsplit_cool.sh jet ${YEAR} ${ENERGY} 400 STAR${YEAR}_61_V6A
		  ./FTNAtempsplit_cool.sh jet ${YEAR} ${ENERGY} 500 STAR${YEAR}_61_V5A
		  ./FTNAtempsplit_cool.sh jet ${YEAR} ${ENERGY} 600 STAR${YEAR}_61_V5A
		  ./FTNAtempsplit_cool.sh jet ${YEAR} ${ENERGY} 700 STAR${YEAR}_61_V5A
	  done
  done
#!/bin/bash

for YEAR in 19
  do
    for ENERGY in 30 #50 100
      do
	    for PART in {1..20}
		  do
		    ./FStemp_digireco.sh jet ${YEAR} ${ENERGY} 0 1E34 STAR${YEAR}_61_V1A ${PART}
		    ./FStemp_digireco.sh jet ${YEAR} ${ENERGY} 100 1E34 STAR${YEAR}_61_V1A ${PART}
		    ./FStemp_digireco.sh jet ${YEAR} ${ENERGY} 200 1E34 STAR${YEAR}_61_V1A ${PART}
		    ./FStemp_digireco.sh jet ${YEAR} ${ENERGY} 300 1E34 STAR${YEAR}_61_V6A ${PART}
		    ./FStemp_digireco.sh jet ${YEAR} ${ENERGY} 400 1E34 STAR${YEAR}_61_V6A ${PART}
		    ./FStemp_digireco.sh jet ${YEAR} ${ENERGY} 500 1E34 STAR${YEAR}_61_V5A ${PART}
		    ./FStemp_digireco.sh jet ${YEAR} ${ENERGY} 600 5E34 STAR${YEAR}_61_V5A ${PART}
		    ./FStemp_digireco.sh jet ${YEAR} ${ENERGY} 700 5E34 STAR${YEAR}_61_V5A ${PART}
		  done
	  done
  done
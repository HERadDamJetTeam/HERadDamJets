#!/bin/bash

for YEAR in 17 19
  do
    for ENERGY in 30
      do
	    for PART in {1..20}
		  do
		    ./FStemp_gensim.sh jet ${YEAR} ${ENERGY} ${PART}
		  done
	  done
  done
#!/bin/bash

high_symsim () {
		Rscript high_symsim.R --ratio_ind2cond $1 \
						--scale_in_diffg $2 \
						--scale_in_nondiffg $3 \
						--ngene $4 \
						--rpt $5 \
						--algorithm $6 \
						--eta $7 \
						--adapt_engaged $8 \
						--adapt_iter $9 \
						> high_ngene_${4}-ric_${1}-sid_${2}-sin_${3}-rpt_${5}-ag_${6}-eta_${7}-ae_${8}-ai_${9}.log 2>&1
}


"$@"

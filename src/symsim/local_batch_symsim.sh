## Run symsim locally

ngene=200
nindpcond=5
brnlen=0.5
bimod=1
sigma=0.6

Rscript symsim.R \
        --ngene ${ngene} \
        --nind_per_cond ${nindpcond} \
        --brn_len ${brnlen} \
        --bimod ${bimod} \
        --sigma ${sigma} \
        > symsim_${nindpcond}subpop_${brnlen}brnlen_${bimod}bimod_0.2alpha_{ngene}gene.log 2>&1

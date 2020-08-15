mcheck_tostan_script := ${root}/src/modelcheck/modelcheck_symsim.R
mcheck_tostan_datadir := ${root}/${local_data_dir}/nobatch
mcheck_stan_outdir := ${root}/${exps}/modelcheck
mcheck_stanvi_outdir := ${mcheck_stan_outdir}/vi

mcheckseed := 1
mcheck_gwise_data := ${mcheck_tostan_datadir}/${mcheckseed}.rdump

# * generate simulation data without batch effect
${mcheck_gwise_data}:
	-mkdir -p ${mcheck_tostan_datadir}
	@Rscript ${mcheck_tostan_script}

.PHONY: simdata_mcheck
simdata_mcheck: ${mcheck_gwise_data}

# * run stan
# ** vi
mcheck_vi: ${mcheck_gwise_data}
	-mkdir -p ${mcheck_stanvi_outdir}
	${cmdstan_dir}/${vi}.sh ${stan_bin} $^ ${mcheck_stanvi_outdir}

# * summary
# ** vi



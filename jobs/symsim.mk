symsim := symsim
symsim_tostan_script := ${pipline_dir}/3-1_symsim_to_bagwiff.R
symsim_tostan_datadir := ${root}/${local_data_dir}/${symsim}/twostage_be_sysmsim/data
symsim_stan_outdir := ${root}/${exps}/${symsim}/stan
symsimseed ?= 1
symsim_gwise_data := ${symsim_to_stan_datadir}/symsim_2be_${symsimseed}.rdump

# * generate simulation data
${symsim_gwise_data}:
	@Rscript ${symsim_tostan_script} --myseed ${symsimseed}

.PHONY: simdata_symsim
simdata_symsim: ${symsim_gwise_data}

# * run stan
# ** mc
symsim_stanmc_outdir := ${symsim_stan_outdir}/${symsimseed}/mc
symsim_onemc_res := ${symsim_stanmc_outdir}/${model_version}_chain_1.csv

${symsim_onemc_res} : ${symsim_gwise_data}
	-mkdir -p ${symsim_stanmc_outdir}
	${cmdstan_dir}/${mc}.sh ${stan_bin} $< ${symsim_stanmc_outdir}

.PHONY: symsim_mc
symsim_sample : ${symsim_onemc_res}
symsim_combchanins : ${symsim_onemc_res}
	${cmdstan_dir}/${cpc}.sh ${symsim_onemc_res}

# ** vi
symsim_stanvi_outdir := ${symsim_stan_outdir}/${symsimseed}/vi

symsim_vi: ${symsim_gwise_data}
	-mkdir -p ${symsim_stanvi_outdir}
	${cmdstan_dir}/${vi}.sh ${stan_bin} $^ ${symsim_stanvi_outdir}

# * summary stan


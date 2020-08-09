# * stan compile
cmdstan_dir := ${root}/src/cmdstan
stanmodel_dir := ${root}/src/stan
stan_bin := ${curdir}/${model_version}

mc := mc
# comine stan multiple parallel chains script
cpc := cpc
vi := vi


${stan_bin}: ${stanmodel_dir}/${model_version}.stan
	cp $< ${stan_bin}.stan; \
	cd ${HOME}/softwares/cmdstan-2.23.0 ;\
	make -j4 STANCFLAGS="--include_paths=${stanmodel_dir}" $@; \
	cd - ;\
	rm ${stan_bin}.stan

.PHONY: compie
compile : ${stan_bin}

.PHONY: clean_stan
clean_stan:
	-rm ${curdir}/*.d ${curdir}/*.hpp ${curdir}/*.o
	-rm ${stan_bin}

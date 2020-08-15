root := $(abspath ..)
modelnm := mssc
model_version ?= v1-1
method ?= vi
cmdstan_dir := ${root}/src/cmdstan
stanmodel_dir := ${root}/src/stan
stan_version := 2.23.0
# stan_version := 2.24.1

curdir := $(shell pwd)
stan_bin := ${curdir}/${model_version}

# * load stan-related codes
.PHONY: compile clean_stan
${stan_bin}: ${stanmodel_dir}/${model_version}.stan
	cp $< ${stan_bin}.stan; \
	cd ${HOME}/softwares/cmdstan-${stan_version} ;\
	make -j4 STANCFLAGS="--include_paths=${stanmodel_dir}" $@; \
	cd - ;\
	rm ${stan_bin}.stan

compile: ${stan_bin}
clean_stan:
	-rm ${curdir}/*.d ${curdir}/*.hpp ${curdir}/*.o
	-rm ${stan_bin}

# * options: tasks
local_data_dir := data
exps := exps
pipline_dir := ${root}/src/pipline
mc := mc
# comine stan multiple parallel chains script
cpc := cpc
vi := vi

include uvm.mk
include pseudobulk.mk
include symsim.mk
include modelcheck.mk

.PHONY: clean
clean: clean_stan clean_uvm

#!make

# use absolute path for project root dir
# it can be defined outside
root ?= /Users/beyondpie/git-recipes/sDESeq

.PHONY: bagwiff_data rm_bagwiff_data compile sample combchains vi

local_data_dir := data
sub := UM
cell := Malignant

curdir := $(shell pwd)
mydir := ${root}/${local_data_dir}/${sub}

# * prepare stan data
pipline_dir := ${root}/src/pipline
bagwiff_r := 1-1_scRNAseq_to_bagwiff.R
bagwiff_py := 1-2_scRNAseq_to_bagwiff.py

sc_seurat := UVM_GSE139829_harmony.rds
sc_condf := genders.csv
genef := tcga_bulk_ensembl2symbol.rds

sc_rdata := UVM_scRNAseq.RData
sc_rdump := UVM_scRNAseq.rdump

${mydir}/${sc_rdata}: ${mydir}/${sc_seurat} ${mydir}/${sc_condf} ${mydir}/${genef}
	Rscript ${pipline_dir}/${bagwiff_r} --data_dir ${local_data_dir} \
                                      --sub ${sub}\
																			--celltype ${cell} \
																			--sc_file ${sc_seurat} \
																			--condf ${sc_condf} \
																			--genef ${genef} \
                                      --output ${sc_rdata}

${mydir}/${sc_rdump}: ${mydir}/${sc_rdata}
	python ${pipline_dir}/${bagwiff_py} --data_dir ${local_data_dir} \
                                      --sub ${sub} \
                                      --infl ${sc_rdata}\
                                      --outfl ${sc_rdump}


bagwiff_data: ${mydir}/${sc_rdata} ${mydir}/${sc_rdump}

rm_bagwiff_data:
	-rm ${mydir}/${sc_rdata}
	-rm ${mydir}/${sc_rdump}


# * stan process
cmdstan_dir := ${root}/src/cmdstan
stanmodel_dir := ${root}/src/stan

mc := mc
# comine stan multiple parallel chains script
cpc := cpc
vi := vi

result_dir := ${mydir}/stan
mcdir := ${result_dir}/${mc}
vidir := ${result_dir}/${vi}

MODEL ?= v1-1
model_bin := ${curdir}/${MODEL}
mybagwiffdata := ${mydir}/${sc_rdump}

onemcresult := ${mcdir}/${MODEL}_chain_1.csv

${model_bin}: ${stanmodel_dir}/${MODEL}.stan
	cd ${HOME}/softwares/cmdstan-2.23.0 ;\
	make -j4 STANCFLAGS="--include_paths=${stanmodel_dir}" $@; \
	cd -

compile : ${model_bin}

${onemcresult} : ${cmdstan_dir}/${mc}.sh ${model_bin} ${mybagwiffdata}
	-mkdir -p ${mcdir}
	$^ ${mcdir}

sample : ${onemcresult}

combchains: ${cmdstan_dir}/${cpc}.sh ${onemcresult}
	$^

vi: ${cmdstan_dir}/${vi}.sh ${model_bin} ${mybagwiffdata}
	-mkdir -p ${vidir}
	$^ ${vidir}

clean:
	-rm ${curdir}/*.d ${curdir}/*.hpp ${curdir}/*.o ${curdir}/*.log
	-rm ${model_bin}










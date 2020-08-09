#!make
# use gmake

# root ?= /home/labszu/git-recipes/mssc
.PHONY: install pybagwiff_data rm_bagwiff_data sample combchains vi \
        test_filter debulk rmbulk rdumpbagwiff \
        visum sc2cellpops pseudo_samplede pseudo_cellpops
sub := UM
cell := Malignant
mydir := ${root}/${local_data_dir}/${sub}
uvm_stan_outdir := ${root}/${exps}/${sub}/stan

# * install packages
installpackges := 0-0_install_packages.R
install: ${pipline_dir}/${installpackges}
	Rscript $<

# * prepare stan data
de_bulkRNAseq_r := 1-0_diff_bulkRNAseq_data.R
bagwiff_r := 1-1_scRNAseq_to_bagwiff.R
bagwiff_py := 1-2_scRNAseq_to_bagwiff.py

gfilteratio := 0.1
sc_seurat := UVM_GSE139829_harmony.rds
sc_condf := genders.csv

genef := tcga_bulk_gsymbol.rds
godeg := tcga_goenrich_diffexp_genes.rds
deg := tcga_diffexp_genes.rds
fpdeg := tcga_fp_diffexp_genes.rds
tndeg := tcga_tn_diffexp_genes.rds

dea := $(addprefix ${mydir}/, ${godeg} ${deg} ${fpdeg} ${tndeg})

ncell := 200
ngene := 140

sc_rdata := UVM_scRNAseq_c${ncell}g${ngene}.RData
sc_pydump := UVM_scRNAseqpy_c${ncell}g${ngene}.rdump
sc_rdump := UVM_scRNAseq_c${ncell}g${ngene}.rdump

${mydir}/${genef} ${dea} &:
	@Rscript ${pipline_dir}/${de_bulkRNAseq_r}

debulk : ${mydir}/${genef} ${dea}
rmbulk :
	-rm ${mydir}/${genef} ${dea}
	-rm *.png

${mydir}/${sc_rdump}: ${mydir}/${genef} ${dea}
	@Rscript ${pipline_dir}/${bagwiff_r} \
    --data_dir ${local_data_dir}\
		--sub ${sub}\
		--celltype ${cell}\
		--sc_file ${sc_seurat}\
		--condf ${sc_condf}\
		--genef ${genef}\
		--deg ${deg}\
    --fpdeg ${fpdeg}\
    --tndeg ${tndeg}\
		--gfilteratio ${gfilteratio} \
    --ncell ${ncell} \
    --ngene ${ngene} \
    --rdump \
		--output $(notdir $@)

rdumpbagwiff: ${mydir}/${sc_rdump}

${mydir}/${sc_rdata}: ${mydir}/${sc_seurat} ${mydir}/${sc_condf} \
                      ${mydir}/${genef} ${mydir}/${deg}
	@Rscript ${pipline_dir}/${bagwiff_r} \
    --data_dir ${local_data_dir}\
		--sub ${sub}\
		--celltype ${cell}\
		--sc_file ${sc_seurat}\
		--condf ${sc_condf}\
		--genef ${genef}\
		--deg ${deg}\
		--gfilteratio ${gfilteratio} \
    --ncell ${ncell} \
    --ngene ${ngene} \
		--output $(notdir $@)

test_filter: ${mydir}/${sc_rdata}

${mydir}/${sc_pyrump}: ${mydir}/${sc_rdata}
	@python ${pipline_dir}/${bagwiff_py} \
    --data_dir ${local_data_dir} \
		--sub ${sub} \
		--infl ${sc_rdata}\
		--outfl $(notdir $@)


pybagwiff_data: ${mydir}/${sc_rdata} ${mydir}/${sc_pydump}

rm_bagwiff_data:
	-rm ${mydir}/${sc_rdata}
	-rm ${mydir}/*.rdump

# * scRNAseq to different cell populations
sc_to_cellpops_r := 1-3_scRNAseq_to_cellpopulations.R
sc2cellpops:
	@Rscript ${pipline_dir}/${sc_to_cellpops_r}

# * run stan
uvm_mc_outdir := ${uvm_stan_outdir}/mc
uvm_vi_outdir := ${uvm_stan_outdir}/vi

mybagwiffdata := ${mydir}/${sc_rdump}
onemcresult := ${uvm_mc_outdir}/${model_version}_chain_1.csv

${onemcresult} : ${cmdstan_dir}/${mc}.sh ${stan_bin} ${mybagwiffdata}
	-mkdir -p ${uvm_mc_outdir}
	$^ ${uvm_mc_outdir}

sample : ${onemcresult}

combchains: ${cmdstan_dir}/${cpc}.sh ${onemcresult}
	$^

vi: ${cmdstan_dir}/${vi}.sh ${stan_bin} ${mybagwiffdata}
	-mkdir -p ${uvm_vi_outdir}
	$^ ${uvm_vi_outdir}

# * summary stan
stan_analyzer := 2-1_stan_vi_analyzer.R
visum:
	Rscript ${pipline_dir}/${stan_analyzer}

# * pseudo-bulk analysis
pseudeseq_analyzer := 2-2_sc_pseudobulk.R
pseudo_sampled:
	Rscript ${pipline_dir}/${pseudeseq_analyzer} --scdataf sampled_scRNAseq_summary.rds

pseudo_cellpops:
	Rscript ${pipline_dir}/${pseudeseq_analyzer} --scdataf scRNAseq_no_malignant.rds
	Rscript ${pipline_dir}/${pseudeseq_analyzer} --scdataf scRNAseq_allcell.rds

clean_uvm:
	-rm ${uvm_mc_outdir}/*.log
	-rm ${uvm_vi_outdir}/*.log

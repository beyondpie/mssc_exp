local_data_dir := data
sub := UM
cell := malignant
root := ../..

mydir := ${root}/${local_data_dir}/${sub}

pipline_dir := ${root}/src/pipline
bagwiff_r := 1-1_scRNAseq_to_bagwiff.R
bagwiff_py := 1-2_scRNAseq_to_bagwiff.py

sc_seurat := UVM_GSE139829_harmony.rds
sc_condf := genders.csv
genef := tcga_bulk_ensembl2symbol.rds

sc_rdata := UVM_scRNAseq.RData
sc_small_rdump := UVM_small_scRNAseq.rdump
sc_big_rdump := UVM_scRNAseq.rdump

${mydir}/${sc_rdata}: ${mydir}/${sc_seurat} ${mydir}/${sc_condf} ${mydir}/${genef}
	Rscript ${pipline_dir}/${bagwiff_r} --data_dir ${local_data_dir} \
                                      --sub ${sub}\
																			--celltype ${cell} \
																			--sc_file ${sc_seurat} \
																			--condf ${sc_condf} \
																			--genef ${genef} \
                                      --soutput ${sc_small_rdump} \
																			--boutput ${sc_rdata}

${mydir}/${sc_big_rdump}: ${mydir}/${sc_rdata}
	python ${pipline_dir}/${bagwiff_py} --data_dir ${local_data_dir} \
                                      --sub ${sub} \
                                      --infl ${sc_rdata}\
                                      --outfl ${sc_big_rdump}


.PHONY: bagwiff_data
bagwiff_data: ${mydir}/${sc_rdata} ${mydir}/${sc_big_rdump}


.PHONY: rm_bagwiff_data
rm_bagwiff_data:
	-rm ${mydir}/${sc_rdata}
	-rm ${mydir}/${sc_small_rdump}
	-rm ${mydir}/${sc_big_rdump}

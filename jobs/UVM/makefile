local_data_dir := data
sub := UM
cell := Malignant
root := ../..

mydir := ${root}/${local_data_dir}/${sub}

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


.PHONY: bagwiff_data
bagwiff_data: ${mydir}/${sc_rdata} ${mydir}/${sc_rdump}


.PHONY: rm_bagwiff_data
rm_bagwiff_data:
	-rm ${mydir}/${sc_rdata}
	-rm ${mydir}/${sc_rdump}

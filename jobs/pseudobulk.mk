# * pseudo-bulk analysis
uvm_analyzer := 2-2_sc_pseudobulk.R
pseudo_sampled:
	Rscript ${pipline_dir}/${uvm_analyzer} --scdataf sampled_scRNAseq_summary.rds

pseudo_cellpops:
	Rscript ${pipline_dir}/${uvm_analyzer} --scdataf scRNAseq_no_malignant.rds
	Rscript ${pipline_dir}/${uvm_analyzer} --scdataf scRNAseq_allcell.rds


SCRIPT := sbatch_symsim.sh

5ind_groupshift : ${SCRIPT}
	@sbatch $< 5 2 1 3 80 0.3 1 1

5ind_nogroupshift : ${SCRIPT}
	@sbatch $< 5 2 0 3 80 0.3 1 1

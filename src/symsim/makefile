
SCRIPT := sbatch_symsim.sh

nrpt := 3
ngene := 80
ratio_ind2cond := 0.3
scale_in_diffg := 1
scale_in_nondiffg := 1


5ind_groupshift : ${SCRIPT}
	@sbatch $< 5 2 1 ${nrpt} ${ngene} ${ratio_ind2cond} ${scale_in_diffg} ${scale_in_nondiffg}

5ind_nogroupshift : ${SCRIPT}
	@sbatch $< 5 2 0 ${nrpt} ${ngene} ${ratio_ind2cond} ${scale_in_diffg} ${scale_in_nondiffg}

3ind_groupshift : ${SCRIPT}
	@sbatch $< 3 1 1 ${nrpt} ${ngene} ${ratio_ind2cond} ${scale_in_diffg} ${scale_in_nondiffg}

3ind_nogroupshift : ${SCRIPT}
	@sbatch $< 3 1 0 ${nrpt} ${ngene} ${ratio_ind2cond} ${scale_in_diffg} ${scale_in_nondiffg}

10ind_groupshift : ${SCRIPT}
	@sbatch $< 10 3 1 ${nrpt} ${ngene} ${ratio_ind2cond} ${scale_in_diffg} ${scale_in_nondiffg}

10ind_nogroupshift : ${SCRIPT}
	@sbatch $< 10 3 0 ${nrpt} ${ngene} ${ratio_ind2cond} ${scale_in_diffg} ${scale_in_nondiffg}


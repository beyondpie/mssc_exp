
SCRIPT := sbatch_symsim.sh
symsimv2 := sbatch_indeffv2.sh

nrpt := 3
ngene := 100
ratio_ind2cond := 0.3
scale_in_diffg := 1
scale_in_nondiffg := 1

# brn_lens := 0.5 0.8 1.0 1.5
# bimods := 0 0.5 1
# sigmas := 0.4 0.6
# capt_alphas:= 0.05 0.1 0.15 0.2

## test
brn_lens := 0.5
bimods := 0
sigmas := 0.4
capt_alphas:= 0.05



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

3ind_symsim_v2: ${symsimv2}
	for brn_len in ${brn_lens}; do \
    for bimod in ${bimod}; do \
      for sigma in ${sigmas}; do \
        for alpha in ${capt_alphas}; do \
          @sbatch $< 3 $${brn_len} $${bimod} $${sigma} $${alpha} ${ngene} ;\
        done; \
      done; \
    done; \
  done

5ind_symsim_v2: ${symsimv2}
	for brn_len in ${brn_lens}; do \
    for bimod in ${bimod}; do \
      for sigma in ${sigmas}; do \
        for alpha in ${capt_alphas}; do \
          @sbatch $< 5 $${brn_len} $${bimod} $${sigma} $${alpha} ${ngene} ;\
        done; \
      done; \
    done; \
  done

10ind_symsim_v2: ${symsimv2}
	for brn_len in ${brn_lens}; do \
    for bimod in ${bimod}; do \
      for sigma in ${sigmas}; do \
        for alpha in ${capt_alphas}; do \
          @sbatch $< 10 $${brn_len} $${bimod} $${sigma} $${alpha} ${ngene} ;\
        done; \
      done; \
    done; \
  done

clean:
	-rm *.err
	-rm *.out

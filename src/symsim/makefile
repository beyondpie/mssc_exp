symsimv2 := sbatch_symsim.sh
symsimv2_test := sbatch_symsim_test.sh
nrpt := 3
ngene := 200

brn_lens := 0.5 1.0 1.5
bimods := 0 1
sigmas := 0.4 0.6
capt_alphas:= 0.05 0.1 0.15 0.2

3ind_symsim_v2: ${symsimv2}
	for brn_len in ${brn_lens}; do \
    for bimod in ${bimods}; do \
      for sigma in ${sigmas}; do \
        for alpha in ${capt_alphas}; do \
           sbatch $< 3 $${brn_len} $${bimod} $${sigma} $${alpha} ${ngene} ;\
        done; \
      done; \
    done; \
  done

5ind_symsim_v2: ${symsimv2}
	for brn_len in ${brn_lens}; do \
    for bimod in ${bimods}; do \
      for sigma in ${sigmas}; do \
        for alpha in ${capt_alphas}; do \
           sbatch $< 5 $${brn_len} $${bimod} $${sigma} $${alpha} ${ngene} ;\
        done; \
      done; \
    done; \
  done

10ind_symsim_v2: ${symsimv2}
	for brn_len in ${brn_lens}; do \
    for bimod in ${bimods}; do \
      for sigma in ${sigmas}; do \
        for alpha in ${capt_alphas}; do \
           sbatch $< 10 $${brn_len} $${bimod} $${sigma} $${alpha} ${ngene} ;\
        done; \
      done; \
    done; \
  done

test: ${symsimv2_test}
	sbatch $< 3 0.5 0 0.4 0.2 100

clean:
	-rm *.err
	-rm *.out

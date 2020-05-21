dir := $(shell pwd)
# data := my27gene18module.rdump
data := v1_2000.rdump
mc := mc.sh
cpc := cpc.sh

samdir := ./result/sample
vidir := ./result/vi

${MODEL}: ${MODEL}.stan
	cd ${HOME}/softwares/cmdstan-2.23.0 ;\
	make -j4 STANCFLAGS="--include_paths=${dir}" ${dir}/$@; \
	cd -

.PHONY: compile
compile : ${MODEL}

.PHONY: sample
sample : ${MODEL} ${data}
	./${mc} ${MODEL} ${data} ${samdir}

.PHONY: ssample
ssample : ${MODEL} ${data}
	./$(word 1, $^) method=sample  num_samples=1000 num_warmup=1000 \
              adapt random seed=12345 \
              data file=$(word 2, $^) \
              output refresh=5 file=$(word 1, $^).csv

.PHONY: vi
vi : ${MODEL} ${data}
	-mkdir -p ${vidir}
	./$(word 1, $^) variational algorithm=meanfield \
                              iter=10000 grad_samples=1 \
                              elbo_samples=100 eta=1.0 \
                              adapt engaged=1 iter=50 \
                              eval_elbo=100 output_samples=1000\
              random seed=12345 \
              data file=$(word 2, $^) \
              output file=${vidir}/$(word 1, $^).csv \
                     diagnostic_file=${vidir}/$(word 1,$^)_diag.csv \
              > ${vidir}/${word 1, $^}.log &


.PHONY: merge
merge: ${MODEL}1.csv
	./${cpc} ${MODEL}


.PHONE: clean
clean:
	-rm ${MODEL} *.d *.hpp *.o *.log


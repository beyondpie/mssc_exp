dir := $(shell pwd)
data := model_v4-20-2.rdump
mc := mc.sh
cpc := cpc.sh

${MODEL}: ${MODEL}.stan
	cd ${HOME}/softwares/cmdstan-2.23.0 ;\
	make ${dir}/$@; \
	cd -

.PHONY: compile
compile : ${MODEL}

.PHONY: sample
sample : ${MODEL} ${data}
	export STAN_NUM_THREADS=2
	./${mc} ${MODEL} ${data}

.PHONY: ssample
ssample : ${MODEL} ${data}
	./$(word 1, $^) method=sample  num_samples=1000 num_warmup=1000 \
              adapt random seed=12345 \
              data file=$(word 2, $^) \
              output refresh=5 file=$(word 1, $^).csv

.PHONY: merge
merge: ${MODEL}1.csv
	./${cpc} ${MODEL}


.PHONE: clean
clean:
	-rm ${MODEL} *.d *.hpp *.o *.log


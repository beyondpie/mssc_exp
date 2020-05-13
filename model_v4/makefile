dir := $(shell pwd)
data := model_v4.rdump

${MODEL}: ${MODEL}.stan
	cd ${HOME}/softwares/cmdstan-2.23.0 ;\
	make ${dir}/$@; \
	cd -

.PHONY: compile
	compile: ${MODEL}

# ${data}:
	# Rscript data_prepare_v3.R

.PHONY: generate
generate: ${data}

.PHONY: sample
sample : ${MODEL} ${data}
	for i in {1..1}; do \
     ./$(word 1, $^) method=$@  num_samples=6000 num_warmup=1500 \
              adapt delta=0.9 random seed=12345 \
              id=$${i} \
              data file=$(word 2, $^) \
              output file=$(word 1, $^)$${i}.csv ;\
  done;

.PHONY: ssample
ssample : ${MODEL} ${data}
	./$(word 1, $^) method=sample  num_samples=1000 num_warmup=1000 \
              adapt random seed=12345 \
              data file=$(word 2, $^) \
              output refresh=5 file=$(word 1, $^).csv


.PHONE: clean
clean:
	-rm ${MODEL} *.d *.hpp *.o *.log


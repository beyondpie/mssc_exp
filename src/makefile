dir := $(shell pwd)

model_i: model_i.stan
	cd ${HOME}/softwares/cmdstan-2.23.0 ;\
	make ${dir}/$@; \
	cd -

.PHONY: compile
compile: model_i

model_i.rdump:
	Rscript model_i.R

.PHONY: generate
generate: model_i.rdump

.PHONY: sample
sample : model_i model_i.rdump
	for i in {1..1}; do \
     ./$(word 1, $^) method=$@ adapt delta=0.85 random seed=12345 \
              id=$${i} \
              data file=$(word 2, $^) \
              output file=$(word 1, $^)$${i}.csv ;\
  done;

.PHONE: clean
clean:
	-rm model_i *.d *.hpp *.o *.csv

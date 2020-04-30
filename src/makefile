dir := $(shell pwd)

${MODEL}: ${MODEL}.stan
	cd ${HOME}/softwares/cmdstan-2.23.0 ;\
	make ${dir}/$@; \
	cd -

.PHONY: compile
compile: ${MODEL}

scHBB.rdump:
	Rscript scHBB.R

.PHONY: generate
generate: ${MODEL}.rdump

.PHONY: sample
sample : ${MODEL} scHBB.rdump
	for i in {1..1}; do \
     ./$(word 1, $^) method=$@ adapt delta=0.85 random seed=12345 \
              id=$${i} \
              data file=$(word 2, $^) \
              output file=$(word 1, $^)$${i}.csv ;\
  done;

.PHONE: clean
clean:
	-rm ${MODEL} *.d *.hpp *.o *.csv

dir := $(shell pwd)

${MODEL}: ${MODEL}.stan
	cd ${HOME}/softwares/cmdstan-2.23.0 ;\
	make ${dir}/$@; \
	cd -

.PHONY: compile
compile: ${MODEL}

# scHBB.rdump:
# Rscript scHBB.R

# .PHONY: generate
# generate: scHBB.rdump

.PHONY: sample
sample : ${MODEL} sc${GENE}.rdump
	for i in {1..1}; do \
     ./$(word 1, $^) method=$@  num_samples=6000 num_warmup=1500 thin=5 \
              adapt delta=0.9 random seed=12345 \
              id=$${i} \
              data file=$(word 2, $^) \
              output file=$(word 1, $^)$${i}_$(word 2,$^).csv ;\
  done;

.PHONY: ssample
ssample : ${MODEL} sc${GENE}.rdump
	./$(word 1, $^) method=sample  num_samples=6000 num_warmup=1500 thin=5 \
              adapt delta=0.9 random seed=12345 \
              data file=$(word 2, $^) \
              output file=$(word 1, $^)_${GENE}.csv


.PHONE: clean
clean:
	-rm ${MODEL} *.d *.hpp *.o


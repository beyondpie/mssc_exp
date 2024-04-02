dir := $(shell pwd)
data := data.rdump

${MODEL}: ${MODEL}.stan
	cd ${HOME}/softwares/cmdstan-2.23.0 ;\
	make ${dir}/$@; \
	cd -

.PHONY: compile
compile: ${MODEL}

.PHONY: sample
sample : ${MODEL} ${data}
	./${MODEL} method=sample  num_samples=1000 num_warmup=1000 save_warmup=0\
              adapt random seed=12345 \
              data file=${data} \
              output refresh=5 file=${MODEL}.csv


.PHONE: clean
clean:
	-rm test *.d *.hpp *.o *.log

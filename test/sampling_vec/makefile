dir := $(shell pwd)
MODEL := test

test: test.stan
	cd ${HOME}/softwares/cmdstan-2.23.0 ;\
	make ${dir}/$@; \
	cd -

.PHONY: compile
compile: test

.PHONY: sample
sample : ${MODEL}
	./$(word 1, $^) method=sample  num_samples=1 num_warmup=22 save_warmup=0\
              adapt random seed=12345


.PHONE: clean
clean:
	-rm test *.d *.hpp *.o *.log

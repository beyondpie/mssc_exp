dir := $(shell pwd)

test: test.stan
	cd ${HOME}/softwares/cmdstan-2.23.0 ;\
	make ${dir}/$@; \
	cd -

.PHONY: compile
mycompile: test

.PHONY: clean
clean:
	-rm test *.d *.hpp *.o *.log


.PHONY: sample
sample : test
	./$(word 1, $^) method=sample  num_samples=10 num_warmup=1000 save_warmup=0\
              adapt random seed=12345 

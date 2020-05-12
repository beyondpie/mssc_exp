dir := $(shell pwd)

test: test.stan
	cd ${HOME}/softwares/cmdstan-2.23.0 ;\
	make ${dir}/$@; \
	cd -

.PHONY: compile
	compile: test

.PHONE: clean
clean:
	-rm test *.d *.hpp *.o *.log

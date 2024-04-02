dir := $(shell pwd)
mc := mc.sh
cpc := cpc.sh
MODEL := test

${MODEL}: ${MODEL}.stan
	cd ${HOME}/softwares/cmdstan-2.23.0 ;\
	make ${dir}/$@; \
	cd -

.PHONY: compile
compile : ${MODEL}

.PHONY: sample
sample : ${MODEL}
	./${mc} ${MODEL}

.PHONY: merge
merge: ${MODEL}1.csv
	./${cpc} ${MODEL}

.PHONE: clean
clean:
	-rm ${MODEL} *.d *.hpp *.o *.log


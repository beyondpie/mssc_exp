root := $(abspath ..)
modelnm := mssc
model_version ?= v1-1

curdir := $(shell pwd)
local_data_dir := data
exps := exps
pipline_dir := ${root}/src/pipline

include stan.mk
include uvm.mk
include pseudobulk.mk
include symsim.mk

.PHONY: clean
clean: clean_stan clean_uvm

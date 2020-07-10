#!/bin/bash

model_name=$(basename ${1})

${1} variational algorithm=meanfield \
    iter=10000 grad_samples=1 \
    elbo_samples=100 eta=1.0 \
    adapt engaged=1 iter=50 \
    eval_elbo=100 output_samples=1000\
    random seed=12345 \
    data file=${2} \
    output file=${3}/${model_name}.csv \
    diagnostic_file=${3}/${model_name}_diagnose.csv \
    > ${3}/${model_name}.log &

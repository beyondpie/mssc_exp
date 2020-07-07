#!/bin/bash

${1} variational algorithm=meanfield \
    iter=10000 grad_samples=1 \
    elbo_samples=100 eta=1.0 \
    adapt engaged=1 iter=50 \
    eval_elbo=100 output_samples=1000\
    random seed=12345 \
    data file=${2} \
    output file=${3} \
    diagnostic_file=${4} \
    > ${5} &

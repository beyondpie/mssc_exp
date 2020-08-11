#!/bin/bash

# combine stan mutiple-chain results
# cpc: combine parallel chains

# select model name, basename get v1-1_chain_1 like pattern
# then only use v1-1
tmp=$(basename ${1})
model_name="${tmp/_*/}"
to_dir=$(dirname ${1})
prefix=${model_name}_chain
resultf=${to_dir}/cpc_${model_name}.csv

grep lp__ ${to_dir}/${prefix}_1.csv > ${resultf}
sed '/^[#l]/d' ${to_dir}/${model_name}*.csv >> ${resultf}

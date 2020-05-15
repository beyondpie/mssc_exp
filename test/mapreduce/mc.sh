#!/bin/bash

# multiple chains

for i in {1..1}
do
    ./${1} method=sample num_samples=1000 num_warmup=1000\
         adapt random seed=355113 \
         id=${i} data file=${2} \
         output refresh=5 file=${1}${i}.csv >${1}${i}.log &
done


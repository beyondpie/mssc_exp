#!/bin/bash

# multiple chains

for i in {1..4}
do
    ./${1} sample random seed=355113 \
         id=${i} data file=${2} \
         output file=${1}${i}.csv >${1}${i}.log &
done


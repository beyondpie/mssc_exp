#!/bin/bash

for i in 4
do
    model="model_ii_${i}"
    make MODEL=${model} compile
    for g in HBB LYZ
    do
        make MODEL=${model} GENE=${g} ssample >${model}_${g}_sample.log &
    done
done


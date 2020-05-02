#!/bin/bash

for i in 2 3 4 5
do
    model="model_i_${i}"
    make MODEL=${model} compile
    for g in HBB LYZ GPX1
    do
        make MODEL=${model} GENE=${g} ssample >${model}_${g}_sample.log &
    done
done


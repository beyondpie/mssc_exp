#!/bin/bash

for i in 3
do
    model="model_ii_${i}"
    make MODEL=${model} compile
    for g in HBA2 HBA1 CXCL8 LYZ GPX1 CST3
    do
        make MODEL=${model} GENE=${g} ssample >${model}_${g}_sample.log &
    done
done


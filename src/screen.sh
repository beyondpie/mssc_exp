#!/bin/bash

for i in 2 3 4 5
do
    for g in HBB LYZ GPX1
    do
        model="model_i_${i}"
        echo "Model: ${model} For GENE: ${g}"
        make MODEL=${model} GENE=${gene} compile
        make MODEL=${model} GENE=${gene} ssample
    done
done


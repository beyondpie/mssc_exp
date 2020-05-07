#!/bin/bash

for i in 1 2 3 4 5
do
    model="v3_${i}"
    make MODEL=${model} compile
    make MODEL=${model} ssample >${model}_sample.log &
done


#!/bin/bash

for i in 6 7
do
    model="v3_${i}"
    make MODEL=${model} compile
    make MODEL=${model} ssample >${model}_sample.log &
done


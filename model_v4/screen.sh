#!/bin/bash

for i in 1
do
    model="v4-${i}"
    make MODEL=${model} compile
    make MODEL=${model} ssample >${model}_sample.log &
done


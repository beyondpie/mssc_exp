#!/bin/bash

# combine stan mutiple-chain results
# cpc: combine parallel chains

grep lp__ ${1}1.csv > cpc_${1}.csv
sed '/^[#l]/d' ${1}*.csv >> cpc_${1}.csv

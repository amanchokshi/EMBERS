#!/bin/bash

mkdir -p ./../../outputs/beam-pointings

for i in {01..$1}
do
    wget "http://ws.mwatelescope.org/metadata/find?mintime=$2&maxtime=$3&extended=1&page=$i" -O ./../../outputs/beam-pointings/pointings_$i.json
    sleep 66
done 

#!/bin/bash

mkdir -p ./../../outputs/beam_pointings

for i in {001..$1}
do
    wget "http://ws.mwatelescope.org/metadata/find?mintime=$2&maxtime=$3&extended=1&page=$i" -O ./../../outputs/beam_pointings/pointings_$i.json
    sleep 66
done 

#!/bin/bash

mkdir ./../../outputs/beam-pointings

for i in {1..48}
do
    wget "http://ws.mwatelescope.org/metadata/find?mintime=1252195218&maxtime=1255219218&extended=1&page=$i" -O ./../../outputs/beam-pointings/pointings_$i.json
    sleep 77
done 

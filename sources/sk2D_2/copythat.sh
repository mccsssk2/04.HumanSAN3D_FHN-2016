#!/bin/bash

mkdir -p png/

for i in `seq 1 24`
do
mkdir -p png/png_${i}
cp dir_${i}/*.png png/png_${i}/
done

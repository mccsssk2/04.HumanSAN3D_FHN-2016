#!/bin/bash

for ((i = 0 ; i < 5000; i++ ))
do
	pvbatch --use-offscreen-rendering hSAN2D.py $i
done

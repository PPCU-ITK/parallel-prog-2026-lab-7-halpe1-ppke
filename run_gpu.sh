#!/bin/bash

OUTPATH="./results"

NS=(1 4 8 16)

for i in ${!NS[@]};
do
	n=(${NS[i]})
	echo "Running x${n}"
	./cfd_euler_gpu ${n} > ${OUTPATH}/gpu_x${n}.txt
done

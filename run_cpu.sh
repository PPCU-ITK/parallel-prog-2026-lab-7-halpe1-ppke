#!/bin/bash

OUTPATH="./results"

NS=(1 4 8 16)

for i in ${!NS[@]};
do
	n=(${NS[i]})
	echo "Running x${n}"
	./cfd_euler_cpu ${n} > ${OUTPATH}/cpu_x${n}.txt
done

#!/bin/bash


# necessary to create directory results/
rm -rf results/
mkdir results

#export OMP_NUM_THREADS=4


mpirun -np 4 \
      ../build/V_LQCD.x \
      -quark_ud ud.yaml \
      -quark_s s.yaml \
      -f ./in.conf_list	



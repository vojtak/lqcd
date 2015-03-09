#!/bin/bash


# necessary to create directory results/
rm -rf results/
mkdir results

export OMP_NUM_THREADS=1


mpirun -np 8 \
      ../build/V_LQCD.x \
      -quark_ud ud.yaml \
      -quark_s s.yaml \
      -f ./in.conf_list	


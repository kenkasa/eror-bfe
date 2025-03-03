#!/usr/bin/env bash

export exe=../../../packages/eror/intel/bin/ermod
export OMP_NUM_THREADS=1

mpiexec.hydra -bootstrap rsh -np 28 $exe >& ermod.log
